"""W2 measurements on a loaded CHILmesh: adjacency, point location, station sampling, render.

Each function times one thing. ``python -m benchmarks.competitors.w2.measures <op> ...``
runs one op in a fresh process and prints a JSON line, so the parent can attach
wall time and peak RSS with ``measure.run_measured``.
"""
from __future__ import annotations

import argparse
import importlib.metadata
import json
import hashlib
import math
import sys
import tempfile
import time
from pathlib import Path
from typing import Callable, Optional

import numpy as np


def sample_centroids(points, connectivity, n: int, seed: int = 0) -> np.ndarray:
    """``n`` element centroids (with replacement, fixed seed): inside-mesh query points."""
    from benchmarks.competitors.metrics.quality import split_elements

    pts = np.asarray(points, dtype=float)[:, :2]
    tris, quads = split_elements(connectivity)
    cents = np.concatenate([pts[tris].mean(axis=1)] + ([pts[quads].mean(axis=1)] if len(quads) else []))
    idx = np.random.default_rng(seed).integers(0, len(cents), size=n)
    return cents[idx]


def time_point_location(mesh, query_points, *, clock: Callable[[], float] = time.perf_counter) -> dict:
    """Time spatial-tree construction apart from the queries (plan R6, no double counting).

    ``tree_build_s`` re-runs ``mesh._build_spatial_indices()`` (the cKDTree the
    constructor builds). ``query_total_s`` covers only the ``find_element`` calls.
    Points outside the mesh return -1 and are counted in ``n_outside``.
    """
    q = np.asarray(query_points, dtype=float)
    t0 = clock()
    mesh._build_spatial_indices()
    tree_build_s = clock() - t0
    outside = 0
    t0 = clock()
    for p in q:
        if mesh.find_element(p) < 0:
            outside += 1
    query_total_s = clock() - t0
    n = len(q)
    return {"n_queries": n, "n_inside": n - outside, "n_outside": outside,
            "tree_build_s": tree_build_s, "query_total_s": query_total_s,
            "query_per_s": query_total_s / n if n else None}


def time_adjacency(connectivity, points) -> dict:
    """Time mesh construction with adjacencies and layers; report the active backend (``backend_info``).

    In CHILmesh 1.6.0 neither adjacency nor peel is accelerated by the C++ backend;
    both backend timings measure the same Python operations.
    """
    import chilmesh
    from chilmesh import CHILmesh

    t0 = time.perf_counter()
    CHILmesh(connectivity=connectivity, points=points, compute_layers=True, compute_adjacencies=True,
             build_spatial_indices=False, validate=False)
    elapsed = time.perf_counter() - t0
    info = chilmesh.backend_info()
    return {"adjacency_s": elapsed, "backend_selected": info.get("selected"), "backend_info": info}


def full_init(mesh) -> dict:
    """Time full construction from the parsed mesh and report the backend."""
    import chilmesh
    from chilmesh import CHILmesh

    t0 = time.perf_counter()
    CHILmesh(connectivity=mesh.connectivity_list, points=mesh.points,
             compute_layers=True, compute_adjacencies=True)
    return {"full_init_s": time.perf_counter() - t0, "backend_info": chilmesh.backend_info()}


def write_roundtrip(path, mesh) -> dict:
    """Time writing and compare the original and written raw fort.14 content."""
    from chilmesh import read_fort14_raw

    original = read_fort14_raw(path)
    with tempfile.TemporaryDirectory() as directory:
        output = str(Path(directory) / "roundtrip.14")
        t0 = time.perf_counter()
        mesh.write_to_fort14(output)
        elapsed = time.perf_counter() - t0
        written = read_fort14_raw(output)
    differences = []
    if original.node_ids != written.node_ids or original.coords != written.coords:
        differences.append("nodes differ")
    if original.elem_ids != written.elem_ids or original.elements != written.elements:
        differences.append("elements differ")
    if (original.open_boundaries != written.open_boundaries or
            original.flow_boundaries != written.flow_boundaries):
        differences.append("boundaries differ")
    return {"write_fort14_s": elapsed, "roundtrip_ok": not differences, "reason": "; ".join(differences)}


def install_mb() -> float:
    """Size of existing files listed by the installed chilmesh distribution, in MiB."""
    dist = importlib.metadata.distribution("chilmesh")
    if not dist.files:
        raise ValueError("installed chilmesh distribution has no file list")
    return sum(p.stat().st_size for f in dist.files
               if (p := Path(dist.locate_file(f))).is_file()) / (1024 * 1024)


def station_sample_max_diff(mesh, values, stations, reference) -> dict:
    """Sample nodal ``values`` at ``stations`` and return the max absolute difference from ``reference``.

    Each station is located with ``mesh.find_element``. The sampled value is the inverse-distance
    weighted mean of the nodal values at the containing element's vertices. ``reference`` holds the
    ``points.cwl.nc`` values, one per station.

    Stations outside the mesh are counted in ``n_outside`` and skipped, never silently dropped.
    """
    values = np.asarray(values, dtype=float)
    reference = np.asarray(reference, dtype=float)
    diffs, outside = [], 0
    for p, ref in zip(np.asarray(stations, dtype=float), reference):
        eid = mesh.find_element(p)
        if eid < 0:
            outside += 1
            continue
        row = mesh.connectivity_list[eid]
        verts = list(dict.fromkeys(int(v) for v in row))
        xy = mesh.points[verts, :2]
        w = 1.0 / np.maximum(np.linalg.norm(xy - p, axis=1), 1e-12)
        diffs.append(abs(float(np.sum(w * values[verts]) / np.sum(w)) - float(ref)))
    return {"max_abs_diff": max(diffs) if diffs else None, "n_sampled": len(diffs), "n_outside": outside,
            "method": "inverse-distance over the containing element's vertices"}


def station_query_points(stations_xy, n: int) -> np.ndarray:
    """``n`` query points from station coordinates, cycled when there are fewer than ``n`` stations."""
    xy = np.asarray(stations_xy, dtype=float)[:, :2]
    return xy[np.arange(n) % len(xy)]


def locate_queries(mesh, n: int, *, aux_dir=None, seed: int = 0):
    """Query points for point location and their source label.

    CO-OPS station coordinates (from ``noaa.load_station_coords``) when the station list is in
    ``aux_dir``; otherwise element centroids, which are always inside the mesh.
    """
    from benchmarks.competitors import noaa

    stations_path, _ = noaa.aux_paths(aux_dir)
    xy = noaa.load_station_coords(stations_path) if stations_path else None
    if xy is not None and len(xy):
        return station_query_points(xy, n), "coops-stations"
    return sample_centroids(mesh.points, mesh.connectivity_list, n, seed=seed), "element-centroids"


def load_maxele_field(mesh, aux_dir, *, maxele_reader=None):
    """Nodal max-elevation field from ``aux_dir``, or None when the maxele file is absent.

    Raises ``noaa.MaxeleError`` when the file lacks ``zeta_max`` or its length differs from the node count.
    Mesh nodal z (bathymetry) is never a fallback.
    """
    from benchmarks.competitors import noaa

    path = noaa.maxele_path(aux_dir)
    if path is None:
        return None
    vals = noaa.read_maxele(path, reader=maxele_reader)
    if len(vals) != len(mesh.points):
        raise noaa.MaxeleError(f"{noaa.MAXELE_FILE} has {len(vals)} values but the mesh has "
                               f"{len(mesh.points)} nodes")
    return vals


def render_values(mesh, aux_dir, *, maxele_reader=None):
    """Field to render: the max-elevation field, or None when it is not available."""
    return load_maxele_field(mesh, aux_dir, maxele_reader=maxele_reader)


def run_station_check(mesh, aux_dir, *, reader=None, maxele_reader=None) -> dict:
    """Sample the max-elevation field at the ``points.cwl.nc`` stations; compare with its peak zeta (T023).

    The reference per station is the max over time of ``zeta`` in ``points.cwl.nc``. Returns
    ``status: not-run`` with "NOAA inputs not fetched" unless both the CO-OPS station list and
    ``points.cwl.nc`` are in ``aux_dir``, and with "NOAA inputs not fetched (maxele)" when
    ``maxele.63.nc`` is absent. A missing ``zeta_max`` variable gives ``status: failed`` and a reason.
    """
    from benchmarks.competitors import noaa

    stations_path, points_path = noaa.aux_paths(aux_dir)
    if stations_path is None or points_path is None:
        return {"status": "not-run", "reason": noaa.NOAA_INPUTS_REASON}
    try:
        field = load_maxele_field(mesh, aux_dir, maxele_reader=maxele_reader)
    except noaa.MaxeleError as exc:
        return {"status": "failed", "reason": str(exc)}
    if field is None:
        return {"status": "not-run", "reason": noaa.MAXELE_REASON}
    xy, ref = noaa.read_points_cwl(points_path, reader=reader)
    return {"status": "ok", **station_sample_max_diff(mesh, field, xy, ref)}


def run_render_check(mesh, aux_dir, *, maxele_reader=None, out_png: Optional[str] = None) -> dict:
    """Time rendering of the max-elevation field; ``not-run`` when it is absent, ``failed`` when unreadable."""
    from benchmarks.competitors import noaa

    try:
        field = render_values(mesh, aux_dir, maxele_reader=maxele_reader)
    except noaa.MaxeleError as exc:
        return {"status": "failed", "reason": str(exc), "render_s": None}
    if field is None:
        return {"status": "not-run", "reason": noaa.MAXELE_REASON, "render_s": None}
    return {"status": "ok", **time_render(mesh, field, out_png=out_png)}


def time_render(mesh, values, *, out_png: Optional[str] = None) -> dict:
    """Time ``plot_field`` with a non-interactive backend."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from chilmesh.chilplotting import plot_field

    t0 = time.perf_counter()
    fig, ax = plt.subplots()
    plot_field(mesh, values, ax=ax)
    if out_png:
        fig.savefig(out_png, dpi=100)
    else:
        fig.canvas.draw()
    elapsed = time.perf_counter() - t0
    plt.close(fig)
    return {"render_s": elapsed}


def _load(path: str, *, mesh_class=None):
    if mesh_class is None:
        from chilmesh import CHILmesh as mesh_class

    return mesh_class.read_from_fort14(path, compute_layers=False, compute_adjacencies=False)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="benchmarks.competitors.w2.measures")
    ap.add_argument("op", choices=["parse", "adjacency", "locate", "render", "stations", "write", "full-init", "adjacency-cpp", "adjacency-python", "locate-tree-build", "locate-query"])
    ap.add_argument("--mesh", required=True, help="path to a fort.14 file")
    ap.add_argument("--n-queries", type=int, default=1000)
    ap.add_argument("--stations", help="optional .npy of (n,2) query points; default: CO-OPS stations "
                                       "from --aux-dir, else element centroids")
    ap.add_argument("--aux-dir", help="directory with the fetched CO-OPS station list, points.cwl.nc and "
                                      "maxele.63.nc")
    ap.add_argument("--fair", action="store_true")
    ap.add_argument("--output", help="persistent fair operation output artifact")
    args = ap.parse_args(argv)
    if args.fair:
        return fair_main(args)

    # Resolve the lazy class import before timing, as competitor adapters do.
    from chilmesh import CHILmesh

    t0 = time.perf_counter()
    mesh = _load(args.mesh, mesh_class=CHILmesh)
    parse_s = time.perf_counter() - t0
    if args.op == "parse":
        out = {"parse_s": parse_s, "n_nodes": mesh.n_verts, "n_elems": mesh.n_elems, "input_format": "fort.14"}
    elif args.op == "adjacency":
        out = time_adjacency(mesh.connectivity_list, mesh.points)
    elif args.op == "write":
        out = write_roundtrip(args.mesh, mesh)
    elif args.op == "full-init":
        out = full_init(mesh)
    elif args.op == "locate":
        if args.stations:
            q, source = np.load(args.stations), "npy-file"
        else:
            q, source = locate_queries(mesh, args.n_queries, aux_dir=args.aux_dir)
        out = {**time_point_location(mesh, q), "query_source": source}
    elif args.op == "stations":
        out = run_station_check(mesh, args.aux_dir)
    else:
        out = run_render_check(mesh, args.aux_dir)
    print(json.dumps(out))
    return 0



def chilmesh_output(mesh):
    """Expose actual records in row order; never invent original IDs here.

    The current reader discards file IDs. The child reports its rows and declares
    record_order='row-order'. Only the controller may map rows to IDs, and only when
    the independent reference proves the file IDs are 1..N in file order
    (verification.verify_output, operator ruling 2026-11-12).
    """
    ids = getattr(mesh, 'node_ids', None)
    elem_ids = getattr(mesh, 'element_ids', None)
    out = dict(n_nodes=int(mesh.n_verts), n_elems=int(mesh.n_elems))
    rows = [list(map(int, row)) for row in mesh.connectivity_list]
    # CHILmesh's mixed triangles use [a,b,c,a], never sorted vertices.
    rows = [row[:3] if len(row) == 4 and row[-1] == row[0] else row for row in rows]
    if ids is None or elem_ids is None:
        return dict(out, coordinates=mesh.points.tolist(), connectivity=rows, connectivity_indexing='row-offset',
                    record_order='row-order', node_ids=None, element_ids=None, mapping_provenance=None)
    out.update(node_ids=list(map(int, ids)), element_ids=list(map(int, elem_ids)),
               coordinates=mesh.points.tolist(), connectivity=rows,
               connectivity_indexing='row-offset', index_to_node_id=list(map(int, ids)),
               mapping_provenance={'kind': 'explicit materialized node_ids/element_ids', 'triangle_padding': 'repeated first vertex'})
    return out


def prepare_fair_operation(operation, mesh_path, *, mesh_class=None, n_queries=1000,
                           aux_dir=None, stations=None, output_path=None):
    """Resolve imports and declared non-parse prerequisites before READY.

    Each call prepares one process invocation. Parse never occurs in preparation
    for the parse operation. locate-query has its own prebuilt tree; tree-build
    performs no queries. Existing legacy helpers remain unchanged.
    """
    from benchmarks.competitors import noaa
    if mesh_class is None:
        from chilmesh import CHILmesh as mesh_class
    load = lambda: _load(mesh_path, mesh_class=mesh_class)
    state = {'mesh': None, 'operation': operation}
    if operation == 'parse':
        def work():
            state['mesh'] = load()
            return state['mesh']
    else:
        mesh = state['mesh'] = load()
        if operation in ('adjacency-cpp', 'adjacency-python', 'full-init'):
            def work():
                return mesh_class(connectivity=mesh.connectivity_list, points=mesh.points,
                                  compute_layers=True, compute_adjacencies=True,
                                  build_spatial_indices=operation == 'full-init', validate=False)
        elif operation in ('locate-tree-build', 'locate-query'):
            if operation == 'locate-tree-build':
                # The public reader builds both trees unconditionally. Drop
                # them before READY; this cell measures new tree construction.
                for attribute in ('_vertex_tree', '_centroid_tree'):
                    if hasattr(mesh, attribute): delattr(mesh, attribute)
                def work():
                    mesh._build_spatial_indices()
                    return {'tree_built': getattr(mesh, '_centroid_tree', None) is not None}
            else:
                q, source = (np.load(stations), 'npy-file') if stations else locate_queries(mesh, n_queries, aux_dir=aux_dir)
                if len(q) == 0: raise ValueError('empty query contract')
                state['query_source'] = source
                state['query_sha256'] = hashlib.sha256(q.tobytes()).hexdigest()
                mesh._build_spatial_indices()
                def work():
                    return {'element_indices': [int(mesh.find_element(p)) for p in q], 'n_queries': len(q)}
        elif operation == 'write':
            if output_path is None: raise ValueError('write requires persistent output path')
            if Path(output_path).exists(): raise FileExistsError(output_path)
            def work():
                result = mesh.write_to_fort14(output_path)
                return {'written': result is not False, 'path': str(output_path)}
        elif operation in ('render', 'stations'):
            field = load_maxele_field(mesh, aux_dir)
            if field is None: raise ValueError(noaa.MAXELE_REASON)
            state['field_sha256'] = hashlib.sha256(field.tobytes()).hexdigest()
            if operation == 'stations':
                station_path, points_path = noaa.aux_paths(aux_dir)
                if station_path is None or points_path is None: raise ValueError(noaa.NOAA_INPUTS_REASON)
                xy, reference = noaa.read_points_cwl(points_path)
                def work(): return station_sample_max_diff(mesh, field, xy, reference)
            else:
                import matplotlib
                matplotlib.use('Agg')
                import matplotlib.pyplot as plt
                from chilmesh.chilplotting import plot_field
                if output_path is None: raise ValueError('render requires persistent output path')
                if Path(output_path).exists(): raise FileExistsError(output_path)
                def work():
                    fig, ax = plt.subplots()
                    try:
                        plot_field(mesh, field, ax=ax)
                        fig.savefig(output_path, dpi=100)
                    finally:
                        plt.close(fig)
                    return {'path': str(output_path)}
        else:
            raise ValueError('unsupported fair operation')
    state['operation'] = work
    state['materialize'] = lambda output: output
    return state


def verify_fair_operation(operation, prepared, output):
    """Check completion after END; file artifacts stay available for the parent."""
    evidence = {'verification_state': 'incomplete', 'reason': 'operation output contract incomplete'}
    if operation == 'parse':
        return {'verification_state': 'verified', 'reason': 'materialized parse; parent checks independent canonical content'}
    if operation in ('adjacency-cpp', 'adjacency-python', 'full-init'):
        if getattr(output, 'adjacencies', None) is not None and getattr(output, 'layers', None) is not None:
            evidence = {'verification_state': 'verified', 'reason': '', 'n_nodes': int(output.n_verts), 'n_elems': int(output.n_elems)}
    elif operation == 'locate-tree-build':
        evidence = dict(output, verification_state='verified' if output.get('tree_built') else 'failed', reason='tree materialization check')
    elif operation == 'locate-query':
        indices = output.get('element_indices', [])
        valid = len(indices) == output.get('n_queries') and len(indices) > 0 and all(-1 <= i < prepared['mesh'].n_elems for i in indices)
        evidence = dict(output, verification_state='verified' if valid else 'failed', reason='query completion/range check', query_sha256=prepared['query_sha256'])
    elif operation == 'stations':
        valid = output.get('n_sampled', 0) > 0 and output.get('max_abs_diff') is not None and math.isfinite(output['max_abs_diff'])
        evidence = dict(output, verification_state='verified' if valid else 'failed', reason='nonempty finite station comparison')
    elif operation in ('write', 'render'):
        path = Path(output['path'])
        if path.is_file() and path.stat().st_size > 0:
            with path.open('rb') as stream: digest = hashlib.file_digest(stream, 'sha256').hexdigest()
            evidence.update(path=str(path), size_bytes=path.stat().st_size, sha256=digest)
            if operation == 'render':
                with path.open('rb') as stream: valid = stream.read(8) == b'\x89PNG\r\n\x1a\n'
                evidence.update(verification_state='verified' if valid else 'failed', reason='PNG signature/completed render')
            else:
                # Independent content verification of the written output; no
                # source reparse or CHILmesh roundtrip can certify the reference.
                evidence['reason'] = 'written artifact retained; independent parent roundtrip reference required'
    return evidence


def fair_main(args):
    import chilmesh
    import chilmesh.mesh_topology  # noqa: F401  loaded lazily by adjacency; import cost stays outside the timed region
    from benchmarks.competitors.protocol import protocol_from_environment
    version = importlib.metadata.version('chilmesh')
    backend = chilmesh.backend_info()
    implementations = {'parse': 'python fort.14 parser plus scipy.spatial.cKDTree',
                       'locate-tree-build': 'scipy.spatial.cKDTree',
                       'locate-query': 'python query loop plus scipy.spatial.cKDTree',
                       'adjacency-cpp': 'python adjacency and layers',
                       'adjacency-python': 'python adjacency and layers',
                       'full-init': 'python adjacency/layers plus scipy.spatial.cKDTree',
                       'render': 'matplotlib', 'stations': 'inverse-distance station sampling', 'write': 'python fort.14 writer'}
    prepared = prepare_fair_operation(args.op, args.mesh, n_queries=args.n_queries, aux_dir=args.aux_dir,
                                      stations=args.stations, output_path=args.output)
    protocol = protocol_from_environment()
    if protocol is None: raise ValueError('fair worker requires inherited control/result channels')
    protocol.ready()
    def verify(output):
        mesh = output if args.op in ('parse', 'full-init', 'adjacency-cpp', 'adjacency-python') else prepared['mesh']
        return dict(tool='chilmesh', version=version, materialization=True, output=chilmesh_output(mesh),
                    actual_backend={'selected': backend.get('selected'), 'operation_implementation': implementations[args.op], 'backend_info': backend},
                    operation_output=verify_fair_operation(args.op, prepared, output))
    protocol.run(prepared['operation'], materialize=prepared['materialize'], verify=verify)
    return 0

if __name__ == "__main__":
    sys.exit(main())
