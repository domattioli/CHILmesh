"""W2 measurements on a loaded CHILmesh: adjacency, point location, station sampling, render.

Each function times one thing. ``python -m benchmarks.competitors.w2.measures <op> ...``
runs one op in a fresh process and prints a JSON line, so the parent can attach
wall time and peak RSS with ``measure.run_measured``.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
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

    Adjacency building is pure Python on every backend. The C++ backend only accelerates
    layerization, so layers must be built (``compute_layers=True``) for the C++ and Python
    backend timings (``adj_cpp_s``, ``adj_python_s``) to differ.
    """
    import chilmesh
    from chilmesh import CHILmesh

    t0 = time.perf_counter()
    CHILmesh(connectivity=connectivity, points=points, compute_layers=True, compute_adjacencies=True,
             build_spatial_indices=False, validate=False)
    elapsed = time.perf_counter() - t0
    info = chilmesh.backend_info()
    return {"adjacency_s": elapsed, "backend_selected": info.get("selected"), "backend_info": info}


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


def _load(path: str):
    from chilmesh import CHILmesh

    return CHILmesh.read_from_fort14(path, compute_layers=False, compute_adjacencies=False)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="benchmarks.competitors.w2.measures")
    ap.add_argument("op", choices=["parse", "adjacency", "locate", "render", "stations"])
    ap.add_argument("--mesh", required=True, help="path to a fort.14 file")
    ap.add_argument("--n-queries", type=int, default=1000)
    ap.add_argument("--stations", help="optional .npy of (n,2) query points; default: CO-OPS stations "
                                       "from --aux-dir, else element centroids")
    ap.add_argument("--aux-dir", help="directory with the fetched CO-OPS station list, points.cwl.nc and "
                                      "maxele.63.nc")
    args = ap.parse_args(argv)

    t0 = time.perf_counter()
    mesh = _load(args.mesh)
    parse_s = time.perf_counter() - t0
    if args.op == "parse":
        out = {"parse_s": parse_s, "n_nodes": mesh.n_verts, "n_elems": mesh.n_elems, "input_format": "fort.14"}
    elif args.op == "adjacency":
        out = time_adjacency(mesh.connectivity_list, mesh.points)
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


if __name__ == "__main__":
    sys.exit(main())
