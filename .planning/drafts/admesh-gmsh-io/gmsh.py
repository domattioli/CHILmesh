"""Gmsh ASCII .msh I/O for admesh (feature 008).

Defined by specs/008-gmsh-io-integration/ — contracts/python-api.md
(signatures), data-model.md (entity shapes and algorithms),
research.md (grammar decisions R1-R10).

Hand-parsed ASCII only — the gmsh PyPI package is NOT a dependency (FR-005).
"""
from __future__ import annotations

import os
import warnings
from typing import Iterator

import numpy as np

from admesh.api import BoundarySegment, Mesh
from admesh.boundary_types import BoundaryType

__all__ = ["GmshParseError", "PHYSICAL_GROUP_MAP", "read_msh", "write_msh"]

PHYSICAL_GROUP_MAP: dict[str, BoundaryType] = {
    "open": BoundaryType.OPEN,
    "mainland": BoundaryType.MAINLAND,
    "island": BoundaryType.ISLAND,
    "mainland_flux": BoundaryType.MAINLAND_FLUX,
}

_BOUNDARY_TYPE_TO_NAME: dict[int, str] = {
    int(v): k for k, v in PHYSICAL_GROUP_MAP.items()
}


class GmshParseError(ValueError):
    def __init__(self, detail: str, *, line_number: int | None = None) -> None:
        self.detail = detail
        self.line_number = line_number
        loc = f" (line {line_number})" if line_number is not None else ""
        super().__init__(f"Gmsh parse error{loc}: {detail}")


class _Cursor:
    __slots__ = ("_lines", "_pos")
    def __init__(self, lines: list[str]) -> None:
        self._lines = lines
        self._pos = 0

    @property
    def line_no(self) -> int:
        return self._pos

    def next(self) -> str | None:
        if self._pos < len(self._lines):
            line = self._lines[self._pos]
            self._pos += 1
            return line
        return None

    def require(self, expected: str = "") -> str:
        line = self.next()
        if line is None:
            raise GmshParseError(f"expected {expected!r}, got <EOF>", line_number=self._pos)
        return line

    def fail(self, detail: str) -> GmshParseError:
        return GmshParseError(detail, line_number=self._pos)


def _skip_block(cursor: _Cursor, end_tag: str) -> None:
    while True:
        line = cursor.next()
        if line is None:
            return
        if line.strip() == end_tag:
            return


# --- Internal data ---
class _PhysGroup:
    __slots__ = ("dim", "tag", "name")
    def __init__(self, dim, tag, name):
        self.dim, self.tag, self.name = dim, tag, name


class _Elem:
    __slots__ = ("etype", "node_ids", "phys_tag")
    def __init__(self, etype, node_ids, phys_tag):
        self.etype, self.node_ids, self.phys_tag = etype, node_ids, phys_tag


class _NodeDataView:
    __slots__ = ("name", "values")
    def __init__(self, name, values):
        self.name, self.values = name, values


class _GmshFile:
    __slots__ = ("version", "nodes", "elements", "physical_groups", "node_data_views")
    def __init__(self, version, nodes, elements, physical_groups, node_data_views):
        self.version = version
        self.nodes = nodes
        self.elements = elements
        self.physical_groups = physical_groups
        self.node_data_views = node_data_views


def _parse_physical_names(cursor: _Cursor) -> dict[int, _PhysGroup]:
    groups: dict[int, _PhysGroup] = {}
    n = int(cursor.require("count").strip())
    for _ in range(n):
        raw = cursor.require("PhysicalNames entry").strip()
        parts = raw.split()
        if len(parts) < 3:
            raise cursor.fail(f"malformed PhysicalNames: {raw!r}")
        dim = int(parts[0])
        tag_id = int(parts[1])
        if '"' in raw:
            name = raw.split('"', 1)[1].rsplit('"', 1)[0]
        else:
            name = parts[2]
        groups[tag_id] = _PhysGroup(dim, tag_id, name)
    cursor.require("$EndPhysicalNames")
    return groups


def _parse_node_data_block(cursor: _Cursor, n_total_nodes: int) -> _NodeDataView | None:
    n_str = int(cursor.require("NodeData numStringTags").strip())
    name = ""
    for i in range(n_str):
        raw = cursor.require("NodeData string tag").strip()
        if i == 0:
            name = raw.strip('"')
    n_real = int(cursor.require("NodeData numRealTags").strip())
    for _ in range(n_real):
        cursor.require("NodeData real tag")
    n_int = int(cursor.require("NodeData numIntegerTags").strip())
    int_tags = [int(cursor.require("NodeData int tag").strip()) for _ in range(n_int)]
    n_vals = int_tags[2] if len(int_tags) > 2 else n_total_nodes
    is_bathy = name.lower() == "bathymetry"
    values_dict: dict[int, float] = {}
    for _ in range(n_vals):
        raw = cursor.require("NodeData value").strip().split()
        tag_id = int(raw[0])
        val = float(raw[1])
        if is_bathy:
            values_dict[tag_id] = val
    cursor.require("$EndNodeData")
    if not is_bathy or not values_dict:
        return None
    return _NodeDataView(name=name, values=values_dict)


def _parse_v22(cursor: _Cursor) -> _GmshFile:
    nodes: dict[int, tuple] = {}
    elements: list[_Elem] = []
    physical_groups: dict[int, _PhysGroup] = {}
    node_data_views: list[_NodeDataView] = []
    while True:
        line = cursor.next()
        if line is None:
            break
        tag = line.strip()
        if not tag or tag.startswith("//"):
            continue
        if tag == "$PhysicalNames":
            physical_groups = _parse_physical_names(cursor)
        elif tag == "$Nodes":
            n = int(cursor.require("Nodes count").strip())
            for _ in range(n):
                raw = cursor.require("node").strip().split()
                if len(raw) < 4:
                    raise cursor.fail(f"malformed node: {raw}")
                tid = int(raw[0])
                nodes[tid] = (float(raw[1]), float(raw[2]), float(raw[3]))
            cursor.require("$EndNodes")
        elif tag == "$Elements":
            n = int(cursor.require("Elements count").strip())
            for _ in range(n):
                raw = cursor.require("element").strip().split()
                if len(raw) < 4:
                    raise cursor.fail(f"malformed element: {raw}")
                etype = int(raw[1])
                ntags = int(raw[2])
                phys_tag = int(raw[3]) if ntags >= 1 else None
                npt = {1: 2, 2: 3, 15: 1}.get(etype, -1)
                if etype not in (1, 2, 15):
                    raise GmshParseError(
                        f"higher-order or unsupported element type {etype}",
                        line_number=cursor.line_no,
                    )
                if etype == 15:
                    continue
                ns = 3 + ntags
                nids = tuple(int(raw[ns + k]) for k in range(npt))
                elements.append(_Elem(etype, nids, phys_tag))
            cursor.require("$EndElements")
        elif tag == "$NodeData":
            v = _parse_node_data_block(cursor, len(nodes))
            if v is not None:
                node_data_views.append(v)
        elif tag.startswith("$End"):
            pass
        else:
            _skip_block(cursor, "$End" + tag[1:])
    return _GmshFile("2.2", nodes, elements, physical_groups, node_data_views)


def _parse_v41(cursor: _Cursor) -> _GmshFile:
    nodes: dict[int, tuple] = {}
    elements: list[_Elem] = []
    physical_groups: dict[int, _PhysGroup] = {}
    entity_phys: dict[tuple, list] = {}
    node_data_views: list[_NodeDataView] = []
    while True:
        line = cursor.next()
        if line is None:
            break
        tag = line.strip()
        if not tag or tag.startswith("//"):
            continue
        if tag == "$PhysicalNames":
            physical_groups = _parse_physical_names(cursor)
        elif tag == "$Entities":
            hdr = cursor.require("Entities header").strip().split()
            n_pts, n_curves, n_surfs, n_vols = int(hdr[0]), int(hdr[1]), int(hdr[2]), int(hdr[3])
            for dim, count in [(0, n_pts), (1, n_curves), (2, n_surfs), (3, n_vols)]:
                for _ in range(count):
                    raw = cursor.require(f"entity dim={dim}").strip().split()
                    etag = int(raw[0])
                    # bbox = 6 floats (indices 1-6), then numPhysTags at index 7
                    n_phys = int(raw[7])
                    phys_tags = [int(raw[8 + k]) for k in range(n_phys)]
                    entity_phys[(dim, etag)] = phys_tags
            cursor.require("$EndEntities")
        elif tag == "$Nodes":
            hdr = cursor.require("Nodes header").strip().split()
            n_blocks = int(hdr[0])
            for _ in range(n_blocks):
                bh = cursor.require("node block header").strip().split()
                n_in_block = int(bh[3])
                tags = [int(cursor.require("node tag").strip()) for _ in range(n_in_block)]
                for tid in tags:
                    coord = cursor.require("node coords").strip().split()
                    nodes[tid] = (float(coord[0]), float(coord[1]), float(coord[2]))
            cursor.require("$EndNodes")
        elif tag == "$Elements":
            hdr = cursor.require("Elements header").strip().split()
            n_blocks = int(hdr[0])
            for _ in range(n_blocks):
                bh = cursor.require("element block header").strip().split()
                edim, etag_id, etype, n_in = int(bh[0]), int(bh[1]), int(bh[2]), int(bh[3])
                phys_list = entity_phys.get((edim, etag_id), [])
                phys_tag = phys_list[0] if phys_list else None
                if len(phys_list) > 1:
                    warnings.warn(
                        f"entity ({edim},{etag_id}) in multiple physical groups; using {phys_tag}",
                        UserWarning, stacklevel=4,
                    )
                npt = {1: 2, 2: 3, 15: 1}.get(etype, -1)
                if etype not in (1, 2, 15):
                    raise GmshParseError(
                        f"higher-order or unsupported element type {etype}",
                        line_number=cursor.line_no,
                    )
                for _ in range(n_in):
                    raw = cursor.require("element").strip().split()
                    if etype == 15:
                        continue
                    nids = tuple(int(raw[1 + k]) for k in range(npt))
                    elements.append(_Elem(etype, nids, phys_tag))
            cursor.require("$EndElements")
        elif tag == "$NodeData":
            v = _parse_node_data_block(cursor, len(nodes))
            if v is not None:
                node_data_views.append(v)
        elif tag.startswith("$End"):
            pass
        else:
            _skip_block(cursor, "$End" + tag[1:])
    return _GmshFile("4.1", nodes, elements, physical_groups, node_data_views)


def _resolve_bc(phys_tag, physical_groups):
    if phys_tag is None:
        return BoundaryType.MAINLAND
    pg = physical_groups.get(phys_tag)
    if pg is None:
        return int(phys_tag)
    if pg.name is None:
        return int(phys_tag)
    mapped = PHYSICAL_GROUP_MAP.get(pg.name.lower())
    if mapped is not None:
        return mapped
    warnings.warn(
        f"physical group name {pg.name!r} not in PHYSICAL_GROUP_MAP; "
        f"preserving numeric tag {phys_tag}",
        UserWarning, stacklevel=5,
    )
    return int(phys_tag)


def _pick_next_node(cur, prev, candidates, nodes):
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]
    if prev is None:
        return candidates[0]
    p_prev = nodes[prev]
    p_cur = nodes[cur]
    v_in = p_cur - p_prev
    best, best_angle = None, float("inf")
    for c in candidates:
        v_out = nodes[c] - p_cur
        angle = float(np.arctan2(
            v_in[0] * v_out[1] - v_in[1] * v_out[0],
            v_in[0] * v_out[0] + v_in[1] * v_out[1],
        ))
        if angle < best_angle:
            best_angle, best = angle, c
    return best


def _walk_rings(edges, nodes):
    from collections import defaultdict
    next_nodes: dict[int, list[int]] = defaultdict(list)
    for a, b in edges:
        next_nodes[a].append(b)
    used: set[tuple] = set()
    rings: list[list[int]] = []
    for start in sorted(next_nodes.keys()):
        for first in list(next_nodes[start]):
            if (start, first) in used:
                continue
            ring = [start]
            prev, cur = start, first
            used.add((start, first))
            guard = len(next_nodes) * 2 + 4
            while cur != start and guard > 0:
                ring.append(cur)
                cands = [v for v in next_nodes.get(cur, []) if (cur, v) not in used]
                nxt = _pick_next_node(cur, prev, cands, nodes)
                if nxt is None:
                    break
                used.add((cur, nxt))
                prev, cur = cur, nxt
                guard -= 1
            if cur == start and len(ring) >= 2:
                rings.append(ring)
            elif len(ring) >= 2:
                # Non-closing chain (open boundary) — still return as-is
                rings.append(ring)
    return rings


def _build_boundary_segments(line_elems, tag_to_idx, physical_groups, nodes):
    from collections import defaultdict
    from admesh.api import _ring_area
    tag_edges: dict = defaultdict(list)
    for e in line_elems:
        a0 = tag_to_idx.get(e.node_ids[0])
        b0 = tag_to_idx.get(e.node_ids[1])
        if a0 is None or b0 is None:
            continue
        tag_edges[e.phys_tag].append((a0, b0))
    segments: list[BoundarySegment] = []
    for phys_tag, edges in tag_edges.items():
        bc = _resolve_bc(phys_tag, physical_groups)
        rings = _walk_rings(edges, nodes)
        for ring in rings:
            is_open = isinstance(bc, BoundaryType) and bc == BoundaryType.OPEN
            segments.append(BoundarySegment(
                node_ids=np.array(ring, dtype=np.int64),
                bc_type=bc,
                is_open=is_open,
            ))
    segments.sort(key=lambda s: _ring_area(s.node_ids.tolist(), nodes), reverse=True)
    return segments


def _to_mesh(gf: _GmshFile) -> Mesh:
    sorted_tags = sorted(gf.nodes.keys())
    tag_to_idx = {t: i for i, t in enumerate(sorted_tags)}
    n = len(sorted_tags)
    nodes = np.empty((n, 2), dtype=np.float64)
    for i, t in enumerate(sorted_tags):
        x, y, z = gf.nodes[t]
        if abs(z) >= 1e-12:
            raise GmshParseError(f"non-planar mesh; node {t} has z={z}")
        nodes[i, 0] = x
        nodes[i, 1] = y

    tris = [e for e in gf.elements if e.etype == 2]
    elements = np.empty((len(tris), 3), dtype=np.int64)
    for i, e in enumerate(tris):
        for j, nid in enumerate(e.node_ids):
            elements[i, j] = tag_to_idx[nid]

    lines = [e for e in gf.elements if e.etype == 1]
    boundaries = _build_boundary_segments(lines, tag_to_idx, gf.physical_groups, nodes)

    bathymetry = None
    for v in gf.node_data_views:
        if v.name.lower() == "bathymetry":
            vals = v.values  # dict tag -> float
            arr = np.zeros(n, dtype=np.float64)
            for t, val in vals.items():
                idx = tag_to_idx.get(t)
                if idx is not None:
                    arr[idx] = val
            bathymetry = -arr  # sign flip: positive-down → positive-up
            break

    mesh = Mesh(nodes=nodes, elements=elements, boundaries=tuple(boundaries),
                bathymetry=bathymetry, quality=None, title="")
    if not mesh.boundaries and mesh.elements.size > 0:
        from admesh.api import _derive_boundary_segments
        bds = _derive_boundary_segments(mesh.elements, mesh.nodes, default_bc=BoundaryType.MAINLAND)
        mesh = Mesh(nodes=nodes, elements=elements, boundaries=bds,
                    bathymetry=bathymetry, quality=None, title="")
    return mesh


def read_msh(path: str | os.PathLike) -> Mesh:
    """Parse a Gmsh ASCII .msh file (v2.2 or v4.1) and return a Mesh."""
    path = os.fspath(path)
    with open(path, encoding="utf-8") as fh:
        lines = [line.rstrip("\n") for line in fh]
    cursor = _Cursor(lines)
    # Find $MeshFormat
    while True:
        line = cursor.next()
        if line is None:
            raise GmshParseError("empty or missing $MeshFormat", line_number=1)
        s = line.strip()
        if s == "$MeshFormat":
            break
        if s:
            raise GmshParseError(f"expected $MeshFormat, got {s!r}", line_number=cursor.line_no)
    fmt = cursor.require("version file-type data-size").strip().split()
    if len(fmt) < 3:
        raise GmshParseError(f"malformed $MeshFormat: {fmt}", line_number=cursor.line_no)
    version, file_type = fmt[0], int(fmt[1])
    if file_type != 0:
        raise GmshParseError("binary .msh files are not supported", line_number=cursor.line_no)
    if version not in ("2.2", "4.1"):
        raise GmshParseError(f"unsupported Gmsh version {version!r}", line_number=cursor.line_no)
    cursor.require("$EndMeshFormat")
    if version == "2.2":
        gf = _parse_v22(cursor)
    else:
        gf = _parse_v41(cursor)
    if not gf.physical_groups:
        warnings.warn(
            "no $PhysicalNames block found in .msh file — "
            "falling back to single MAINLAND boundary ring",
            UserWarning, stacklevel=2,
        )
    return _to_mesh(gf)


# --- Writer helpers ---

def _collect_boundary_tags(mesh):
    bc_to_tag: dict[int, int] = {}
    tag_to_name: dict[int, str] = {}
    seg_tags: list = []
    next_tag = 1
    for seg in mesh.boundaries:
        bc_int = int(seg.bc_type)
        if bc_int not in bc_to_tag:
            name = _BOUNDARY_TYPE_TO_NAME.get(bc_int, f"ibtype_{bc_int}")
            bc_to_tag[bc_int] = next_tag
            tag_to_name[next_tag] = name
            next_tag += 1
        seg_tags.append((seg, bc_to_tag[bc_int]))
    return tag_to_name, seg_tags


def _write_node_data(mesh: Mesh, out) -> None:
    n = mesh.n_nodes
    out.write("$NodeData\n1\n\"bathymetry\"\n1\n0.0\n3\n0\n1\n")
    out.write(f"{n}\n")
    for i in range(n):
        out.write(f"{i + 1} {-float(mesh.bathymetry[i]):.15g}\n")
    out.write("$EndNodeData\n")


def _write_v22(mesh: Mesh, out) -> None:
    tag_to_name, seg_tags = _collect_boundary_tags(mesh)
    out.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
    surf_tag = (max(tag_to_name.keys()) + 1) if tag_to_name else 1
    out.write("$PhysicalNames\n")
    out.write(f"{len(tag_to_name) + 1}\n")
    for tag, name in tag_to_name.items():
        out.write(f'1 {tag} "{name}"\n')
    out.write(f'2 {surf_tag} "domain"\n')
    out.write("$EndPhysicalNames\n")
    out.write(f"$Nodes\n{mesh.n_nodes}\n")
    for i in range(mesh.n_nodes):
        out.write(f"{i+1} {float(mesh.nodes[i,0]):.15g} {float(mesh.nodes[i,1]):.15g} 0\n")
    out.write("$EndNodes\n")
    n_line = sum(len(seg.node_ids) for seg, _ in seg_tags)
    total = n_line + mesh.n_elements
    out.write(f"$Elements\n{total}\n")
    eid = 1
    for seg, ptag in seg_tags:
        ring = seg.node_ids
        n = len(ring)
        for k in range(n):
            a, b = int(ring[k]) + 1, int(ring[(k+1)%n]) + 1
            out.write(f"{eid} 1 2 {ptag} 1 {a} {b}\n")
            eid += 1
    for i in range(mesh.n_elements):
        n0, n1, n2 = int(mesh.elements[i,0])+1, int(mesh.elements[i,1])+1, int(mesh.elements[i,2])+1
        out.write(f"{eid} 2 2 {surf_tag} 1 {n0} {n1} {n2}\n")
        eid += 1
    out.write("$EndElements\n")
    if mesh.bathymetry is not None:
        _write_node_data(mesh, out)


def _write_v41(mesh: Mesh, out) -> None:
    tag_to_name, seg_tags = _collect_boundary_tags(mesh)
    out.write("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n")
    n_segs = len(seg_tags)
    surf_ent_tag = n_segs + 1
    surf_phys_tag = (max(tag_to_name.keys()) + 1) if tag_to_name else 1
    out.write("$PhysicalNames\n")
    out.write(f"{len(tag_to_name) + 1}\n")
    for tag, name in tag_to_name.items():
        out.write(f'1 {tag} "{name}"\n')
    out.write(f'2 {surf_phys_tag} "domain"\n')
    out.write("$EndPhysicalNames\n")
    out.write("$Entities\n")
    out.write(f"0 {n_segs} 1 0\n")
    for ci, (seg, ptag) in enumerate(seg_tags):
        sn = mesh.nodes[seg.node_ids]
        xmn, ymn = float(sn[:,0].min()), float(sn[:,1].min())
        xmx, ymx = float(sn[:,0].max()), float(sn[:,1].max())
        ctag = ci + 1
        out.write(f"{ctag} {xmn:.15g} {ymn:.15g} 0 {xmx:.15g} {ymx:.15g} 0 1 {ptag} 0\n")
    xmn, ymn = float(mesh.nodes[:,0].min()), float(mesh.nodes[:,1].min())
    xmx, ymx = float(mesh.nodes[:,0].max()), float(mesh.nodes[:,1].max())
    curve_list = " ".join(str(i+1) for i in range(n_segs))
    out.write(f"{surf_ent_tag} {xmn:.15g} {ymn:.15g} 0 {xmx:.15g} {ymx:.15g} 0 1 {surf_phys_tag} {n_segs} {curve_list}\n")
    out.write("$EndEntities\n")
    n_nodes = mesh.n_nodes
    out.write("$Nodes\n")
    out.write(f"1 {n_nodes} 1 {n_nodes}\n")
    out.write(f"2 {surf_ent_tag} 0 {n_nodes}\n")
    for i in range(n_nodes):
        out.write(f"{i+1}\n")
    for i in range(n_nodes):
        out.write(f"{float(mesh.nodes[i,0]):.15g} {float(mesh.nodes[i,1]):.15g} 0\n")
    out.write("$EndNodes\n")
    n_line = sum(len(seg.node_ids) for seg, _ in seg_tags)
    total = n_line + mesh.n_elements
    n_blocks = n_segs + 1
    out.write("$Elements\n")
    out.write(f"{n_blocks} {total} 1 {total}\n")
    eid = 1
    for ci, (seg, ptag) in enumerate(seg_tags):
        ctag = ci + 1
        ring = seg.node_ids
        n = len(ring)
        out.write(f"1 {ctag} 1 {n}\n")
        for k in range(n):
            a, b = int(ring[k])+1, int(ring[(k+1)%n])+1
            out.write(f"{eid} {a} {b}\n")
            eid += 1
    out.write(f"2 {surf_ent_tag} 2 {mesh.n_elements}\n")
    for i in range(mesh.n_elements):
        n0, n1, n2 = int(mesh.elements[i,0])+1, int(mesh.elements[i,1])+1, int(mesh.elements[i,2])+1
        out.write(f"{eid} {n0} {n1} {n2}\n")
        eid += 1
    out.write("$EndElements\n")
    if mesh.bathymetry is not None:
        _write_node_data(mesh, out)


def write_msh(mesh: Mesh, path: str | os.PathLike, *, version: str = "4.1") -> None:
    """Serialize a Mesh as a Gmsh ASCII .msh file."""
    if version not in ("2.2", "4.1"):
        raise ValueError(f"version must be '2.2' or '4.1', got {version!r}")
    with open(os.fspath(path), "w", encoding="utf-8") as out:
        if version == "2.2":
            _write_v22(mesh, out)
        else:
            _write_v41(mesh, out)
