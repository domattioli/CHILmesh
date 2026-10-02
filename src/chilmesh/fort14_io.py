"""ADCIRC fort.14 (mesh) byte-faithful I/O for CHILmesh.

``CHILmesh.read_from_fort14`` is lossy by design for identity-preserving
consumers: it discards the fort.14 node/element **id columns** (positional
only), normalizes element **winding** to CCW in the constructor, and **pads**
mixed tri/quad connectivity to 4 columns. This module is the faithful
counterpart — it preserves id columns, original winding, true per-element
arity (no padding), and the full NOPE/NBOU boundary block including
barrier/weir columns — so a downstream identity/registry layer (Valence #214)
can re-derive byte-identical digests. The CCW-normalizing ``CHILmesh``
constructor is left untouched.

Additive module on the ``fort13_io`` / ``fort15_io`` precedent.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

# ADCIRC NBOU flow-boundary node-line column layouts, by IBTYPE:
#   1-sided (external) barrier {3,13,23}: NBVV BARLANHT BARLANCFSP
#   2-sided (internal) barrier {4,24,64} and pipe barriers {5,25}:
#       NBVV IBCONN BARINHT BARINCFSB BARINCFSP [PIPEHT PIPECOEF PIPEDIAM]
_IBTYPE_1SIDED_BARRIER = frozenset({3, 13, 23})
_IBTYPE_2SIDED_BARRIER = frozenset({4, 24, 64, 5, 25})


def _lead_int(tokens, idx):
    """Return int(tokens[idx]) if that token is an integer literal, else None.

    Tolerates ADCIRC headers with trailing prose ('24 = Number of nodes ...')
    by extracting only the leading integer token and discarding the rest.
    """
    if idx >= len(tokens):
        return None
    tok = tokens[idx]
    try:
        return int(tok)
    except (TypeError, ValueError):
        try:
            return int(float(tok))
        except (TypeError, ValueError):
            return None


class Fort14ParseError(ValueError):
    """Raised when parsing a fort.14 file encounters an error."""


@dataclass
class OpenBoundary:
    """An ADCIRC NOPE open-boundary segment (faithful)."""
    nodes: list                     # NBDV node ids, as written (1-based)
    ibtype: Optional[int] = None    # IBTYPEE (usually absent)


@dataclass
class FlowBoundary:
    """An ADCIRC NBOU flow/land-boundary segment (faithful).

    ``nodes`` is always the primary node column (NBVV). For internal-barrier
    ibtypes, ``back_nodes`` is the paired node (IBCONN) and ``heights`` the
    crest height; ``coeffs`` holds the remaining per-node coefficient columns.
    ``node_extra`` retains the raw trailing tokens per node verbatim so the
    segment round-trips losslessly regardless of ibtype.
    """
    ibtype: Optional[int]
    nodes: list
    back_nodes: Optional[list] = None
    heights: Optional[list] = None
    coeffs: Optional[list] = None
    node_extra: list = field(default_factory=list)


@dataclass
class Fort14Raw:
    """Faithful in-memory record of an ADCIRC fort.14 file."""
    grid_name: str
    n_nodes: int
    n_elems: int
    node_ids: list                  # column-0 node ids, in file order
    coords: dict                    # node id -> (x, y, depth)
    elem_ids: list                  # column-0 element ids, in file order
    elements: dict                  # element id -> tuple of node ids (true arity, file winding)
    open_boundaries: list = field(default_factory=list)
    flow_boundaries: list = field(default_factory=list)
    boundaries_present: bool = False  # True iff a NOPE/NBOU block was physically parsed (#259)


class _Cursor:
    """Mutable line index shared by the boundary-block helpers.

    The ``near line N`` text in :class:`Fort14ParseError` must report the
    position at the moment of failure, so helpers advance one shared index
    instead of returning a new one.
    """

    __slots__ = ("i",)

    def __init__(self, i: int) -> None:
        self.i = i


def _read_node_rows(lines: list, i: int, n_nodes: int, path) -> tuple:
    """Parse the node block starting at line ``i``.

    Returns
    -------
    tuple
        ``(node_ids, coords, next_line_index)``.

    Raises
    ------
    Fort14ParseError
        On a short or non-numeric node row.
    """
    node_ids: list = []
    coords: dict = {}
    try:
        for _ in range(n_nodes):
            p = lines[i].split(); i += 1
            nid = int(float(p[0]))
            x, y = float(p[1]), float(p[2])
            dp = float(p[3]) if len(p) > 3 else 0.0
            node_ids.append(nid)
            coords[nid] = (x, y, dp)
    except (IndexError, ValueError) as e:
        raise Fort14ParseError(f"{path}: bad node row near line {i}: {e}") from e
    return node_ids, coords, i


def _read_element_rows(lines: list, i: int, n_elems: int, path) -> tuple:
    """Parse the element block starting at line ``i``.

    Returns
    -------
    tuple
        ``(elem_ids, elements, next_line_index)``.

    Raises
    ------
    Fort14ParseError
        On a short row, non-numeric token, or an arity mismatch.
    """
    elem_ids: list = []
    elements: dict = {}
    try:
        for _ in range(n_elems):
            p = lines[i].split(); i += 1
            eid = int(float(p[0])); nhy = int(float(p[1]))
            verts = tuple(int(float(v)) for v in p[2:2 + nhy])
            if len(verts) != nhy:
                raise ValueError(f"element {eid} declares {nhy} nodes, found {len(verts)}")
            elem_ids.append(eid)
            elements[eid] = verts
    except (IndexError, ValueError) as e:
        raise Fort14ParseError(f"{path}: bad element row near line {i}: {e}") from e
    return elem_ids, elements, i


def _parse_open_segment(lines: list, cur: _Cursor) -> OpenBoundary:
    """Parse one NOPE segment: a header row then NVDLL node rows.

    Stops early (without error) if the file ends mid-segment, mirroring the
    tolerant legacy behaviour.
    """
    parts = lines[cur.i].split(); cur.i += 1
    nvdll = int(parts[0])
    ibtypee = _lead_int(parts, 1)
    nodes = []
    for _ in range(nvdll):
        if cur.i >= len(lines):
            break
        nodes.append(int(lines[cur.i].split()[0])); cur.i += 1
    return OpenBoundary(nodes=nodes, ibtype=ibtypee)


def _parse_open_boundaries(lines: list, cur: _Cursor, nope: int) -> list:
    """Parse ``nope`` open-boundary segments, stopping at end of file."""
    segments = []
    for _ in range(nope):
        if cur.i >= len(lines):
            break
        segments.append(_parse_open_segment(lines, cur))
    return segments


def _accumulate_flow_node(toks: list, ibtype, nodes, back_nodes, heights, coeffs, node_extra) -> None:
    """Append one flow-boundary node row to the per-segment accumulators.

    Barrier ibtypes carry extra numeric columns that are decoded into
    ``back_nodes`` / ``heights`` / ``coeffs``; the raw trailing tokens are
    always kept in ``node_extra`` so the row round-trips verbatim.
    """
    nodes.append(int(toks[0]))
    rest = toks[1:]
    node_extra.append(rest)
    if ibtype in _IBTYPE_2SIDED_BARRIER:
        if len(rest) >= 1:
            back_nodes.append(int(float(rest[0])))
        if len(rest) >= 2:
            heights.append(float(rest[1]))
        coeffs.append([float(t) for t in rest[2:]])
    elif ibtype in _IBTYPE_1SIDED_BARRIER:
        if len(rest) >= 1:
            heights.append(float(rest[0]))
        coeffs.append([float(t) for t in rest[1:]])


def _parse_flow_segment(lines: list, cur: _Cursor) -> FlowBoundary:
    """Parse one NBOU segment: a header row then NVELL node rows."""
    parts = lines[cur.i].split(); cur.i += 1
    nvell = int(parts[0])
    ibtype = _lead_int(parts, 1)
    nodes: list = []
    back_nodes: list = []
    heights: list = []
    coeffs: list = []
    node_extra: list = []
    for _ in range(nvell):
        if cur.i >= len(lines):
            break
        toks = lines[cur.i].split(); cur.i += 1
        _accumulate_flow_node(toks, ibtype, nodes, back_nodes, heights, coeffs, node_extra)
    return FlowBoundary(
        ibtype=ibtype,
        nodes=nodes,
        back_nodes=back_nodes or None,
        heights=heights or None,
        coeffs=coeffs or None,
        node_extra=node_extra,
    )


def _parse_flow_boundaries(lines: list, cur: _Cursor) -> list:
    """Parse the NBOU block; returns ``[]`` if the file ended after NOPE."""
    segments: list = []
    # Guard: if file truncated after open boundaries, skip flow boundaries
    if cur.i >= len(lines):
        return segments
    nbou = int(lines[cur.i].split()[0]); cur.i += 1
    if cur.i < len(lines):
        cur.i += 1  # NVEL (total flow nodes)
    for _ in range(nbou):
        if cur.i >= len(lines):
            break
        segments.append(_parse_flow_segment(lines, cur))
    return segments


def _read_boundary_block(lines: list, i: int, path) -> tuple:
    """Parse the NOPE/NBOU block starting at line ``i``.

    Returns
    -------
    tuple
        ``(open_boundaries, flow_boundaries)``.

    Raises
    ------
    Fort14ParseError
        On a malformed boundary row.
    """
    cur = _Cursor(i)
    try:
        nope = int(lines[cur.i].split()[0]); cur.i += 1
        cur.i += 1  # NETA (total open nodes) — recomputed on write
        open_boundaries = _parse_open_boundaries(lines, cur, nope)
        flow_boundaries = _parse_flow_boundaries(lines, cur)
    except (IndexError, ValueError) as e:
        raise Fort14ParseError(f"{path}: malformed boundary block near line {cur.i}: {e}") from e
    return open_boundaries, flow_boundaries


def read_fort14_raw(filename, parse_boundaries: bool = True) -> Fort14Raw:
    """Read an ADCIRC fort.14 preserving ids, winding, arity, and boundaries.

    Parameters:
        filename: path to the fort.14 file.
        parse_boundaries: if False, skip the boundary block and return empty
            ``open_boundaries`` / ``flow_boundaries``. Default True preserves
            all boundary data for byte-faithful I/O.

    Raises:
        Fort14ParseError: on a malformed header, node/element block, or
            boundary block. A legacy mesh with no boundary block is accepted
            (empty ``open_boundaries`` / ``flow_boundaries``).

    Note:
        ``boundaries_present`` is ``True`` only when a NOPE/NBOU block was
        physically parsed, distinguishing a present-but-empty (0/0) boundary
        section from an absent one (#259). It stays ``False`` when
        ``parse_boundaries=False`` (the block was skipped, not inspected).
    """
    path = Path(filename)
    with open(path, encoding="utf-8") as fh:
        lines = fh.readlines()
    if len(lines) < 2:
        raise Fort14ParseError(f"{path}: too short to be a fort.14 file")

    i = 0
    grid_name = lines[i].rstrip("\n"); i += 1
    try:
        hdr = lines[i].split(); i += 1
        n_elems, n_nodes = int(float(hdr[0])), int(float(hdr[1]))
    except (IndexError, ValueError) as e:
        raise Fort14ParseError(f"{path}: bad NE/NP header: {e}") from e

    node_ids, coords, i = _read_node_rows(lines, i, n_nodes, path)
    elem_ids, elements, i = _read_element_rows(lines, i, n_elems, path)

    open_boundaries: list = []
    flow_boundaries: list = []
    boundaries_present = False
    if parse_boundaries and i < len(lines) and lines[i].strip():
        boundaries_present = True
        open_boundaries, flow_boundaries = _read_boundary_block(lines, i, path)

    return Fort14Raw(
        grid_name=grid_name, n_nodes=n_nodes, n_elems=n_elems,
        node_ids=node_ids, coords=coords, elem_ids=elem_ids, elements=elements,
        open_boundaries=open_boundaries, flow_boundaries=flow_boundaries,
        boundaries_present=boundaries_present,
    )


def write_fort14_raw(raw: Fort14Raw, filename) -> bool:
    """Write a :class:`Fort14Raw` back to ADCIRC fort.14 format.

    Node/element id columns, winding, and arity are emitted as stored; flow
    boundary node lines are reproduced from ``node_extra`` (raw tokens), so a
    record obtained from :func:`read_fort14_raw` round-trips losslessly.

    Returns ``True`` on success.
    """
    path = Path(filename)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(f"{raw.grid_name}\n")
        f.write(f"{raw.n_elems} {raw.n_nodes}\n")
        for nid in raw.node_ids:
            x, y, dp = raw.coords[nid]
            f.write(f"{nid} {x:.10f} {y:.10f} {dp:.10f}\n")
        for eid in raw.elem_ids:
            verts = raw.elements[eid]
            f.write(f"{eid} {len(verts)} {' '.join(str(v) for v in verts)}\n")

        # NOPE open boundaries
        f.write(f"{len(raw.open_boundaries)}\n")
        f.write(f"{sum(len(b.nodes) for b in raw.open_boundaries)}\n")
        for b in raw.open_boundaries:
            hdr = f"{len(b.nodes)}" + (f" {b.ibtype}" if b.ibtype is not None else "")
            f.write(hdr + "\n")
            for n in b.nodes:
                f.write(f"{n}\n")

        # NBOU flow boundaries
        f.write(f"{len(raw.flow_boundaries)}\n")
        f.write(f"{sum(len(b.nodes) for b in raw.flow_boundaries)}\n")
        for b in raw.flow_boundaries:
            hdr = f"{len(b.nodes)}" + (f" {b.ibtype}" if b.ibtype is not None else "")
            f.write(hdr + "\n")
            for j, n in enumerate(b.nodes):
                extra = b.node_extra[j] if (b.node_extra and j < len(b.node_extra)) else []
                if extra:
                    f.write(f"{n} " + " ".join(str(t) for t in extra) + "\n")
                else:
                    f.write(f"{n}\n")
    return True


__all__ = [
    "Fort14Raw", "OpenBoundary", "FlowBoundary",
    "read_fort14_raw", "write_fort14_raw", "Fort14ParseError",
]
