"""Ordered vertex paths along the outer boundary of a skeletonization layer.

Idiomatic Python port of MATLAB ``PathsOnOV.m`` from the legacy
``QuADMesh-MATLAB/02_QuADMESH_Library/03_Layer_Paths`` library.

Graph-theoretic summary
-----------------------
For a given mesh layer ``k``, the algorithm operates on the subgraph
``G = (V, E)`` where

* ``V = layers["OV"][k]`` -- vertices on the layer's outer ring, and
* ``E`` -- edges of the layer's sub-mesh (elements ``OE[k] ∪ IE[k]``) that
  appear in exactly one layer-element AND whose both endpoints lie in ``V``.

``G`` is the layer's boundary ring graph. For a simply-connected annular
layer it is a single cycle; for a donut it is two disjoint cycles; for
layers with pinch points it is a multigraph with degree-≥3 "junction"
vertices where the medial axis branches.

The function decomposes ``G`` into an ordered set of closed walks
(``Path[0]``, ``Path[1]``, …) covering every edge of ``G`` exactly once.
This is a greedy Eulerian-style decomposition with two heuristics applied
when more than one outgoing edge is available at a junction:

1. Prefer the next edge that shares a layer-element with the previous
   path-edge (keeps consecutive path-edges on the same face).
2. Failing that, prefer the smallest turning angle (locally-straight path).

These match the intent of the MATLAB original; the order in which paths
are seeded across components is intentionally simpler (junction-first,
then any remaining vertex) but the resulting cover is graph-theoretically
equivalent.

Downstream heuristics (e.g. the Tri2Quad routine) consume the returned
paths to walk vertices in a consistent order around each layer.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Dict, List, Set, Tuple

import numpy as np

if TYPE_CHECKING:  # pragma: no cover -- type-only import to avoid circularity
    from .CHILmesh import CHILmesh


__all__ = ["paths_on_outer_vertices"]


def paths_on_outer_vertices(mesh: "CHILmesh", layer_idx: int) -> List[np.ndarray]:
    """Return ordered vertex paths covering a layer's outer-boundary edges.

    Args:
        mesh: A :class:`CHILmesh` whose layers have been computed.
        layer_idx: Zero-based layer index (0 = outermost).

    Returns:
        A list of 1-D ``ndarray``s of global vertex IDs, one per path. Each
        path is a sequence ``v_0, v_1, ..., v_k`` where consecutive vertices
        are joined by an edge of the layer's outer-vertex subgraph. Closed
        walks have ``v_0 == v_k``; the duplicate is intentional so the path
        can be interpreted as either a vertex sequence or an edge sequence.

    Raises:
        IndexError: If ``layer_idx`` is outside ``[0, n_layers)``.
    """
    if layer_idx < 0 or layer_idx >= mesh.n_layers:
        raise IndexError(
            f"layer_idx {layer_idx} out of range [0, {mesh.n_layers})"
        )

    adj, edge_layer_elems = _build_outer_vertex_subgraph(mesh, layer_idx)
    if not edge_layer_elems:
        return []

    return _decompose_into_paths(adj, edge_layer_elems, mesh.points)


def _build_outer_vertex_subgraph(
    mesh: "CHILmesh", layer_idx: int
) -> Tuple[Dict[int, List[Tuple[int, int]]], Dict[int, Set[int]]]:
    """Build the outer-vertex subgraph for layer ``layer_idx``.

    Returns
    -------
    adj
        ``{vertex_id: [(neighbor_id, edge_id), ...]}`` symmetric adjacency.
    edge_layer_elems
        ``{edge_id: {layer_element_id}}`` -- the (single) layer element on
        the inside of each graph edge. Used as the "shares-face" key for
        the junction heuristic.
    """
    ov = np.asarray(mesh.layers["OV"][layer_idx], dtype=int)
    if ov.size == 0:
        return {}, {}
    ov_set: Set[int] = {int(v) for v in ov}

    oe = np.asarray(mesh.layers["OE"][layer_idx], dtype=int)
    ie = np.asarray(mesh.layers["IE"][layer_idx], dtype=int)
    layer_elems: Set[int] = {int(e) for e in np.concatenate([oe, ie])}

    edge2vert = mesh.adjacencies["Edge2Vert"]
    edge2elem = mesh.adjacencies["Edge2Elem"]

    # Vectorised filter: keep edges with both endpoints in OV.
    both_in_ov = np.isin(edge2vert[:, 0], ov) & np.isin(edge2vert[:, 1], ov)
    candidate_eids = np.where(both_in_ov)[0]

    adj: Dict[int, List[Tuple[int, int]]] = {int(v): [] for v in ov}
    edge_layer_elems: Dict[int, Set[int]] = {}

    for eid in candidate_eids:
        eid_int = int(eid)
        e2e = edge2elem[eid_int]
        layer_neighbours = {int(x) for x in e2e if x >= 0 and int(x) in layer_elems}
        # An edge belongs to the layer's boundary iff exactly one of its
        # incident elements is in this layer (the other is in a different
        # layer or outside the mesh).
        if len(layer_neighbours) != 1:
            continue
        u, v = int(edge2vert[eid_int, 0]), int(edge2vert[eid_int, 1])
        adj[u].append((v, eid_int))
        adj[v].append((u, eid_int))
        edge_layer_elems[eid_int] = layer_neighbours

    return adj, edge_layer_elems


def _decompose_into_paths(
    adj: Dict[int, List[Tuple[int, int]]],
    edge_layer_elems: Dict[int, Set[int]],
    points: np.ndarray,
) -> List[np.ndarray]:
    """Greedy Eulerian-style decomposition of the outer-vertex subgraph.

    Each iteration picks a seed vertex (junctions first), walks unvisited
    edges using the two-tier heuristic, and stops when the walk closes back
    on its seed or runs out of unvisited outgoing edges. Repeat until every
    edge has been visited.
    """
    visited: Set[int] = set()
    n_edges = len(edge_layer_elems)
    paths: List[np.ndarray] = []

    def degree_remaining(v: int) -> int:
        return sum(1 for _, e in adj[v] if e not in visited)

    while len(visited) < n_edges:
        seed = _pick_seed(adj, visited, degree_remaining)
        if seed is None:  # pragma: no cover -- guarded by while condition
            break

        path = _walk(seed, adj, edge_layer_elems, visited, points)
        if len(path) <= 1:  # pragma: no cover -- pick_seed guarantees degree>0
            break
        paths.append(np.asarray(path, dtype=int))

    return paths


def _pick_seed(
    adj: Dict[int, List[Tuple[int, int]]],
    visited: Set[int],
    degree_remaining,
) -> int | None:
    """Choose a seed vertex: prefer junctions (degree > 2 in remaining graph)."""
    best_any: int | None = None
    best_junction: int | None = None
    best_junction_deg = -1
    for v, edges in adj.items():
        deg = degree_remaining(v)
        if deg == 0:
            continue
        if best_any is None:
            best_any = v
        if deg > 2 and deg > best_junction_deg:
            best_junction = v
            best_junction_deg = deg
    return best_junction if best_junction is not None else best_any


def _walk(
    start: int,
    adj: Dict[int, List[Tuple[int, int]]],
    edge_layer_elems: Dict[int, Set[int]],
    visited: Set[int],
    points: np.ndarray,
) -> List[int]:
    """Walk unvisited edges from ``start`` until closure or dead-end."""
    path: List[int] = [start]
    prev_edge: int | None = None
    prev_v: int | None = None
    cur = start

    while True:
        outgoing = [(n, e) for n, e in adj[cur] if e not in visited]
        if not outgoing:
            break
        nbr, eid = _pick_next(cur, prev_v, prev_edge, outgoing,
                               edge_layer_elems, points)
        visited.add(eid)
        path.append(nbr)
        prev_v = cur
        prev_edge = eid
        cur = nbr
        if cur == start:
            break
    return path


def _pick_next(
    cur: int,
    prev_v: int | None,
    prev_edge: int | None,
    outgoing: List[Tuple[int, int]],
    edge_layer_elems: Dict[int, Set[int]],
    points: np.ndarray,
) -> Tuple[int, int]:
    """Apply the two-tier heuristic to choose the next edge at ``cur``."""
    if len(outgoing) == 1 or prev_edge is None:
        return outgoing[0]

    # Heuristic 1: prefer next edge sharing a layer-element with prev_edge.
    prev_elems = edge_layer_elems[prev_edge]
    same_face = [(n, e) for n, e in outgoing if edge_layer_elems[e] & prev_elems]
    candidates = same_face if same_face else outgoing
    if len(candidates) == 1:
        return candidates[0]

    # Heuristic 2: minimise turning angle |∠(prev_v → cur → nbr)|.
    p_cur = points[cur]
    a = p_cur - points[prev_v]
    a_norm = float(np.linalg.norm(a))
    if a_norm == 0.0:  # pragma: no cover -- degenerate coincident points
        return candidates[0]

    def turning(nbr: int) -> float:
        b = points[nbr] - p_cur
        b_norm = float(np.linalg.norm(b))
        if b_norm == 0.0:  # pragma: no cover -- degenerate coincident points
            return np.pi
        cos_ang = float(np.dot(a, b) / (a_norm * b_norm))
        return float(np.arccos(np.clip(cos_ang, -1.0, 1.0)))

    return min(candidates, key=lambda ne: turning(ne[0]))
