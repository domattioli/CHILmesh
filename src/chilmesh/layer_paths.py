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

    Per-layer cost: ``O(|OE|+|IE|+|E_layer|)`` -- scans only the layer's
    elements and their incident edges, not the full mesh. Summed across
    all layers this is ``O(m)`` because the skeletonization invariant
    guarantees each element belongs to exactly one layer.

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

    oe = np.asarray(mesh.layers["OE"][layer_idx], dtype=int)
    ie = np.asarray(mesh.layers["IE"][layer_idx], dtype=int)
    layer_elems_arr = np.concatenate([oe, ie])

    edge2vert = mesh.adjacencies["Edge2Vert"]
    edge2elem = mesh.adjacencies["Edge2Elem"]

    # O(1) membership bitmaps replace set lookups.
    ov_mask = np.zeros(mesh.n_verts, dtype=bool)
    ov_mask[ov] = True
    layer_elem_mask = np.zeros(mesh.n_elems, dtype=bool)
    layer_elem_mask[layer_elems_arr] = True

    # Step 1: scope candidate edges to those touching a layer element.
    # Mixed meshes pad triangle rows of Elem2Edge with edge id 0; the
    # downstream "exactly one layer neighbour" check below filters any
    # spurious inclusions, so we don't need to scrub padding here.
    layer_eids = np.unique(mesh.elem2edge(layer_elems_arr).ravel())

    # Step 2: bitmap-filter to edges with BOTH endpoints in OV.
    e2v_subset = edge2vert[layer_eids]
    both_in_ov = ov_mask[e2v_subset[:, 0]] & ov_mask[e2v_subset[:, 1]]
    candidate_eids = layer_eids[both_in_ov]
    if candidate_eids.size == 0:
        return {int(v): [] for v in ov}, {}

    # Step 3: vectorised layer-element neighbour count.
    e2e_subset = edge2elem[candidate_eids]
    in_layer = np.zeros_like(e2e_subset, dtype=bool)
    valid = e2e_subset >= 0
    in_layer[valid] = layer_elem_mask[e2e_subset[valid]]
    is_layer_boundary = in_layer.sum(axis=1) == 1

    final_eids = candidate_eids[is_layer_boundary]
    final_e2v = edge2vert[final_eids]
    # For each kept edge, pick the (unique) layer-element neighbour.
    final_e2e = e2e_subset[is_layer_boundary]
    final_in_layer = in_layer[is_layer_boundary]
    layer_neighbour_ids = final_e2e[
        np.arange(final_e2e.shape[0]), final_in_layer.argmax(axis=1)
    ]

    adj: Dict[int, List[Tuple[int, int]]] = {int(v): [] for v in ov}
    edge_layer_elems: Dict[int, Set[int]] = {}
    for i in range(final_eids.shape[0]):
        eid_int = int(final_eids[i])
        u, v = int(final_e2v[i, 0]), int(final_e2v[i, 1])
        adj[u].append((v, eid_int))
        adj[v].append((u, eid_int))
        edge_layer_elems[eid_int] = {int(layer_neighbour_ids[i])}

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
