"""Invariant tests for :func:`chilmesh.paths_on_outer_vertices`.

Validates the graph-theoretic contract of the layer-paths decomposition on
every small fixture (annulus, donut, structured, quad_2x2). ``block_o`` is
excluded for runtime; the invariants are not topology-specific so the small
fixtures already cover the simply-connected, multiply-connected, and
boundary-junction cases.
"""
from __future__ import annotations

import numpy as np
import pytest

from chilmesh import examples, paths_on_outer_vertices


FIXTURE_NAMES = ["annulus", "donut", "structured", "quad_2x2"]


def _layer_graph(mesh, layer_idx):
    """Recompute the outer-vertex subgraph the function should be covering.

    Returns
    -------
    ov_set : set[int]
        Outer vertices of the layer.
    graph_edges : set[int]
        Edge IDs in the layer's outer-vertex subgraph.
    """
    ov = mesh.layers["OV"][layer_idx]
    ov_set = {int(v) for v in ov}
    oe = mesh.layers["OE"][layer_idx]
    ie = mesh.layers["IE"][layer_idx]
    layer_elems = {int(e) for e in np.concatenate([oe, ie])}

    e2v = mesh.adjacencies["Edge2Vert"]
    e2e = mesh.adjacencies["Edge2Elem"]
    graph_edges = set()
    for eid in range(e2v.shape[0]):
        u, v = int(e2v[eid, 0]), int(e2v[eid, 1])
        if u not in ov_set or v not in ov_set:
            continue
        in_layer = [int(x) for x in e2e[eid] if x >= 0 and int(x) in layer_elems]
        if len(in_layer) == 1:
            graph_edges.add(eid)
    return ov_set, graph_edges


def _path_edges(mesh, path):
    """Convert a vertex-sequence path into the set of edge IDs it traverses."""
    e2v = mesh.adjacencies["Edge2Vert"]
    endpoint_pairs = {(min(int(a), int(b)), max(int(a), int(b))): eid
                      for eid, (a, b) in enumerate(e2v)}
    edges = []
    for i in range(len(path) - 1):
        a, b = int(path[i]), int(path[i + 1])
        key = (min(a, b), max(a, b))
        edges.append(endpoint_pairs[key])
    return edges


@pytest.fixture(params=FIXTURE_NAMES)
def small_mesh(request):
    return getattr(examples, request.param)()


def test_every_layer_edge_covered_exactly_once(small_mesh):
    """Each layer-boundary edge between two OV verts is in exactly one path."""
    for k in range(small_mesh.n_layers):
        _, expected_edges = _layer_graph(small_mesh, k)
        paths = small_mesh.paths_on_outer_vertices(k)

        seen: list[int] = []
        for p in paths:
            seen.extend(_path_edges(small_mesh, p))

        assert len(seen) == len(set(seen)), (
            f"layer {k}: edge visited more than once ({len(seen)} traversals, "
            f"{len(set(seen))} unique)"
        )
        assert set(seen) == expected_edges, (
            f"layer {k}: covered {len(set(seen))} edges but graph has "
            f"{len(expected_edges)}"
        )


def test_paths_traverse_only_outer_vertices(small_mesh):
    """Every vertex in every path must be in the layer's outer-vertex set."""
    for k in range(small_mesh.n_layers):
        ov_set, _ = _layer_graph(small_mesh, k)
        paths = small_mesh.paths_on_outer_vertices(k)
        for i, p in enumerate(paths):
            assert set(int(v) for v in p).issubset(ov_set), (
                f"layer {k} path {i}: contains non-outer vertices "
                f"{set(int(v) for v in p) - ov_set}"
            )


def test_consecutive_path_vertices_share_an_edge(small_mesh):
    """Every (v_i, v_{i+1}) pair in a path must be a real layer-graph edge."""
    for k in range(small_mesh.n_layers):
        _, graph_edges = _layer_graph(small_mesh, k)
        paths = small_mesh.paths_on_outer_vertices(k)
        e2v = small_mesh.adjacencies["Edge2Vert"]
        endpoint_pairs = {
            (min(int(a), int(b)), max(int(a), int(b))): eid
            for eid, (a, b) in enumerate(e2v)
        }
        for i, p in enumerate(paths):
            for j in range(len(p) - 1):
                a, b = int(p[j]), int(p[j + 1])
                key = (min(a, b), max(a, b))
                assert key in endpoint_pairs, (
                    f"layer {k} path {i}: ({a},{b}) is not a mesh edge"
                )
                assert endpoint_pairs[key] in graph_edges, (
                    f"layer {k} path {i}: edge ({a},{b}) not in layer graph"
                )


def test_closed_walks_when_topology_allows(small_mesh):
    """Paths that exit a vertex with even remaining degree should close.

    Concretely: for a connected component of the layer-graph with all
    even degrees (i.e. a single cycle, the common case), the decomposition
    must return one closed path. We check the simplest case: each path's
    start has degree >= 2 in the graph, so we expect closure unless the
    walk was terminated by a junction dead-end.
    """
    for k in range(small_mesh.n_layers):
        paths = small_mesh.paths_on_outer_vertices(k)
        adj = _adjacency_from_graph(small_mesh, k)
        for i, p in enumerate(paths):
            start, end = int(p[0]), int(p[-1])
            # If the start node has degree 2 in the full graph and all
            # connected component edges are in this path, it must be closed.
            if len(adj[start]) == 2 and start != end:
                # Allow non-closure only when the walk dead-ended at a
                # junction (which can happen on the second pass for non-
                # Eulerian components).
                assert len(adj[end]) > 2, (
                    f"layer {k} path {i}: simple-cycle start has degree 2 "
                    f"but path did not close ({start} -> {end})"
                )


def test_returns_list_of_ndarrays(small_mesh):
    """Public API contract: list of np.ndarray with integer dtype."""
    for k in range(small_mesh.n_layers):
        paths = small_mesh.paths_on_outer_vertices(k)
        assert isinstance(paths, list)
        for p in paths:
            assert isinstance(p, np.ndarray)
            assert p.dtype.kind in ("i", "u")
            assert p.ndim == 1


def test_invalid_layer_index_raises(small_mesh):
    """Out-of-range layer index must raise IndexError."""
    with pytest.raises(IndexError):
        small_mesh.paths_on_outer_vertices(small_mesh.n_layers)
    with pytest.raises(IndexError):
        small_mesh.paths_on_outer_vertices(-1)


def test_module_level_and_method_agree(small_mesh):
    """The free function and the bound method must produce equal output."""
    for k in range(small_mesh.n_layers):
        via_method = small_mesh.paths_on_outer_vertices(k)
        via_function = paths_on_outer_vertices(small_mesh, k)
        assert len(via_method) == len(via_function)
        for a, b in zip(via_method, via_function):
            np.testing.assert_array_equal(a, b)


def test_donut_layer_zero_has_two_disjoint_cycles():
    """Topology check: donut layer 0 must produce two disconnected paths
    (the outer perimeter and the inner hole boundary)."""
    mesh = examples.donut()
    paths = mesh.paths_on_outer_vertices(0)
    assert len(paths) >= 2, (
        f"donut layer 0 should have >=2 paths (outer + inner ring), got {len(paths)}"
    )
    # No vertex should appear in two different paths if the graph really
    # has two disjoint components (no junctions at L0 for donut).
    seen_verts = set()
    for p in paths:
        verts_in_path = {int(v) for v in p}
        overlap = seen_verts & verts_in_path
        assert not overlap, (
            f"donut L0: paths share vertex {overlap} despite expected "
            "disjoint cycles"
        )
        seen_verts |= verts_in_path


def _adjacency_from_graph(mesh, layer_idx):
    """Helper: rebuild adjacency for invariants (not relying on internals)."""
    ov = mesh.layers["OV"][layer_idx]
    _, graph_edges = _layer_graph(mesh, layer_idx)
    e2v = mesh.adjacencies["Edge2Vert"]
    adj: dict[int, list[int]] = {int(v): [] for v in ov}
    for eid in graph_edges:
        u, v = int(e2v[eid, 0]), int(e2v[eid, 1])
        adj[u].append(v)
        adj[v].append(u)
    return adj
