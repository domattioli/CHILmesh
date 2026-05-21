"""Tests for the public ``CHILmesh.submesh`` API (#138).

The sub-mesh restricts a parent mesh to a subset of elements, remaps vertex
indices to a compact [0, k) numbering, and returns a fully-initialized
``CHILmesh``. Used downstream by quadmesh-matlab's ``two_part_smoother`` to
apply different smoothing to boundary vs interior mesh regions.
"""
from __future__ import annotations

import numpy as np
import pytest

import chilmesh


class TestSubmeshBasics:
    def test_returns_chilmesh(self, mesh):
        elem_ids = np.arange(min(5, mesh.n_elems))
        sub = mesh.submesh(elem_ids)
        assert isinstance(sub, chilmesh.CHILmesh)

    def test_element_count_matches_unique_input(self, mesh):
        elem_ids = np.array([0, 1, 2, 1, 0])  # duplicates intentional
        sub = mesh.submesh(elem_ids)
        assert sub.n_elems == 3

    def test_full_selection_round_trips_topology(self, mesh):
        """Selecting every element produces a mesh with the same n_elems and
        n_verts as the parent (no vertices are dropped because every vertex is
        referenced by at least one element)."""
        sub = mesh.submesh(np.arange(mesh.n_elems))
        assert sub.n_elems == mesh.n_elems
        assert sub.n_verts == mesh.n_verts

    def test_vertex_count_matches_unique_referenced_vertices(self, mesh):
        elem_ids = np.arange(min(10, mesh.n_elems))
        n_unique = np.unique(mesh.connectivity_list[elem_ids]).size
        sub = mesh.submesh(elem_ids)
        assert sub.n_verts == n_unique

    def test_grid_name_suffixed(self, mesh):
        sub = mesh.submesh(np.array([0]))
        if mesh.grid_name is not None:
            assert sub.grid_name == f"{mesh.grid_name}_submesh"
        else:
            assert sub.grid_name is None

    def test_parent_mesh_unchanged(self, fresh_mesh):
        parent_n_elems = fresh_mesh.n_elems
        parent_n_verts = fresh_mesh.n_verts
        parent_conn = fresh_mesh.connectivity_list.copy()
        fresh_mesh.submesh(np.arange(min(5, fresh_mesh.n_elems)))
        assert fresh_mesh.n_elems == parent_n_elems
        assert fresh_mesh.n_verts == parent_n_verts
        np.testing.assert_array_equal(fresh_mesh.connectivity_list, parent_conn)


class TestSubmeshGeometry:
    def test_points_preserved_for_referenced_vertices(self, mesh):
        elem_ids = np.arange(min(5, mesh.n_elems))
        unique_verts = np.unique(mesh.connectivity_list[elem_ids])
        sub = mesh.submesh(elem_ids)
        # Sub-mesh point i must equal parent point unique_verts[i].
        np.testing.assert_allclose(sub.points, mesh.points[unique_verts])

    def test_remapped_connectivity_points_at_correct_vertices(self, mesh):
        """For each element in the sub-mesh, the (x,y,z) of every vertex must
        match the (x,y,z) of the corresponding vertex in the parent element."""
        elem_ids = np.arange(min(5, mesh.n_elems))
        sub = mesh.submesh(elem_ids)
        for sub_elem_idx, parent_elem_id in enumerate(np.sort(elem_ids)):
            parent_row = mesh.connectivity_list[parent_elem_id]
            sub_row = sub.connectivity_list[sub_elem_idx]
            # Geometry must match position-for-position regardless of remapping.
            np.testing.assert_allclose(
                sub.points[sub_row], mesh.points[parent_row]
            )

    def test_signed_areas_match_parent(self, mesh):
        elem_ids = np.arange(min(10, mesh.n_elems))
        parent_areas = mesh.signed_area(elem_ids)
        sub = mesh.submesh(elem_ids)
        sub_areas = sub.signed_area()
        # Ordering: submesh sorts elem_ids, so re-sort parent_areas to match.
        order = np.argsort(elem_ids)
        np.testing.assert_allclose(
            np.sort(sub_areas), np.sort(parent_areas[order])
        )


class TestSubmeshTopology:
    def test_adjacencies_built_by_default(self, mesh):
        sub = mesh.submesh(np.arange(min(5, mesh.n_elems)))
        assert sub.adjacencies != {}
        assert "Edge2Vert" in sub.adjacencies

    def test_compute_layers_false_skips_layers(self, mesh):
        sub = mesh.submesh(
            np.arange(min(5, mesh.n_elems)), compute_layers=False
        )
        assert sub.n_layers == 0
        assert sub.layers["OE"] == []

    def test_compute_adjacencies_false_skips_adjacencies(self, mesh):
        sub = mesh.submesh(
            np.arange(min(5, mesh.n_elems)),
            compute_layers=False,
            compute_adjacencies=False,
        )
        assert sub.adjacencies == {}

    def test_interior_edge_becomes_boundary_when_partner_excluded(self):
        """If element A and B share an edge in the parent but B is excluded
        from the sub-mesh, that edge must be a boundary edge in the sub-mesh.
        """
        # Two adjacent triangles sharing edge (1,2).
        points = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 0.0],
            ]
        )
        connectivity = np.array([[0, 1, 2], [1, 3, 2]], dtype=int)
        parent = chilmesh.CHILmesh(connectivity=connectivity, points=points)

        # Sub-mesh containing only element 0 — the shared edge must now be on
        # the boundary (Edge2Elem second column == -1 for every edge).
        sub = parent.submesh(np.array([0]))
        edge2elem = sub.adjacencies["Edge2Elem"]
        assert np.all(edge2elem[:, 1] == -1), (
            "All edges in a single-triangle submesh must be boundary edges"
        )


class TestSubmeshValidation:
    def test_empty_input_raises(self, mesh):
        with pytest.raises(ValueError, match="at least one element"):
            mesh.submesh(np.array([], dtype=int))

    def test_negative_index_raises(self, mesh):
        with pytest.raises(ValueError, match="out of range"):
            mesh.submesh(np.array([-1]))

    def test_oob_index_raises(self, mesh):
        with pytest.raises(ValueError, match="out of range"):
            mesh.submesh(np.array([mesh.n_elems]))

    def test_accepts_list_input(self, mesh):
        sub = mesh.submesh([0, 1, 2])
        assert sub.n_elems == 3

    def test_accepts_tuple_input(self, mesh):
        sub = mesh.submesh((0, 1, 2))
        assert sub.n_elems == 3


class TestSubmeshLayerExtraction:
    """The primary motivating use case (#138) — extract a skeletonization layer
    as a standalone mesh for differential smoothing."""

    def test_outer_layer_extraction(self, mesh):
        if mesh.n_layers == 0:
            pytest.skip("Mesh has no layers")
        outer_layer_elems = mesh.elements_in_layer(0)
        sub = mesh.submesh(outer_layer_elems)
        assert sub.n_elems == outer_layer_elems.size
        # Vertex count must equal the number of distinct vertices on the layer.
        expected_n_verts = np.unique(
            mesh.connectivity_list[outer_layer_elems]
        ).size
        assert sub.n_verts == expected_n_verts
