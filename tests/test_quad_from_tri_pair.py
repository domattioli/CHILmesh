"""Tests for pure tri-pair → quad connectivity helpers."""
from __future__ import annotations

import pytest
import numpy as np

import chilmesh
from chilmesh.mesh_topology import quad_from_tri_pair, quads_from_tri_pairs


class TestQuadFromTriPair:
    """Tests for quad_from_tri_pair scalar form."""

    def test_basic_adjacent_tris_ccw(self):
        """Two CCW triangles sharing one edge → returns 4-vert CCW quad."""
        # Simple two-triangle quad in the (x,y) plane.
        points = np.array([
            [0.0, 0.0],  # 0: bottom-left
            [1.0, 0.0],  # 1: bottom-right
            [0.5, 1.0],  # 2: top-left (apex of tri_a)
            [0.5, -1.0],  # 3: bottom (apex of tri_b)
        ])

        # tri_a: [0, 1, 2] CCW
        tri_a = np.array([0, 1, 2])

        # tri_b: [0, 1, 3] CCW; shares edge (0, 1)
        tri_b = np.array([0, 1, 3])

        quad = quad_from_tri_pair(points, tri_a, tri_b)

        # Expected: [apex_a=2, s0=0, apex_b=3, s1=1] or [apex_a=2, s1=1, apex_b=3, s0=0]
        # depending on CCW winding. Both are valid.
        assert quad.shape == (4,)
        assert set(quad) == {0, 1, 2, 3}

        # Verify quad is CCW by checking signed area > 0.
        quad_pts = points[quad]
        x, y = quad_pts[:, 0], quad_pts[:, 1]
        signed_area = 0.5 * (np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))
        assert signed_area > 0, f"Quad {quad} has negative area {signed_area}"

    def test_ccw_correction(self):
        """Two tris that would give CW quad are auto-corrected to CCW."""
        # Triangle pair where the naive [apex_a, s0, apex_b, s1] ordering
        # would be CW; expect auto-correction to [apex_a, s1, apex_b, s0].
        points = np.array([
            [0.0, 0.0],  # 0: bottom-left
            [1.0, 0.0],  # 1: bottom-right
            [0.5, -1.0],  # 2: bottom (apex of tri_a, CCW with shared edge 0-1)
            [0.5, 1.0],  # 3: top (apex of tri_b, CCW with shared edge 0-1)
        ])

        tri_a = np.array([0, 1, 2])
        tri_b = np.array([0, 1, 3])

        quad = quad_from_tri_pair(points, tri_a, tri_b)

        # Verify result is CCW.
        quad_pts = points[quad]
        x, y = quad_pts[:, 0], quad_pts[:, 1]
        signed_area = 0.5 * (np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))
        assert signed_area > 0, f"Quad {quad} has negative area {signed_area}"

    def test_padded_input_arrays(self):
        """Padded [a,b,c,c] or [a,b,c,-1] inputs work correctly."""
        points = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 1.0],
            [0.5, -1.0],
        ])

        # Padded triangle [0, 1, 2, 2]
        tri_a_padded = np.array([0, 1, 2, 2])

        # Padded triangle [0, 1, 3, -1]
        tri_b_padded = np.array([0, 1, 3, -1])

        quad = quad_from_tri_pair(points, tri_a_padded, tri_b_padded)

        assert quad.shape == (4,)
        assert set(quad) == {0, 1, 2, 3}

    def test_error_non_adjacent_zero_shared(self):
        """Non-adjacent triangles (0 shared verts) → ValueError."""
        points = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 1.0],
            [2.0, 0.0],
            [2.5, 1.0],
            [2.0, 1.0],
        ])

        tri_a = np.array([0, 1, 2])
        tri_b = np.array([3, 4, 5])  # No shared vertices

        with pytest.raises(ValueError, match="must share exactly 2 vertices"):
            quad_from_tri_pair(points, tri_a, tri_b)

    def test_error_non_adjacent_one_shared(self):
        """Non-adjacent triangles (1 shared vert) → ValueError."""
        points = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 1.0],
            [1.0, 0.5],
        ])

        tri_a = np.array([0, 1, 2])
        tri_b = np.array([1, 3, 0])  # Wait, this shares two: 0 and 1. Let me fix.
        tri_b = np.array([1, 3, 4])  # Now only shares vertex 1

        # But vertex 4 doesn't exist. Use a simpler setup.
        points = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 1.0],
            [2.0, 0.5],
            [2.0, 1.5],
        ])

        tri_a = np.array([0, 1, 2])
        tri_b = np.array([1, 3, 4])  # Shares only vertex 1

        with pytest.raises(ValueError, match="must share exactly 2 vertices"):
            quad_from_tri_pair(points, tri_a, tri_b)

    def test_error_identical_tris(self):
        """Identical tri_a == tri_b → ValueError (3 shared vertices)."""
        points = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 1.0],
        ])

        tri = np.array([0, 1, 2])

        with pytest.raises(ValueError, match="must share exactly 2 vertices"):
            quad_from_tri_pair(points, tri, tri)

    def test_error_invalid_quad_input(self):
        """Input with <3 unique valid vertices → ValueError."""
        points = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 1.0],
        ])

        # Only 2 unique verts
        tri_invalid = np.array([0, 0, 1])

        tri_valid = np.array([0, 1, 2])

        with pytest.raises(ValueError, match="must have exactly 3 unique valid vertices"):
            quad_from_tri_pair(points, tri_invalid, tri_valid)

    def test_purity(self):
        """Input arrays unchanged (deep copy check before/after)."""
        points = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 1.0],
            [0.5, -1.0],
        ])

        tri_a = np.array([0, 1, 2])
        tri_b = np.array([0, 1, 3])

        points_copy = points.copy()
        tri_a_copy = tri_a.copy()
        tri_b_copy = tri_b.copy()

        quad = quad_from_tri_pair(points, tri_a, tri_b)

        # Verify inputs unchanged.
        np.testing.assert_array_equal(points, points_copy)
        np.testing.assert_array_equal(tri_a, tri_a_copy)
        np.testing.assert_array_equal(tri_b, tri_b_copy)

        assert quad is not points
        assert quad is not tri_a
        assert quad is not tri_b

    def test_parity_with_mutator(self):
        """Build 2-tri mesh, call quad_from_tri_pair, call mesh.mutations.merge_elements."""
        # Create a simple two-triangle mesh.
        points = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 1.0],
            [0.5, -1.0],
        ], dtype=float)

        conn = np.array([
            [0, 1, 2],  # elem 0: tri_a
            [0, 1, 3],  # elem 1: tri_b
        ])

        mesh = chilmesh.CHILmesh(connectivity=conn, points=points)

        # Pure helper: call quad_from_tri_pair directly.
        quad_pure = quad_from_tri_pair(mesh.points, conn[0, :], conn[1, :])

        # Mutator: call merge_elements to transform the same pair.
        mesh_mut = mesh.copy()
        mutator = chilmesh.MutableMesh(mesh_mut)

        # Find the shared edge between elem 0 and 1.
        # Both share vertices 0 and 1. merge_elements finds it automatically.
        mutator.merge_elements(elem_a=0, elem_b=1)

        # Extract the quad from the mutated mesh.
        quad_mutated = mesh_mut.connectivity_list[0, :]

        # Compare: both should represent the same quad (possibly in different
        # vertex orderings due to CCW correction, but the same 4 vertices).
        assert set(quad_pure) == set(quad_mutated[:4])

    def test_with_3d_points(self):
        """Works with (n, 3) points; only first 2 cols used."""
        points = np.array([
            [0.0, 0.0, 0.1],
            [1.0, 0.0, 0.2],
            [0.5, 1.0, 0.3],
            [0.5, -1.0, 0.4],
        ])

        tri_a = np.array([0, 1, 2])
        tri_b = np.array([0, 1, 3])

        quad = quad_from_tri_pair(points, tri_a, tri_b)

        assert quad.shape == (4,)
        assert set(quad) == {0, 1, 2, 3}


class TestQuadsFromTriPairs:
    """Tests for quads_from_tri_pairs batch form."""

    def test_batch_form_single_pair(self):
        """Batch form on single pair matches scalar call."""
        points = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 1.0],
            [0.5, -1.0],
        ])

        conn = np.array([
            [0, 1, 2],
            [0, 1, 3],
        ])

        pairs = np.array([[0, 1]])

        quad_scalar = quad_from_tri_pair(points, conn[0, :], conn[1, :])
        quads_batch = quads_from_tri_pairs(points, conn, pairs)

        assert quads_batch.shape == (1, 4)
        np.testing.assert_array_equal(quads_batch[0, :], quad_scalar)

    def test_batch_form_multiple_pairs(self):
        """Batch form processes multiple pairs correctly."""
        points = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 1.0],
            [0.5, -1.0],
            [1.5, 0.5],
            [1.5, 1.5],
        ])

        conn = np.array([
            [0, 1, 2],  # 0
            [0, 1, 3],  # 1
            [1, 4, 5],  # 2
            [1, 4, 2],  # 3
        ])

        pairs = np.array([[0, 1], [2, 3]])

        quads = quads_from_tri_pairs(points, conn, pairs)

        assert quads.shape == (2, 4)

        # Verify each quad matches scalar call.
        for i, (a, b) in enumerate(pairs):
            quad_expected = quad_from_tri_pair(points, conn[a, :], conn[b, :])
            np.testing.assert_array_equal(quads[i, :], quad_expected)

    def test_batch_form_list_input(self):
        """Batch form accepts list input for pairs."""
        points = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 1.0],
            [0.5, -1.0],
        ])

        conn = np.array([
            [0, 1, 2],
            [0, 1, 3],
        ])

        pairs_list = [[0, 1]]

        quads = quads_from_tri_pairs(points, conn, pairs_list)

        assert quads.shape == (1, 4)

    def test_batch_form_error_wrong_shape(self):
        """Batch form raises ValueError for wrong pairs shape."""
        points = np.array([[0.0, 0.0], [1.0, 0.0]])
        conn = np.array([[0, 1, 2]])

        pairs_wrong = np.array([[0, 1, 2]])  # shape (1, 3), not (n, 2)

        with pytest.raises(ValueError, match="pairs must have shape"):
            quads_from_tri_pairs(points, conn, pairs_wrong)

    def test_batch_form_empty(self):
        """Batch form handles empty pairs array."""
        points = np.array([[0.0, 0.0], [1.0, 0.0]])
        conn = np.array([[0, 1, 2]])

        pairs = np.array([], dtype=int).reshape(0, 2)

        quads = quads_from_tri_pairs(points, conn, pairs)

        assert quads.shape == (0, 4)


class TestQuadFromTriPairFixtures:
    """Parametrized tests over mesh fixtures (annulus, donut, structured)."""

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "structured"])
    def test_extract_adjacent_interior_tris(self, fixture_name):
        """Extract two adjacent interior triangles and call quad_from_tri_pair."""
        mesh = chilmesh.examples.__getattribute__(fixture_name)()

        conn = mesh.connectivity_list

        # Find two adjacent interior triangles (no boundary edges).
        interior_tris = []
        for i in range(conn.shape[0]):
            elem = conn[i, :3]
            verts = sorted(set(int(v) for v in elem if int(v) >= 0))
            if len(verts) == 3:
                interior_tris.append((i, elem, verts))

        assert len(interior_tris) > 0, f"No triangles found in {fixture_name}"

        # Find a pair of adjacent tris (sharing 2 vertices).
        pair_found = False
        for i, (idx_a, tri_a, verts_a) in enumerate(interior_tris):
            for idx_b, tri_b, verts_b in interior_tris[i+1:]:
                shared = set(verts_a) & set(verts_b)
                if len(shared) == 2:
                    # Found adjacent pair.
                    pair_found = True
                    quad = quad_from_tri_pair(mesh.points, tri_a, tri_b)

                    # Verify result.
                    assert quad.shape == (4,)
                    assert len(set(quad)) == 4
                    assert all(0 <= v < mesh.n_verts for v in quad)

                    # Verify CCW winding.
                    quad_pts = mesh.points[quad, :2]
                    x, y = quad_pts[:, 0], quad_pts[:, 1]
                    signed_area = 0.5 * (np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))
                    assert signed_area > 0, f"Quad {quad} in {fixture_name} has negative area"

                    break
            if pair_found:
                break

        assert pair_found, f"No adjacent tri pairs found in {fixture_name}"

    @pytest.mark.parametrize("fixture_name", ["annulus", "donut", "structured"])
    def test_batch_form_on_multiple_fixture_pairs(self, fixture_name):
        """Call quads_from_tri_pairs on multiple adjacent pairs from fixture."""
        mesh = chilmesh.examples.__getattribute__(fixture_name)()

        conn = mesh.connectivity_list

        # Find multiple adjacent interior triangle pairs.
        pairs = []
        interior_elems = []
        for i in range(min(20, conn.shape[0])):  # Limit to first 20 for speed
            elem = conn[i, :3]
            verts = sorted(set(int(v) for v in elem if int(v) >= 0))
            if len(verts) == 3:
                interior_elems.append((i, verts))

        for i, (idx_a, verts_a) in enumerate(interior_elems):
            for idx_b, verts_b in interior_elems[i+1:]:
                shared = set(verts_a) & set(verts_b)
                if len(shared) == 2 and len(pairs) < 3:  # Collect up to 3 pairs
                    pairs.append([idx_a, idx_b])

        if len(pairs) > 0:
            pairs_arr = np.array(pairs)
            quads = quads_from_tri_pairs(mesh.points, conn, pairs_arr)

            assert quads.shape == (len(pairs), 4)

            # Verify each quad.
            for i, (a, b) in enumerate(pairs):
                quad_expected = quad_from_tri_pair(mesh.points, conn[a, :], conn[b, :])
                np.testing.assert_array_equal(quads[i, :], quad_expected)

                # Verify CCW winding.
                quad_pts = mesh.points[quads[i], :2]
                x, y = quad_pts[:, 0], quad_pts[:, 1]
                signed_area = 0.5 * (np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))
                assert signed_area > 0, f"Quad {i} in {fixture_name} has negative area"
