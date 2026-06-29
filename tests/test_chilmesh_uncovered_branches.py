"""Tests for currently-untested branches in CHILmesh.py.

Covers three specific uncovered branches:
1. CHILmesh.boundary_edges() fallback path (lines 711-721, no Edge2Elem)
2. CHILmesh.element_quality() angle metric branches (lines 1299, 1327-1332)
3. Module-level _check_fort14 success and except branches (lines 2782-2787)
"""

from __future__ import annotations

import numpy as np
import pytest
from pathlib import Path

from chilmesh import CHILmesh, examples
from chilmesh.CHILmesh import _check_fort14


class TestBoundaryEdgesFallbackPath:
    """Tests for boundary_edges() fallback when Edge2Elem is missing.

    Target lines 711-721: The no-Edge2Elem branch counts edge occurrences
    to identify boundary edges (those appearing exactly once).
    """

    def test_boundary_edges_without_edge2elem(self):
        """Lines 711-721: Fallback path when adjacencies lack Edge2Elem.

        Constructs a simple two-triangle mesh (sharing hypotenuse) without
        building adjacencies, forcing the fallback boundary-edge detection.
        """
        # Simple two right-isoceles triangles: [0,1,2] and [1,3,2] sharing edge (1,2)
        pts = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]], dtype=float)
        conn = np.array([[0, 1, 2], [1, 3, 2]], dtype=int)

        # Create mesh WITHOUT building adjacencies
        mesh = CHILmesh(
            connectivity=conn,
            points=pts,
            compute_layers=False,
            compute_adjacencies=False,
        )

        # Verify Edge2Elem is not present
        assert "Edge2Elem" not in mesh.adjacencies, (
            "Expected Edge2Elem missing from adjacencies"
        )

        # Call boundary_edges() which should trigger the fallback
        boundary_ids = mesh.boundary_edges()

        # Two right triangles sharing hypotenuse (1,2):
        # Triangle 1 [0,1,2]: edges (0,1), (1,2), (0,2) → indices 0, 1, 2
        # Triangle 2 [1,3,2]: edges (1,3), (2,3), (1,2) → indices 3, 4, (1,2 seen)
        # edge_list: [(0,1), (1,2), (0,2), (1,3), (2,3)]
        # edge (1,2) appears 2x (interior), rest appear 1x (boundary)
        # Boundary edge indices: 0, 2, 3, 4
        boundary_ids_sorted = sorted(boundary_ids)
        expected = [0, 2, 3, 4]

        np.testing.assert_array_equal(
            boundary_ids_sorted, expected, err_msg="Boundary edge IDs mismatch"
        )
        assert len(boundary_ids) == 4, f"Expected 4 boundary edges, got {len(boundary_ids)}"

    def test_boundary_edges_single_triangle(self):
        """Lines 711-721: Fallback path with single triangle (all edges boundary).

        Isolate case: single triangle with no internal edges.
        """
        pts = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=float)
        conn = np.array([[0, 1, 2]], dtype=int)

        mesh = CHILmesh(
            connectivity=conn,
            points=pts,
            compute_layers=False,
            compute_adjacencies=False,
        )

        boundary_ids = mesh.boundary_edges()

        # All 3 edges are boundary
        assert len(boundary_ids) == 3, f"Expected 3 boundary edges, got {len(boundary_ids)}"
        np.testing.assert_array_equal(sorted(boundary_ids), [0, 1, 2])


class TestElementQualityAngleBranches:
    """Tests for element_quality() angle metric branches.

    Target lines 1299 (ValueError check), 1327-1332 (min_angle vs max_angle).
    Right-isoceles triangles have angles 45°, 45°, 90° (π/4, π/4, π/2).
    """

    def test_element_quality_min_angle(self):
        """Lines 1327-1332: metric='min_angle' branch.

        For right-isoceles triangle: min angle = 45° = π/4 ≈ 0.7853981633974483.
        """
        pts = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]], dtype=float)
        conn = np.array([[0, 1, 2], [1, 3, 2]], dtype=int)

        mesh = CHILmesh(connectivity=conn, points=pts)

        qualities = mesh.element_quality(metric="min_angle")

        # Both triangles are right-isoceles: min angle = 45° = π/4
        expected_min_angle = np.pi / 4
        np.testing.assert_allclose(
            qualities,
            expected_min_angle,
            rtol=1e-6,
            err_msg="min_angle quality mismatch",
        )

    def test_element_quality_max_angle(self):
        """Lines 1327-1332: metric='max_angle' branch.

        For right-isoceles triangle: max angle = 90° = π/2 ≈ 1.5707963267948966.
        """
        pts = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]], dtype=float)
        conn = np.array([[0, 1, 2], [1, 3, 2]], dtype=int)

        mesh = CHILmesh(connectivity=conn, points=pts)

        qualities = mesh.element_quality(metric="max_angle")

        # Both triangles are right-isoceles: max angle = 90° = π/2
        expected_max_angle = np.pi / 2
        np.testing.assert_allclose(
            qualities,
            expected_max_angle,
            rtol=1e-6,
            err_msg="max_angle quality mismatch",
        )

    def test_element_quality_unknown_metric_raises_valueerror(self):
        """Lines 1298-1299: ValueError for unknown metric.

        Tests the metric validation branch.
        """
        pts = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=float)
        conn = np.array([[0, 1, 2]], dtype=int)

        mesh = CHILmesh(connectivity=conn, points=pts)

        # Should raise ValueError for unknown metric
        with pytest.raises(ValueError, match="Unknown metric"):
            mesh.element_quality(metric="bogus")

    def test_element_quality_aspect_ratio_valid(self):
        """Sanity check: aspect_ratio metric still works (not a fallback branch).

        Ensures we didn't break the primary metric.
        """
        pts = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=float)
        conn = np.array([[0, 1, 2]], dtype=int)

        mesh = CHILmesh(connectivity=conn, points=pts)

        qualities = mesh.element_quality(metric="aspect_ratio")

        # aspect_ratio should be > 0 and <= 1
        assert len(qualities) == 1, "Expected 1 element"
        assert 0 < qualities[0] <= 1, f"aspect_ratio out of range: {qualities[0]}"


class TestCheckFort14Branches:
    """Tests for module-level _check_fort14() function.

    Target lines 2782-2787: Both success (return True) and except (return False)
    branches when reading a fort.14 file.
    """

    def test_check_fort14_valid_file(self):
        """Lines 2782-2784: Success branch when fort.14 is readable.

        Uses a bundled fixture that is known to be valid.
        """
        # Use bundled quad_2x2.fort.14
        valid_path = examples.fixture_path("quad_2x2.fort.14")

        result = _check_fort14(str(valid_path))

        assert result is True, f"Expected _check_fort14 to return True for valid file: {valid_path}"

    def test_check_fort14_nonexistent_file(self):
        """Lines 2785-2787: Except branch when file does not exist.

        Nonexistent path should trigger exception, print error, return False.
        """
        nonexistent_path = "/nonexistent/path/does_not_exist.14"

        result = _check_fort14(nonexistent_path)

        assert result is False, "Expected _check_fort14 to return False for nonexistent file"

    def test_check_fort14_invalid_fort14_content(self):
        """Lines 2785-2787: Except branch when file exists but is invalid fort.14.

        Create a temporary file with garbage content (not valid fort.14 format).
        """
        import tempfile
        import os

        # Create a temporary file with invalid content
        with tempfile.NamedTemporaryFile(mode="w", suffix=".14", delete=False) as f:
            f.write("This is not a valid fort.14 file\n")
            temp_path = f.name

        try:
            result = _check_fort14(temp_path)
            assert result is False, "Expected _check_fort14 to return False for invalid fort.14"
        finally:
            # Clean up
            os.unlink(temp_path)
