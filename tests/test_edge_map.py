"""Unit tests for EdgeMap class."""

import pytest
import numpy as np
from chilmesh.mesh_topology import EdgeMap


class TestEdgeMapBasics:
    """Test basic EdgeMap operations."""

    def test_empty_map(self):
        """Empty map has zero edges."""
        edge_map = EdgeMap()
        assert len(edge_map) == 0
        assert edge_map.find_edge(0, 1) is None

    def test_add_single_edge(self):
        """Can add a single edge."""
        edge_map = EdgeMap()
        edge_id = edge_map.add_edge(0, 1)
        assert edge_id == 0
        assert len(edge_map) == 1
        assert edge_map.find_edge(0, 1) == 0

    def test_add_multiple_edges(self):
        """Can add multiple edges with increasing IDs."""
        edge_map = EdgeMap()
        e0 = edge_map.add_edge(0, 1)
        e1 = edge_map.add_edge(1, 2)
        e2 = edge_map.add_edge(2, 3)
        assert e0 == 0 and e1 == 1 and e2 == 2
        assert len(edge_map) == 3

    def test_canonical_form(self):
        """Edges stored in canonical form (v1 < v2)."""
        edge_map = EdgeMap()
        e1 = edge_map.add_edge(0, 5)
        e2 = edge_map.add_edge(5, 0)  # Reverse order
        assert e1 == e2
        assert len(edge_map) == 1

    def test_idempotent_add(self):
        """Adding same edge twice returns same ID."""
        edge_map = EdgeMap()
        e1 = edge_map.add_edge(2, 7)
        e2 = edge_map.add_edge(2, 7)
        assert e1 == e2
        assert len(edge_map) == 1


class TestEdgeMapFind:
    """Test edge finding operations."""

    def test_find_existing(self):
        """find_edge returns ID for existing edge."""
        edge_map = EdgeMap()
        added_id = edge_map.add_edge(1, 3)
        found_id = edge_map.find_edge(1, 3)
        assert found_id == added_id

    def test_find_nonexistent(self):
        """find_edge returns None for nonexistent edge."""
        edge_map = EdgeMap()
        edge_map.add_edge(0, 1)
        assert edge_map.find_edge(0, 2) is None
        assert edge_map.find_edge(5, 10) is None

    def test_find_respects_canonical_form(self):
        """find_edge finds edges regardless of vertex order."""
        edge_map = EdgeMap()
        edge_map.add_edge(3, 7)
        assert edge_map.find_edge(3, 7) is not None
        assert edge_map.find_edge(7, 3) is not None


class TestEdgeMapRemove:
    """Test edge removal operations."""

    def test_remove_existing(self):
        """remove_edge deletes existing edge."""
        edge_map = EdgeMap()
        edge_map.add_edge(0, 1)
        result = edge_map.remove_edge(0, 1)
        assert result is True
        assert len(edge_map) == 0
        assert edge_map.find_edge(0, 1) is None

    def test_remove_nonexistent(self):
        """remove_edge returns False for nonexistent edge."""
        edge_map = EdgeMap()
        result = edge_map.remove_edge(0, 1)
        assert result is False
        assert len(edge_map) == 0

    def test_remove_respects_canonical_form(self):
        """remove_edge deletes edge regardless of vertex order."""
        edge_map = EdgeMap()
        edge_map.add_edge(2, 5)
        result = edge_map.remove_edge(5, 2)  # Reverse order
        assert result is True
        assert len(edge_map) == 0


class TestEdgeMapContains:
    """Test membership checking."""

    def test_contains_existing(self):
        """__contains__ returns True for existing edges."""
        edge_map = EdgeMap()
        edge_map.add_edge(1, 4)
        assert (1, 4) in edge_map
        assert (4, 1) in edge_map

    def test_contains_nonexistent(self):
        """__contains__ returns False for nonexistent edges."""
        edge_map = EdgeMap()
        assert (0, 1) not in edge_map


class TestEdgeMapConversion:
    """Test conversion to array and list formats."""

    def test_to_array_empty(self):
        """to_array on empty map returns shape (0, 2)."""
        edge_map = EdgeMap()
        arr = edge_map.to_array()
        assert arr.shape == (0, 2)
        assert arr.dtype == int

    def test_to_array_single_edge(self):
        """to_array converts single edge correctly."""
        edge_map = EdgeMap()
        edge_map.add_edge(0, 3)
        arr = edge_map.to_array()
        assert arr.shape == (1, 2)
        assert arr[0, 0] == 0 and arr[0, 1] == 3

    def test_to_array_multiple_edges(self):
        """to_array preserves edge ordering by ID."""
        edge_map = EdgeMap()
        edge_map.add_edge(0, 1)  # ID 0
        edge_map.add_edge(2, 3)  # ID 1
        edge_map.add_edge(4, 5)  # ID 2
        arr = edge_map.to_array()
        assert arr.shape == (3, 2)
        assert arr[0, 0] == 0 and arr[0, 1] == 1
        assert arr[1, 0] == 2 and arr[1, 1] == 3
        assert arr[2, 0] == 4 and arr[2, 1] == 5

    def test_to_array_canonical_form(self):
        """to_array stores edges in canonical form (v1 < v2)."""
        edge_map = EdgeMap()
        edge_map.add_edge(5, 1)  # Added as (1, 5)
        arr = edge_map.to_array()
        assert arr[0, 0] < arr[0, 1]

    def test_to_list_empty(self):
        """to_list on empty map returns empty list."""
        edge_map = EdgeMap()
        lst = edge_map.to_list()
        assert lst == []

    def test_to_list_multiple_edges(self):
        """to_list preserves edge ordering by ID."""
        edge_map = EdgeMap()
        edge_map.add_edge(0, 1)
        edge_map.add_edge(2, 3)
        lst = edge_map.to_list()
        assert len(lst) == 2
        assert lst[0] == (0, 1)
        assert lst[1] == (2, 3)


class TestEdgeMapIntegration:
    """Integration tests simulating mesh edge discovery."""

    def test_triangle_edges(self):
        """EdgeMap handles triangle edges correctly."""
        edge_map = EdgeMap()
        # Triangle with vertices 0, 1, 2
        edges = [(0, 1), (1, 2), (2, 0)]
        for v1, v2 in edges:
            edge_map.add_edge(v1, v2)
        assert len(edge_map) == 3
        for v1, v2 in edges:
            assert edge_map.find_edge(v1, v2) is not None

    def test_quad_edges(self):
        """EdgeMap handles quad edges correctly."""
        edge_map = EdgeMap()
        # Quad with vertices 0, 1, 2, 3
        edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
        for v1, v2 in edges:
            edge_map.add_edge(v1, v2)
        assert len(edge_map) == 4

    def test_shared_edges(self):
        """EdgeMap deduplicates shared edges between elements."""
        edge_map = EdgeMap()
        # Two triangles sharing edge (0, 1)
        # Triangle 1: (0, 1, 2)
        edge_map.add_edge(0, 1)
        edge_map.add_edge(1, 2)
        edge_map.add_edge(2, 0)
        # Triangle 2: (0, 1, 3) - shares edge (0, 1)
        edge_map.add_edge(0, 1)  # Already exists
        edge_map.add_edge(1, 3)
        edge_map.add_edge(3, 0)
        assert len(edge_map) == 5  # 0-1, 1-2, 2-0, 1-3, 3-0

    def test_large_mesh(self):
        """EdgeMap handles large mesh (100 vertices, structured grid)."""
        edge_map = EdgeMap()
        n = 10  # 10x10 grid
        # Add edges for a structured grid
        for i in range(n):
            for j in range(n):
                v = i * n + j
                # Horizontal edges
                if j < n - 1:
                    edge_map.add_edge(v, v + 1)
                # Vertical edges
                if i < n - 1:
                    edge_map.add_edge(v, v + n)
        # 10x9 horizontal + 9x10 vertical = 180 edges
        assert len(edge_map) == 180
