"""Mesh topology utilities for efficient edge lookup and management."""

from typing import Dict, Tuple, Optional
import numpy as np


class EdgeMap:
    """Hash-based edge lookup for O(1) edge ID queries.

    Stores edges in canonical form (v1, v2) where v1 < v2 to ensure
    consistent lookups regardless of vertex order.
    """

    def __init__(self) -> None:
        """Initialize empty edge map."""
        self._map: Dict[Tuple[int, int], int] = {}
        self._next_id: int = 0

    def add_edge(self, v1: int, v2: int) -> int:
        """Add edge and return its ID.

        Edges are stored in canonical form (min, max). If edge already exists,
        returns its existing ID (idempotent).

        Parameters:
            v1, v2: Vertex indices

        Returns:
            Unique edge ID (0-indexed)

        Time Complexity: O(1) amortized
        """
        edge = (min(v1, v2), max(v1, v2))
        if edge not in self._map:
            self._map[edge] = self._next_id
            self._next_id += 1
        return self._map[edge]

    def find_edge(self, v1: int, v2: int) -> Optional[int]:
        """Find edge ID, or None if not found.

        Parameters:
            v1, v2: Vertex indices

        Returns:
            Edge ID if exists, None otherwise

        Time Complexity: O(1) amortized
        """
        edge = (min(v1, v2), max(v1, v2))
        return self._map.get(edge)

    def remove_edge(self, v1: int, v2: int) -> bool:
        """Remove edge, return True if existed.

        Parameters:
            v1, v2: Vertex indices

        Returns:
            True if edge was present and removed, False otherwise

        Time Complexity: O(1) amortized
        """
        edge = (min(v1, v2), max(v1, v2))
        if edge in self._map:
            del self._map[edge]
            return True
        return False

    def __len__(self) -> int:
        """Return total number of edges."""
        return len(self._map)

    def __contains__(self, edge: Tuple[int, int]) -> bool:
        """Check if edge exists (in canonical form).

        Parameters:
            edge: Tuple of (v1, v2)

        Returns:
            True if edge exists
        """
        canonical = (min(edge[0], edge[1]), max(edge[0], edge[1]))
        return canonical in self._map

    def to_array(self) -> np.ndarray:
        """Convert to Edge2Vert ndarray format.

        Returns:
            ndarray of shape (n_edges, 2) where each row is (v1, v2)
            with v1 < v2. Rows are sorted by edge ID.

        Time Complexity: O(n log n) due to sorting
        """
        if len(self._map) == 0:
            return np.array([], dtype=int).reshape(0, 2)

        # Sort by edge ID to get consistent ordering
        sorted_edges = sorted(self._map.items(), key=lambda x: x[1])
        return np.array([edge for edge, _ in sorted_edges], dtype=int)

    def to_list(self) -> list:
        """Convert to list of edges.

        Returns:
            List of (v1, v2) tuples sorted by edge ID

        Time Complexity: O(n log n) due to sorting
        """
        if len(self._map) == 0:
            return []

        sorted_edges = sorted(self._map.items(), key=lambda x: x[1])
        return [edge for edge, _ in sorted_edges]
