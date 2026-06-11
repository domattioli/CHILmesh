"""Mesh topology utilities for efficient edge lookup and management."""
from __future__ import annotations


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


def _polygon_signed_area(pts: np.ndarray) -> float:
    """Shoelace signed area of an (n, 2) polygon; positive when CCW."""
    x = pts[:, 0]
    y = pts[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))


def quad_from_tri_pair(points: np.ndarray, tri_a: np.ndarray, tri_b: np.ndarray) -> np.ndarray:
    """Pure (non-mutating) quad vertex ordering from two adjacent triangles.

    Given two adjacent triangles sharing exactly one edge, construct the quad
    formed by removing the shared edge (interior diagonal) and retaining the
    four distinct vertices in CCW order.

    Parameters
    ----------
    points : np.ndarray
        (n, 2) or (n, 3) array of vertex coordinates. Only first 2 cols used.
    tri_a, tri_b : np.ndarray
        Integer vertex index arrays, shape (3,) or (4,) with padding.
        Padded arrays: [v0, v1, v2, v2] or with -1 sentinel for unused slots.
        Must each contain exactly 3 unique valid (>= 0) vertices when deduplicated.

    Returns
    -------
    np.ndarray
        Shape (4,) integer array [apex_a, s0, apex_b, s1] in CCW order,
        where (s0, s1) are the shared-edge vertices (sorted) and apex_a, apex_b
        are the opposite corners. Quad area computed via shoelace formula is > 0.

    Raises
    ------
    ValueError
        If tri_a and tri_b do not share exactly 2 vertices (non-adjacent or identical).
        If either input is not a valid triangle (< 3 unique valid vertices).

    Notes
    -----
    - Purity: input arrays are not modified; computation is side-effect-free.
    - Padding detection: trailing -1 or duplicate vertices are filtered.
    - Winding: result guaranteed CCW if input vertex coordinates are CCW-ordered.
    """
    # Extract unique valid verts from each triangle (filter -1 and duplicates).
    verts_a = [int(v) for v in tri_a if int(v) >= 0]
    verts_a = sorted(set(verts_a))

    verts_b = [int(v) for v in tri_b if int(v) >= 0]
    verts_b = sorted(set(verts_b))

    # Validate: must be exactly 3 verts each.
    if len(verts_a) != 3:
        raise ValueError(f"tri_a must have exactly 3 unique valid vertices; got {len(verts_a)}")
    if len(verts_b) != 3:
        raise ValueError(f"tri_b must have exactly 3 unique valid vertices; got {len(verts_b)}")

    # Find shared edge: must be exactly 2 common vertices.
    shared = set(verts_a) & set(verts_b)
    if len(shared) != 2:
        raise ValueError(
            f"triangles must share exactly 2 vertices; got {len(shared)}. "
            f"tri_a={verts_a}, tri_b={verts_b}"
        )

    # Opposite corners.
    apex_a = (set(verts_a) - shared).pop()
    apex_b = (set(verts_b) - shared).pop()

    # Shared edge in sorted order.
    s0, s1 = sorted(shared)

    # Build quad: [apex_a, s0, apex_b, s1].
    quad = np.array([apex_a, s0, apex_b, s1], dtype=np.int64)

    # Apply CCW winding correction.
    signed_area = _polygon_signed_area(points[quad, :2])
    if signed_area < 0:
        quad = np.array([apex_a, s1, apex_b, s0], dtype=np.int64)

    return quad


def quads_from_tri_pairs(
    points: np.ndarray, conn: np.ndarray, pairs: np.ndarray
) -> np.ndarray:
    """Batch form of quad_from_tri_pair.

    Parameters
    ----------
    points : np.ndarray
        (n, 2+) vertex coordinates.
    conn : np.ndarray
        Element connectivity array, shape (m, 3) or (m, 4).
        First 3 columns used as triangle vertices.
    pairs : np.ndarray or list
        (n_pairs, 2) array of element-index pairs.

    Returns
    -------
    np.ndarray
        Shape (n_pairs, 4) with one quad per row.

    Raises
    ------
    ValueError
        For any pair that violates quad_from_tri_pair preconditions.
    """
    pairs = np.asarray(pairs)
    if pairs.ndim != 2 or pairs.shape[1] != 2:
        raise ValueError(f"pairs must have shape (n, 2); got {pairs.shape}")

    quads = []
    for elem_a, elem_b in pairs:
        quad = quad_from_tri_pair(points, conn[elem_a, :], conn[elem_b, :])
        quads.append(quad)

    if len(quads) == 0:
        return np.array([], dtype=np.int64).reshape(0, 4)
    return np.array(quads, dtype=np.int64)
