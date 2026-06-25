"""Unit tests for the two-tier junction heuristic in layer_paths._pick_next and _walk dead-end.

Covers:
  - _pick_next early returns (lines 248-249)
  - _pick_next heuristic 1: same-face preference (lines 252-256)
  - _pick_next heuristic 2: turning angle minimization (lines 258-273)
  - _walk dead-end break (line 226)
"""
from __future__ import annotations

import numpy as np
import pytest

from chilmesh.layer_paths import _pick_next, _walk


class TestPickNextEarlyReturn:
    """Test _pick_next early returns: single outgoing or no prev_edge."""

    def test_pick_next_single_outgoing_returns_it(self):
        """Single outgoing edge → return it regardless of prev_edge."""
        outgoing = [(2, 10)]  # (neighbor_id, edge_id)
        edge_layer_elems = {10: {5}}
        points = np.array([[0, 0], [1, 0], [2, 0]], dtype=float)

        result = _pick_next(
            cur=1,
            prev_v=0,
            prev_edge=None,
            outgoing=outgoing,
            edge_layer_elems=edge_layer_elems,
            points=points,
        )

        assert result == (2, 10), f"Expected (2, 10), got {result}"

    def test_pick_next_prev_edge_none_returns_first(self):
        """Multiple outgoing but prev_edge is None → return outgoing[0]."""
        outgoing = [(2, 10), (3, 11), (4, 12)]
        edge_layer_elems = {10: {5}, 11: {5}, 12: {9}}
        points = np.array(
            [[0, 0], [1, 0], [2, 0], [3, 0], [4, 0]], dtype=float
        )

        result = _pick_next(
            cur=1,
            prev_v=None,
            prev_edge=None,  # Critical: None triggers early return
            outgoing=outgoing,
            edge_layer_elems=edge_layer_elems,
            points=points,
        )

        assert result == (2, 10), f"Expected first candidate (2, 10), got {result}"

    def test_pick_next_single_outgoing_ignores_geometry(self):
        """Single outgoing → early return before any geometry checks."""
        # Even with wild geometry (zero norm, degenerate), single returns it.
        outgoing = [(2, 100)]
        edge_layer_elems = {100: {999}}
        points = np.array([[0, 0], [0, 0], [0, 0]], dtype=float)  # Degenerate

        result = _pick_next(
            cur=1,
            prev_v=0,
            prev_edge=99,
            outgoing=outgoing,
            edge_layer_elems=edge_layer_elems,
            points=points,
        )

        assert result == (2, 100), "Single outgoing must return immediately"


class TestPickNextHeuristic1SameFace:
    """Test heuristic 1: prefer edges sharing a layer-element with prev_edge."""

    def test_heuristic1_unique_same_face_returns_it(self):
        """Two candidates, exactly one shares a face → return that one."""
        # prev_edge=10 has element {5}
        # Candidate A (nbr=2, eid=11) shares element {5}
        # Candidate B (nbr=3, eid=12) has disjoint element {9}
        # Expected: return A
        prev_edge = 10
        outgoing = [(2, 11), (3, 12)]
        edge_layer_elems = {
            10: {5},      # prev_edge
            11: {5},      # same_face candidate
            12: {9},      # disjoint
        }
        points = np.array([[0, 0], [1, 0], [2, 0], [3, 0]], dtype=float)

        result = _pick_next(
            cur=1,
            prev_v=0,
            prev_edge=prev_edge,
            outgoing=outgoing,
            edge_layer_elems=edge_layer_elems,
            points=points,
        )

        # Should return the one with same_face (nbr=2, eid=11)
        assert result == (2, 11), f"Expected (2, 11), got {result}"

    def test_heuristic1_multiple_same_face_falls_through(self):
        """Multiple candidates share face → falls through to heuristic 2 (turning angle)."""
        # prev_edge=10 has element {5}
        # Both candidates share {5}, so same_face has 2 elements
        # len(candidates) > 1 → fall through to heuristic 2 (turning angle)
        prev_edge = 10
        outgoing = [(2, 11), (3, 12)]
        edge_layer_elems = {
            10: {5},  # prev_edge
            11: {5},  # shares
            12: {5},  # also shares
        }
        # Geometry: prev_v at (0) = (-1, 0), cur at (1) = (0, 0), so incoming direction a=(1, 0)
        # Neighbor 2 at (1, 0): turning angle = 0 (straight)
        # Neighbor 3 at (0, 1): turning angle = π/2 (90° turn)
        # Expect neighbor 2 (smaller angle)
        points = np.array(
            [[-1, 0], [0, 0], [1, 0], [0, 1]], dtype=float
        )

        result = _pick_next(
            cur=1,
            prev_v=0,
            prev_edge=prev_edge,
            outgoing=outgoing,
            edge_layer_elems=edge_layer_elems,
            points=points,
        )

        # Heuristic 2 picks min turning angle → (2, 11) with angle ~0
        assert result == (2, 11), f"Expected straight path (2, 11), got {result}"

    def test_heuristic1_none_match_fallback_to_all(self):
        """No candidates share a face → candidates = all outgoing."""
        # prev_edge=10 has element {5}
        # Neither candidate shares {5}: 11 has {7}, 12 has {9}
        # same_face is empty → candidates = outgoing (all 2)
        # Fall through to heuristic 2
        prev_edge = 10
        outgoing = [(2, 11), (3, 12)]
        edge_layer_elems = {
            10: {5},
            11: {7},  # disjoint
            12: {9},  # disjoint
        }
        points = np.array(
            [[-1, 0], [0, 0], [1, 0], [0, 1]], dtype=float
        )

        result = _pick_next(
            cur=1,
            prev_v=0,
            prev_edge=prev_edge,
            outgoing=outgoing,
            edge_layer_elems=edge_layer_elems,
            points=points,
        )

        # Heuristic 2: (2, 11) at (1, 0) has angle ~0, (3, 12) at (0, 1) has angle ~π/2
        # Expect (2, 11)
        assert result == (2, 11), f"Expected smallest angle (2, 11), got {result}"


class TestPickNextHeuristic2TurningAngle:
    """Test heuristic 2: minimize turning angle when >1 candidates."""

    def test_heuristic2_straight_wins_over_turn(self):
        """Two candidates with same-face: straight path beats 90° turn."""
        prev_edge = 10
        outgoing = [(2, 11), (3, 12)]
        edge_layer_elems = {
            10: {5},
            11: {5},  # both share element {5}
            12: {5},
        }
        # Geometry: prev_v=(0)=(-1,0), cur=(1)=(0,0) → incoming dir a=(1,0)
        # Neighbor 2 at (1,0) → outgoing dir b=(1,0) → angle=0 (straight)
        # Neighbor 3 at (0,1) → outgoing dir b=(0,1) → angle=π/2 (90° turn)
        points = np.array(
            [[-1, 0], [0, 0], [1, 0], [0, 1]], dtype=float
        )

        result = _pick_next(
            cur=1,
            prev_v=0,
            prev_edge=prev_edge,
            outgoing=outgoing,
            edge_layer_elems=edge_layer_elems,
            points=points,
        )

        assert result == (2, 11), f"Expected straight (2, 11), got {result}"

    def test_heuristic2_left_turn_vs_right_turn(self):
        """Two candidates, one left-turn 45°, one right-turn 45°."""
        # Both have same turning magnitude but opposite directions.
        # cos_ang for left-turn 45°: cos(π/4) ≈ 0.707
        # cos_ang for right-turn 45°: cos(π/4) ≈ 0.707 (same magnitude)
        # Both have the same angle, so ties go to min which picks the first.
        # (In Python, ties on min pick the first in iteration order.)
        prev_edge = 10
        # Candidate 2: 45° left turn
        # Candidate 3: 45° right turn (same magnitude, opposite direction)
        # Both at distance ~sqrt(2) from cur to make angles exact.
        outgoing = [(2, 11), (3, 12)]
        edge_layer_elems = {
            10: {5},
            11: {5},
            12: {5},
        }
        # prev_v=(0)=(-1,0), cur=(1)=(0,0), a=(1,0)
        # Neighbor 2 at (0.707, 0.707): left-turn 45°
        # Neighbor 3 at (0.707, -0.707): right-turn 45° (symmetric)
        # Both should have angle ≈ π/4, so min picks the first: (2, 11)
        points = np.array(
            [[-1, 0], [0, 0], [0.707, 0.707], [0.707, -0.707]], dtype=float
        )

        result = _pick_next(
            cur=1,
            prev_v=0,
            prev_edge=prev_edge,
            outgoing=outgoing,
            edge_layer_elems=edge_layer_elems,
            points=points,
        )

        # Ties go to min, which picks the first: (2, 11)
        assert result == (2, 11), f"Expected first of ties (2, 11), got {result}"

    def test_heuristic2_acute_beats_obtuse(self):
        """Acute turning angle beats obtuse (more-than-90°) angle."""
        prev_edge = 10
        outgoing = [(2, 11), (3, 12)]
        edge_layer_elems = {10: {5}, 11: {5}, 12: {5}}
        # prev_v=(0)=(-1,0), cur=(1)=(0,0), a=(1,0)
        # Neighbor 2 at (1, 0.5): acute angle
        # Neighbor 3 at (-1, 0): 180° U-turn (obtuse)
        points = np.array([[-1, 0], [0, 0], [1, 0.5], [-1, 0]], dtype=float)

        result = _pick_next(
            cur=1,
            prev_v=0,
            prev_edge=prev_edge,
            outgoing=outgoing,
            edge_layer_elems=edge_layer_elems,
            points=points,
        )

        assert result == (2, 11), f"Expected acute angle (2, 11), got {result}"

    def test_heuristic2_gentle_left_vs_sharp_right(self):
        """Gentle 30° left-turn beats sharp 89° right-turn."""
        prev_edge = 10
        outgoing = [(2, 11), (3, 12)]
        edge_layer_elems = {10: {5}, 11: {5}, 12: {5}}
        # prev_v=(0)=(-1,0), cur=(1)=(0,0), a=(1,0)
        # Neighbor 2: gentle 30° left-turn → cos_ang=cos(30°)≈0.866
        # Neighbor 3: sharp 89° right-turn → cos_ang=cos(89°)≈0.0175
        # Expect 2 (smaller angle)
        points = np.array(
            [[-1, 0], [0, 0], [0.866, 0.5], [0.0175, -0.9998]], dtype=float
        )

        result = _pick_next(
            cur=1,
            prev_v=0,
            prev_edge=prev_edge,
            outgoing=outgoing,
            edge_layer_elems=edge_layer_elems,
            points=points,
        )

        assert result == (2, 11), f"Expected gentle turn (2, 11), got {result}"


class TestWalkDeadEnd:
    """Test _walk dead-end: when walk runs out of unvisited edges."""

    def test_walk_dead_end_open_path_three_vertex_line(self):
        """Three vertices in a line: 0-1-2, start at middle (1), walk hits dead-end."""
        # Graph: 0 --(eid=100)-- 1 --(eid=101)-- 2
        # Start at 1 (middle), no visited edges.
        # First step: from 1, outgoing edges are (0, 100) and (2, 101).
        # _pick_next with prev_edge=None → returns outgoing[0] = (0, 100)
        # Walk to 0, mark eid 100 as visited.
        # From 0, only unvisited edge is... none (100 is visited, and no other edge).
        # break at line 226 (dead-end).
        # Path should be [1, 0].

        adj = {
            0: [(1, 100)],
            1: [(0, 100), (2, 101)],
            2: [(1, 101)],
        }
        edge_layer_elems = {100: {5}, 101: {6}}
        visited: set[int] = set()
        points = np.array([[0, 0], [1, 0], [2, 0]], dtype=float)

        path = _walk(
            start=1,
            adj=adj,
            edge_layer_elems=edge_layer_elems,
            visited=visited,
            points=points,
        )

        # Should walk from 1 to 0 and stop (dead-end).
        expected_path = [1, 0]
        np.testing.assert_array_equal(
            path, expected_path,
            err_msg=f"Expected path {expected_path}, got {path}"
        )
        # Verify only eid 100 was visited.
        assert visited == {100}, f"Expected visited={{100}}, got {visited}"

    def test_walk_closed_cycle_four_vertices(self):
        """Four vertices in a square cycle: 0-1-2-3-0, walk closes."""
        # Graph (square): 0 --(eid=10)-- 1 --(eid=11)-- 2
        #                 |                              |
        #               eid=13                        eid=12
        #                 |                              |
        #                 3 ----------(eid=14)----------
        # Start at 0, walk all edges, close back to 0.
        # outgoing[0] order matters; construct adj so first is deterministic.
        adj = {
            0: [(1, 10), (3, 13)],
            1: [(0, 10), (2, 11)],
            2: [(1, 11), (3, 12)],
            3: [(2, 12), (0, 13)],
        }
        edge_layer_elems = {10: {5}, 11: {5}, 12: {5}, 13: {5}}
        visited: set[int] = set()
        points = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], dtype=float)

        path = _walk(
            start=0,
            adj=adj,
            edge_layer_elems=edge_layer_elems,
            visited=visited,
            points=points,
        )

        # Closed walk: [0, 1, 2, 3, 0]
        expected_path = [0, 1, 2, 3, 0]
        np.testing.assert_array_equal(path, expected_path)
        # All 4 edges visited.
        assert visited == {10, 11, 12, 13}

    def test_walk_dead_end_two_vertex_edge(self):
        """Simplest case: single edge 0-1, start at 1."""
        adj = {
            0: [(1, 10)],
            1: [(0, 10)],
        }
        edge_layer_elems = {10: {5}}
        visited: set[int] = set()
        points = np.array([[0, 0], [1, 0]], dtype=float)

        path = _walk(
            start=1,
            adj=adj,
            edge_layer_elems=edge_layer_elems,
            visited=visited,
            points=points,
        )

        # Walk: 1 → 0 → dead-end (no unvisited outgoing from 0).
        expected_path = [1, 0]
        np.testing.assert_array_equal(path, expected_path)
        assert visited == {10}

    def test_walk_dead_end_y_shaped_graph(self):
        """Y-shaped graph: 0-1-2, 1-3. Start at leaf 2."""
        # Graph:     2
        #            |
        #          (11)
        #            |
        #      0 --(10)-- 1 --(12)-- 3
        # Start at leaf 2. Walk: 2 → 1. From 1, pick next among (0, 10) and (3, 12).
        # With prev_edge=11 (edge from 2 to 1), heuristic 1 may or may not prefer one.
        # For simplicity, construct edges so prev_edge=11 with element {5}:
        #   edge 10: {5}, edge 11: {5}, edge 12: {7}
        # Then heuristic 1 prefers 0 (via edge 10, shares {5}).
        # Walk: 2 → 1 → 0 → dead-end.
        adj = {
            0: [(1, 10)],
            1: [(0, 10), (2, 11), (3, 12)],
            2: [(1, 11)],
            3: [(1, 12)],
        }
        edge_layer_elems = {10: {5}, 11: {5}, 12: {7}}
        visited: set[int] = set()
        # Colinear points: 2, 1, 0 on x-axis; 3 off-axis.
        points = np.array([[0, 0], [1, 0], [2, 0], [1, 1]], dtype=float)

        path = _walk(
            start=2,
            adj=adj,
            edge_layer_elems=edge_layer_elems,
            visited=visited,
            points=points,
        )

        # Walk: 2 → 1 → 0 (heuristic 1 prefers edge 10 which shares {5}).
        # Then 0 has no unvisited outgoing (only edge 10, already visited).
        # Dead-end break.
        expected_path = [2, 1, 0]
        np.testing.assert_array_equal(path, expected_path)
        assert visited == {11, 10}


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
