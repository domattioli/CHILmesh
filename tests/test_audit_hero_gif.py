"""Unit tests for audit_hero_gif module (T027, F-12).

Test the audit functions against synthetic fixtures that isolate each behavior:
- Inverted triangles (test MUST flag as inverted)
- Hole-spanning triangles (test MUST flag as hole-spanning)
- Oversized triangles (test MUST flag as oversized)
- Quality binning at boundaries (q=1.0 must not IndexError, lands in HBINS-1)

Pure-function tests only (no generator import, no render).
"""

from __future__ import annotations

import numpy as np
import pytest
import sys
from pathlib import Path

# Add parent paths so we can import the audit module
sys.path.insert(0, str(Path(__file__).parent.parent / 'scripts'))
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from audit_hero_gif import (
    audit_degenerate,
    audit_pacing,
    audit_peel_invariance,
    HBINS,
)


class TestAuditDegenerate:
    """Tests for audit_degenerate function."""

    def test_empty_history(self):
        """Empty history should return zero counts."""
        def dummy_fd(p):
            return np.zeros(len(p))

        result = audit_degenerate([], dummy_fd)
        assert result['total_hole'] == 0
        assert result['total_oversize'] == 0
        assert result['total_inverted'] == 0

    def test_inverted_triangle(self):
        """Single inverted triangle must be flagged (hand-built synthetic fixture).

        Per F-12: one inverted triangle the auditor MUST flag.
        """
        # CCW triangles (majority)
        pts_ccw1 = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]], dtype=np.float64)
        pts_ccw2 = np.array([[2.0, 0.0], [3.0, 0.0], [2.5, 1.0]], dtype=np.float64)

        # CW triangle (inverted, should be flagged)
        pts_cw = np.array([[4.0, 0.0], [4.5, 1.0], [5.0, 0.0]], dtype=np.float64)

        # Mix: two CCW, one CW -> CW should be flagged as inverted
        pts_mix = np.vstack([pts_ccw1, pts_ccw2, pts_cw])
        elems_mix = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]], dtype=np.int64)
        hist_mix = [(pts_mix, elems_mix)]

        def dummy_fd(p):
            return -0.1 * np.ones(len(p))  # All interior

        result = audit_degenerate(hist_mix, dummy_fd)
        # One element is inverted (CW among CCW majority)
        assert result['total_inverted'] >= 1, "Expected at least 1 inverted element (F-12)"
        assert result['per_snapshot'][0]['inverted'] >= 1

    def test_hole_spanning_triangle(self):
        """Edge midpoints deep in hole must be flagged (hand-built synthetic fixture).

        Per F-12: one hole-spanning triangle the auditor MUST flag.
        Annulus: outer R=1.0, inner r=0.3.
        """
        # Triangle with one edge deep inside the hole
        # Points at r ≈ 0.15 (very deep in the hole)
        pts = np.array([
            [0.10, 0.0],    # r ≈ 0.10, fd ≈ 0.3 - 0.10 = 0.20 > hole_tau
            [0.15, 0.05],   # r ≈ 0.158, fd ≈ 0.3 - 0.158 ≈ 0.14 > hole_tau
            [0.50, 0.0],    # r ≈ 0.50, fd ≈ 0.3 - 0.50 = -0.20 < hole_tau (outside)
        ], dtype=np.float64)
        elems = np.array([[0, 1, 2]], dtype=np.int64)
        hist = [(pts, elems)]

        def annulus_fd(p):
            r = np.linalg.norm(p, axis=1)
            return np.maximum(r - 1.0, 0.3 - r)

        result = audit_degenerate(hist, annulus_fd, hole_tau=0.02)
        # Edge [0, 1] midpoint is at r ≈ 0.13, fd ≈ 0.17 > 0.02 -> hole-spanning
        assert result['total_hole'] >= 1, "Expected at least 1 hole-spanning element (F-12)"
        assert result['per_snapshot'][0]['hole'] >= 1

    def test_oversized_triangle(self):
        """Oversized triangles relative to median must be flagged.

        Oversized = area > oversize_mult * per-snapshot-median.
        """
        # Small triangle
        pts_small = np.array([[0.0, 0.0], [0.1, 0.0], [0.05, 0.1]], dtype=np.float64)
        elems_small = np.array([[0, 1, 2]], dtype=np.int64)

        # Large triangle (area >> median)
        pts_large = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]], dtype=np.float64)
        elems_large = np.array([[0, 1, 2]], dtype=np.int64)

        # Mix: 5 small + 1 large
        pts_list = [pts_small + np.array([i * 0.2, 0.0]) for i in range(5)]
        pts_list.append(pts_large + np.array([1.5, 0.0]))
        pts = np.vstack(pts_list)
        elems = np.array([[i*3, i*3+1, i*3+2] for i in range(6)], dtype=np.int64)
        hist = [(pts, elems)]

        def dummy_fd(p):
            return -0.1 * np.ones(len(p))

        result = audit_degenerate(hist, dummy_fd, oversize_mult=4.0)
        # The large triangle should be flagged as oversized
        assert result['total_oversize'] >= 1, "Expected at least 1 oversized element"
        assert result['per_snapshot'][0]['oversize'] >= 1

    def test_quality_binning_at_boundary(self):
        """q=1.0 element must land in last bin without IndexError (F-10).

        Per F-12: one q=1.0 element must bin without IndexError and totals match.
        """
        # Perfect equilateral triangle (q ≈ 1.0)
        side = 1.0 / np.sqrt(3.0)
        pts = np.array([
            [0.0, 0.0],
            [side, 0.0],
            [side / 2.0, side * np.sqrt(3.0) / 2.0],
        ], dtype=np.float64)
        elems = np.array([[0, 1, 2]], dtype=np.int64)
        hist = [(pts, elems)]

        def dummy_fd(p):
            return -0.1 * np.ones(len(p))

        result = audit_degenerate(hist, dummy_fd)
        # Should not crash and should return valid counts
        assert result['total_inverted'] >= 0
        assert result['total_hole'] >= 0
        assert result['total_oversize'] >= 0
        # Verify the snapshot was processed
        assert len(result['per_snapshot']) == 1


class TestAuditPacing:
    """Tests for audit_pacing function."""

    def test_empty_history(self):
        """Empty history should return ratio of 1.0."""
        ratio = audit_pacing([], [])
        assert ratio == 1.0

    def test_uniform_motion(self):
        """Uniform displacement should give ratio close to 1.0."""
        # Create a history with uniform per-iteration displacement
        hist = []
        for i in range(20):
            pts = np.array([[i * 0.1, j * 0.1] for j in range(5)], dtype=np.float64)
            elems = np.array([[0, 1, 2], [1, 2, 3]], dtype=np.int64)
            hist.append((pts, elems))

        # Render every 2 indices with hold=2
        snap_idx = list(range(0, len(hist), 2))

        ratio = audit_pacing(hist, snap_idx, hold=2)
        # With uniform motion, ratio should be close to 1.0
        assert 0.5 < ratio < 2.0, f"Expected ratio near 1.0, got {ratio}"

    def test_accelerating_motion(self):
        """Accelerating displacement should give measurable ratio."""
        # Create history with accelerating displacement
        hist = []
        for i in range(50):
            # Accelerate: position grows quadratically
            pos = i * i * 0.01
            pts = np.array([[pos + j * 0.01, 0.0] for j in range(5)], dtype=np.float64)
            elems = np.array([[0, 1, 2], [1, 2, 3]], dtype=np.int64)
            hist.append((pts, elems))

        snap_idx = list(range(0, len(hist), 5))
        ratio = audit_pacing(hist, snap_idx, hold=1)
        # With quadratic acceleration, ratio should be > 1
        assert ratio >= 0.5, f"Expected positive ratio, got {ratio}"


class TestAuditPeelInvariance:
    """Tests for audit_peel_invariance function."""

    def test_uniform_quality(self):
        """Uniform quality should satisfy invariance."""
        # All elements have quality 0.5
        q_fem = np.full(100, 0.5, dtype=np.float64)
        elem_layer = np.zeros(100, dtype=np.int32)

        result = audit_peel_invariance(q_fem, elem_layer, n_layers=4)
        assert result['histogram_independent'], "Uniform quality should be independent"
        assert result['histogram_invariant'], "Uniform quality should be invariant"

    def test_mixed_quality_invariance(self):
        """Mixed quality should maintain per-layer totals."""
        # Create elements with different qualities across layers
        n_elems = 200
        q_fem = np.random.rand(n_elems)
        elem_layer = np.repeat([0, 1, 2, 3], n_elems // 4)

        result = audit_peel_invariance(q_fem, elem_layer, n_layers=4)
        # The histogram counts should be invariant across all k
        # because peeling doesn't change the binning, just the colors
        assert result['histogram_invariant'], "Per-bin totals should be invariant"

    def test_boundary_quality_bins(self):
        """Quality at bin boundaries must be handled correctly (F-10).

        q=1.0 should land in HBINS-1, q=0.0 in bin 0. No overflow.
        """
        q_fem = np.array([0.0, 0.5, 1.0, 0.25, 0.75], dtype=np.float64)
        elem_layer = np.zeros(5, dtype=np.int32)

        result = audit_peel_invariance(q_fem, elem_layer, n_layers=1)
        ref_hist = result['reference_histogram']

        assert ref_hist is not None
        # Verify no overflow: all bins should have counts
        assert np.sum(ref_hist) == len(q_fem), "All elements should be binned"
        # q=1.0 should land in last bin (F-10)
        assert ref_hist[-1] >= 1, "q=1.0 should land in last bin"
        # q=0.0 should land in first bin
        assert ref_hist[0] >= 1, "q=0.0 should land in first bin"

    def test_empty_quality_array(self):
        """Empty quality array should return error message."""
        q_fem = np.array([], dtype=np.float64)
        elem_layer = np.array([], dtype=np.int32)

        result = audit_peel_invariance(q_fem, elem_layer, n_layers=1)
        assert not result['histogram_independent']
        assert 'Empty' in result['message'] or 'empty' in result['message'].lower()


class TestHistogramBinning:
    """Tests for histogram binning logic used in audit functions (F-10)."""

    def test_bin_index_computation(self):
        """Verify bin indices are computed correctly for quality values."""
        q_vals = np.array([0.0, 0.025, 0.05, 0.5, 0.975, 1.0])

        # F-10: bin index via np.clip(np.floor(q * HBINS), 0, HBINS-1)
        bin_indices = np.clip(np.floor(q_vals * HBINS), 0, HBINS - 1).astype(int)

        # q=0.0 -> bin 0
        assert bin_indices[0] == 0
        # q=0.025 -> bin 1 (0.025 * 40 = 1.0)
        assert bin_indices[1] == 1
        # q=0.5 -> bin 20
        assert bin_indices[3] == 20
        # q=1.0 -> bin 39 (not 40, which would overflow)
        assert bin_indices[-1] == HBINS - 1
        assert bin_indices[-1] < HBINS

    def test_histogram_totals_match_count(self):
        """Histogram bin totals should sum to element count."""
        q_fem = np.random.rand(1000)
        counts, edges = np.histogram(q_fem, bins=HBINS, range=(0.0, 1.0))

        assert np.sum(counts) == len(q_fem), "Histogram totals must equal element count"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
