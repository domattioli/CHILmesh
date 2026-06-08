"""Tests for boundary-type seeding in skeletonization (#129).

Verifies that ``_skeletonize(seed_boundary_kinds=..., seed_ibtypes=...)`` limits
layer-0 peeling to edges whose nodes belong to matching boundary segments, while
the default (no filter) remains backward-compatible with existing MATLAB-parity
behaviour.
"""
from __future__ import annotations

import warnings
from copy import deepcopy

import numpy as np
import pytest

from chilmesh import examples, CHILmesh


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _structured_with_segments():
    """Return structured mesh with synthetic boundary_segments and the node sets.

    Two segments:
    - ``flow`` (IBTYPE=0): nodes on bottom edge (y ≈ -1)
    - ``open`` (IBTYPE=None): remaining boundary nodes
    """
    m = deepcopy(examples.structured())
    bottom = np.where(np.isclose(m.points[:, 1], -1.0))[0]
    other = np.array(
        [v for v in m.layers["OV"][0] if v not in set(bottom.tolist())], dtype=int
    )
    m.boundary_segments = [
        {"kind": "flow", "ibtype": 0, "nodes": bottom},
        {"kind": "open", "ibtype": None, "nodes": other},
    ]
    return m, bottom, other


# ---------------------------------------------------------------------------
# Default behaviour (backward compatibility)
# ---------------------------------------------------------------------------

class TestDefaultBackwardCompat:
    def test_default_seeds_all_boundary_nodes(self):
        """seed_boundary_kinds=None must reproduce the current OV[0]."""
        m = examples.structured()
        default_ov0 = set(m.layers["OV"][0].tolist())

        m2 = deepcopy(examples.structured())
        m2._skeletonize()  # explicit call with defaults
        assert set(m2.layers["OV"][0].tolist()) == default_ov0

    def test_default_n_layers_unchanged(self):
        m = deepcopy(examples.structured())
        n = m.n_layers
        m._skeletonize()
        assert m.n_layers == n


# ---------------------------------------------------------------------------
# Kind-based filtering
# ---------------------------------------------------------------------------

class TestKindFilter:
    def test_flow_only_ov0_subset_of_bottom(self):
        """OV[0] when seeding from 'flow' must be a subset of bottom-edge nodes."""
        m, bottom, _ = _structured_with_segments()
        bottom_set = set(bottom.tolist())
        m._skeletonize(seed_boundary_kinds=["flow"])
        assert set(m.layers["OV"][0].tolist()).issubset(bottom_set)

    def test_open_only_ov0_subset_of_other(self):
        """OV[0] when seeding from 'open' must be a subset of non-bottom nodes."""
        m, _, other = _structured_with_segments()
        other_set = set(other.tolist())
        m._skeletonize(seed_boundary_kinds=["open"])
        assert set(m.layers["OV"][0].tolist()).issubset(other_set)

    def test_both_kinds_equals_default(self):
        """Seeding from both kinds must reproduce the default OV[0]."""
        m_default = examples.structured()
        default_ov0 = set(m_default.layers["OV"][0].tolist())

        m, _, _ = _structured_with_segments()
        m._skeletonize(seed_boundary_kinds=["flow", "open"])
        assert set(m.layers["OV"][0].tolist()) == default_ov0

    def test_filtered_ov0_smaller_than_default(self):
        """Single-kind seed must produce strictly fewer OV[0] nodes than unfiltered."""
        m_default = examples.structured()
        default_size = len(m_default.layers["OV"][0])

        m, _, _ = _structured_with_segments()
        m._skeletonize(seed_boundary_kinds=["flow"])
        assert len(m.layers["OV"][0]) < default_size


# ---------------------------------------------------------------------------
# IBTYPE-based filtering
# ---------------------------------------------------------------------------

class TestIbtypeFilter:
    def test_ibtype0_seeds_flow_nodes(self):
        """seed_ibtypes=[0] must seed from nodes of the IBTYPE=0 (flow) segment."""
        m, bottom, _ = _structured_with_segments()
        bottom_set = set(bottom.tolist())
        m._skeletonize(seed_ibtypes=[0])
        assert set(m.layers["OV"][0].tolist()).issubset(bottom_set)

    def test_ibtype_none_seeds_open_nodes(self):
        """seed_ibtypes=[None] must seed from the open (ibtype=None) segment."""
        m, _, other = _structured_with_segments()
        other_set = set(other.tolist())
        m._skeletonize(seed_ibtypes=[None])
        assert set(m.layers["OV"][0].tolist()).issubset(other_set)


# ---------------------------------------------------------------------------
# Combined kind + ibtype (intersection)
# ---------------------------------------------------------------------------

class TestCombinedFilter:
    def test_intersection_narrows_to_matching_segment(self):
        """kind='flow' AND ibtype=0 must match exactly one segment (flow, ibtype=0)."""
        m, bottom, _ = _structured_with_segments()
        bottom_set = set(bottom.tolist())
        m._skeletonize(seed_boundary_kinds=["flow"], seed_ibtypes=[0])
        assert set(m.layers["OV"][0].tolist()).issubset(bottom_set)

    def test_conflicting_combination_raises(self):
        """kind='flow' AND ibtype=None should raise (no segment has both)."""
        m, _, _ = _structured_with_segments()
        with pytest.raises(ValueError, match="No boundary segments match"):
            m._skeletonize(seed_boundary_kinds=["flow"], seed_ibtypes=[None])


# ---------------------------------------------------------------------------
# Error and fallback paths
# ---------------------------------------------------------------------------

class TestErrorAndFallback:
    def test_no_matching_kind_raises(self):
        """Filtering to a kind that no segment has must raise ValueError."""
        m, _, _ = _structured_with_segments()
        with pytest.raises(ValueError, match="No boundary segments match"):
            m._skeletonize(seed_boundary_kinds=["nonexistent"])

    def test_no_matching_ibtype_raises(self):
        m, _, _ = _structured_with_segments()
        with pytest.raises(ValueError, match="No boundary segments match"):
            m._skeletonize(seed_ibtypes=[999])

    def test_empty_boundary_segments_warns_and_uses_all(self):
        """When boundary_segments is empty, fall back to all edges with a warning."""
        m = deepcopy(examples.structured())
        # explicitly empty (default for node+elem-only meshes)
        m.boundary_segments = []
        default_ov0 = set(m.layers["OV"][0].tolist())

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            m._skeletonize(seed_boundary_kinds=["flow"])

        assert any("boundary_segments is empty" in str(w.message) for w in caught)
        assert set(m.layers["OV"][0].tolist()) == default_ov0

    def test_compute_layers_false_accepts_filter_silently(self):
        """seed_boundary_kinds passed to __init__ with compute_layers=False: no error."""
        m = examples.structured()
        pts, conn = m.points[:, :2], m.connectivity_list
        m2 = CHILmesh(
            connectivity=conn,
            points=pts,
            compute_layers=False,
            seed_boundary_kinds=["flow"],
        )
        assert m2.n_layers == 0


# ---------------------------------------------------------------------------
# Public __init__ passthrough
# ---------------------------------------------------------------------------

class TestPublicAPIPassthrough:
    def test_init_with_seed_kinds_produces_filtered_layers(self):
        """CHILmesh(..., seed_boundary_kinds=['flow']) must apply filter at construction."""
        m_base, bottom, _ = _structured_with_segments()
        pts = m_base.points[:, :2]
        conn = m_base.connectivity_list

        m = CHILmesh(connectivity=conn, points=pts, seed_boundary_kinds=["flow"])
        # Without boundary_segments the filter cannot apply — warn + fallback.
        # This tests the code path, not the filter result.
        assert m.n_layers > 0

    def test_init_with_seed_kinds_and_segments(self):
        """CHILmesh constructed then segments injected; _skeletonize re-run with filter."""
        m, bottom, _ = _structured_with_segments()
        bottom_set = set(bottom.tolist())
        # Re-run skeletonize now that segments are set
        m._skeletonize(seed_boundary_kinds=["flow"])
        ov0 = set(m.layers["OV"][0].tolist())
        assert ov0.issubset(bottom_set)
        assert len(ov0) > 0
