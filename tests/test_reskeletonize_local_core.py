"""Tests for reskeletonize_local core path: partial re-peel on deep element mutations.

This test module exercises the UNTESTED core block of reskeletonize_local
(mutations.py lines 604-671), which only runs when start_layer > 0 (i.e.,
affected element in layer >= radius+1). Existing test_mutations.py tests only
exercise the FALLBACK paths (full _skeletonize rebuild) because they mutate
shallow layer 0/1 elements.

See: mutations.py::MutableMesh.reskeletonize_local (lines 552-671)
"""

from __future__ import annotations

import pytest
import numpy as np
from chilmesh import CHILmesh, MutableMesh
from chilmesh import examples


@pytest.fixture(params=["annulus", "structured"])
def deep_mesh_name(request):
    """Parametrize over fixtures deep enough to have layer index >= 3."""
    return request.param


def _fresh(name: str) -> CHILmesh:
    """Load a fresh (uncached) mesh by name."""
    example_fn = getattr(examples, name)
    if hasattr(example_fn, '__wrapped__'):
        return example_fn.__wrapped__()
    return example_fn()


class TestReskeletonizeLocalCore:
    """Tests for reskeletonize_local core partial re-peel path."""

    def test_core_path_runs_deep_element(self, deep_mesh_name):
        """Core block executes when deep element is re-peeled with radius=2.

        Verifies precondition that the selected element forces entry into the
        core partial-peel logic (start_layer > 0), then confirms layer structure
        is valid after re-peel.
        """
        mesh = _fresh(deep_mesh_name)
        mutable = MutableMesh(mesh)

        # Require depth: need n_layers >= 4 to guarantee layer index >= 3
        if mesh.n_layers < 4:
            pytest.skip(f"{deep_mesh_name} has only {mesh.n_layers} layers, need >= 4")

        # Pick a deep element: the first from the deepest layer
        deepest_layer_idx = mesh.n_layers - 1
        deep_elem = int(mesh.layers['OE'][deepest_layer_idx][0])

        # Compute affected_layer by scanning OE/IE
        elem_set = {deep_elem}
        affected_layer = None
        for iL in range(mesh.n_layers):
            oe_set = set(int(e) for e in mesh.layers['OE'][iL])
            ie_set = set(int(e) for e in mesh.layers['IE'][iL])
            if elem_set & (oe_set | ie_set):
                affected_layer = iL
                break

        assert affected_layer is not None, "Element not found in any layer"
        assert affected_layer >= 3, f"Element in layer {affected_layer}, need >= 3"

        # Verify precondition: with radius=2, start_layer = affected_layer - 2
        start_layer = max(0, affected_layer - 2)
        assert start_layer > 0, f"start_layer would be {start_layer}, core block requires > 0"

        # Run reskeletonize_local with radius=2
        mutable.reskeletonize_local(np.array([deep_elem]), radius=2)

        # Assert layer structure is valid
        assert mesh.n_layers >= 4, f"After re-peel, n_layers = {mesh.n_layers}, expected >= 4"
        assert len(mesh.layers['OE']) == mesh.n_layers
        assert len(mesh.layers['IE']) == mesh.n_layers
        assert len(mesh.layers['OV']) == mesh.n_layers
        assert len(mesh.layers['IV']) == mesh.n_layers
        assert len(mesh.layers['bEdgeIDs']) == mesh.n_layers

    def test_parity_with_full_skeletonize_radius2(self, deep_mesh_name):
        """Partial re-peel (radius=2) produces layers identical to full rebuild.

        This is the gold assertion: given unchanged mesh topology, a partial
        re-peel from a deep checkpoint must reproduce full _skeletonize() output.
        """
        # Load two independent fresh meshes from the same name
        mesh_local = _fresh(deep_mesh_name)
        mesh_full = _fresh(deep_mesh_name)

        if mesh_local.n_layers < 4:
            pytest.skip(f"{deep_mesh_name} has only {mesh_local.n_layers} layers, need >= 4")

        # Pick deepest element on mesh_local
        deep_elem = int(mesh_local.layers['OE'][mesh_local.n_layers - 1][0])

        # Verify it's in a deep layer (allows start_layer > 0)
        elem_set = {deep_elem}
        affected_layer = None
        for iL in range(mesh_local.n_layers):
            oe_set = set(int(e) for e in mesh_local.layers['OE'][iL])
            ie_set = set(int(e) for e in mesh_local.layers['IE'][iL])
            if elem_set & (oe_set | ie_set):
                affected_layer = iL
                break

        if affected_layer is None or affected_layer < 3:
            pytest.skip(f"Element in layer {affected_layer}, need >= 3 for core path")

        # On mesh_local: run reskeletonize_local
        mutable_local = MutableMesh(mesh_local)
        mutable_local.reskeletonize_local(np.array([deep_elem]), radius=2)

        # Capture OE, IE, OV, IV as per-layer sets on mesh_local
        local_oe = [set(int(e) for e in mesh_local.layers['OE'][iL])
                    for iL in range(mesh_local.n_layers)]
        local_ie = [set(int(e) for e in mesh_local.layers['IE'][iL])
                    for iL in range(mesh_local.n_layers)]
        local_ov = [set(int(e) for e in mesh_local.layers['OV'][iL])
                    for iL in range(mesh_local.n_layers)]
        local_iv = [set(int(e) for e in mesh_local.layers['IV'][iL])
                    for iL in range(mesh_local.n_layers)]

        # On mesh_full: run full _skeletonize
        mesh_full._skeletonize()

        # Capture the same on mesh_full
        full_oe = [set(int(e) for e in mesh_full.layers['OE'][iL])
                   for iL in range(mesh_full.n_layers)]
        full_ie = [set(int(e) for e in mesh_full.layers['IE'][iL])
                   for iL in range(mesh_full.n_layers)]
        full_ov = [set(int(e) for e in mesh_full.layers['OV'][iL])
                   for iL in range(mesh_full.n_layers)]
        full_iv = [set(int(e) for e in mesh_full.layers['IV'][iL])
                   for iL in range(mesh_full.n_layers)]

        # Assert n_layers equal
        assert mesh_local.n_layers == mesh_full.n_layers, \
            f"n_layers differ: partial={mesh_local.n_layers}, full={mesh_full.n_layers}"

        # Assert per-layer sets are identical for all four layer types
        assert local_oe == full_oe, "OE per-layer sets differ"
        assert local_ie == full_ie, "IE per-layer sets differ"
        assert local_ov == full_ov, "OV per-layer sets differ"
        assert local_iv == full_iv, "IV per-layer sets differ"

    def test_parity_with_full_skeletonize_radius1(self, deep_mesh_name):
        """Partial re-peel (radius=1) produces layers identical to full rebuild.

        Same parity check as radius=2 test, but with smaller safety margin.
        Allows testing both radius values within core-path preconditions.
        """
        mesh_local = _fresh(deep_mesh_name)
        mesh_full = _fresh(deep_mesh_name)

        if mesh_local.n_layers < 4:
            pytest.skip(f"{deep_mesh_name} has only {mesh_local.n_layers} layers, need >= 4")

        deep_elem = int(mesh_local.layers['OE'][mesh_local.n_layers - 1][0])

        # Verify it's in a deep layer (allows start_layer > 0 with radius=1)
        elem_set = {deep_elem}
        affected_layer = None
        for iL in range(mesh_local.n_layers):
            oe_set = set(int(e) for e in mesh_local.layers['OE'][iL])
            ie_set = set(int(e) for e in mesh_local.layers['IE'][iL])
            if elem_set & (oe_set | ie_set):
                affected_layer = iL
                break

        if affected_layer is None or affected_layer < 2:
            pytest.skip(f"Element in layer {affected_layer}, need >= 2 for core path with radius=1")

        # On mesh_local: run reskeletonize_local with radius=1
        mutable_local = MutableMesh(mesh_local)
        mutable_local.reskeletonize_local(np.array([deep_elem]), radius=1)

        # Capture OE per-layer sets
        local_oe = [set(int(e) for e in mesh_local.layers['OE'][iL])
                    for iL in range(mesh_local.n_layers)]

        # On mesh_full: run full _skeletonize
        mesh_full._skeletonize()
        full_oe = [set(int(e) for e in mesh_full.layers['OE'][iL])
                   for iL in range(mesh_full.n_layers)]

        # Assert n_layers and OE parity
        assert mesh_local.n_layers == mesh_full.n_layers, \
            f"n_layers differ: partial={mesh_local.n_layers}, full={mesh_full.n_layers}"
        assert local_oe == full_oe, "OE per-layer sets differ"

    def test_no_negative_ids_after_local_repeel(self, deep_mesh_name):
        """After deep radius=2 re-peel, all layer element/vertex IDs are non-negative.

        Negative IDs indicate incomplete re-peeling or data corruption.
        Also verifies that OV and IV in each layer are disjoint.
        """
        mesh = _fresh(deep_mesh_name)
        mutable = MutableMesh(mesh)

        if mesh.n_layers < 4:
            pytest.skip(f"{deep_mesh_name} has only {mesh.n_layers} layers, need >= 4")

        deep_elem = int(mesh.layers['OE'][mesh.n_layers - 1][0])

        # Verify deep precondition
        elem_set = {deep_elem}
        affected_layer = None
        for iL in range(mesh.n_layers):
            oe_set = set(int(e) for e in mesh.layers['OE'][iL])
            ie_set = set(int(e) for e in mesh.layers['IE'][iL])
            if elem_set & (oe_set | ie_set):
                affected_layer = iL
                break

        if affected_layer is None or affected_layer < 3:
            pytest.skip(f"Element in layer {affected_layer}, need >= 3")

        mutable.reskeletonize_local(np.array([deep_elem]), radius=2)

        # Assert all IDs are non-negative
        for iL in range(mesh.n_layers):
            oe = mesh.layers['OE'][iL]
            ie = mesh.layers['IE'][iL]
            ov = mesh.layers['OV'][iL]
            iv = mesh.layers['IV'][iL]

            assert np.all(oe >= 0), f"Layer {iL} OE contains negative: {oe[oe < 0]}"
            assert np.all(ie >= 0), f"Layer {iL} IE contains negative: {ie[ie < 0]}"
            assert np.all(ov >= 0), f"Layer {iL} OV contains negative: {ov[ov < 0]}"
            assert np.all(iv >= 0), f"Layer {iL} IV contains negative: {iv[iv < 0]}"

            # Assert OV and IV are disjoint
            ov_set = set(int(v) for v in ov)
            iv_set = set(int(v) for v in iv)
            overlap = ov_set & iv_set
            assert len(overlap) == 0, \
                f"Layer {iL} OV and IV overlap at vertices: {overlap}"
