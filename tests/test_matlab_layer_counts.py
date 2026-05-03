"""Regression test for issue #74: pinning per-fixture layer counts to MATLAB-correct values.

These counts were captured after the spec 006 fix that ports the MATLAB QuADMesh+
``meshLayers`` algorithm to Python. They are the source of truth for layer counts
going forward.

External ground-truth from MATLAB runs (provided by maintainer):
- Italy domain: 15 layers
- Lake Erie domain: 17 layers
(Italy and Lake Erie meshes are not bundled with CHILmesh; they live in the
ADMESH-Domains catalog.)
"""
from __future__ import annotations

import pytest

from chilmesh import examples


# Captured 2026-05-03 after spec 006 fix landed.
# Each entry is (n_elems, n_layers).
EXPECTED = {
    "annulus": (580, 4),
    "donut": (276, 3),
    "structured": (660, 5),
    "block_o": (5214, 9),
    "quad_2x2": (4, 1),
}


@pytest.mark.parametrize("fixture_name,expected", list(EXPECTED.items()))
def test_layer_count_matches_matlab_reference(fixture_name: str, expected) -> None:
    """Per-fixture layer count must match the MATLAB ``meshLayers`` reference."""
    expected_n_elems, expected_n_layers = expected
    mesh = getattr(examples, fixture_name)()

    assert mesh.n_elems == expected_n_elems, (
        f"{fixture_name}: expected {expected_n_elems} elements, got {mesh.n_elems}"
    )
    assert mesh.n_layers == expected_n_layers, (
        f"{fixture_name}: expected {expected_n_layers} layers, got {mesh.n_layers}. "
        f"If this changes intentionally, update EXPECTED in this file and document why."
    )
