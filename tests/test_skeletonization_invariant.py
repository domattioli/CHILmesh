"""Regression test for issue #74: skeletonization layer separation invariant.

A vertex appearing in any element of layer k MUST NOT appear in any element
of layer m where |k - m| >= 2. This is the medial-axis layer separation
property that the original MATLAB QuADMesh+ meshLayers function upholds.

This test parametrizes over all four built-in fixtures (annulus, donut,
structured, block_o). It expects 0 violations after the spec 006 fix.
"""
from __future__ import annotations

from typing import Set

import numpy as np
import pytest

from chilmesh import examples


FIXTURES = ["annulus", "donut", "structured", "block_o"]


def _layer_vertices(mesh, layer_idx: int) -> Set[int]:
    """Return the set of vertex IDs in OE[k] ∪ IE[k]."""
    elems = np.concatenate(
        (mesh.layers["OE"][layer_idx], mesh.layers["IE"][layer_idx])
    )
    verts: Set[int] = set()
    for e in elems:
        for v in mesh.connectivity_list[e]:
            vi = int(v)
            if vi >= 0:
                verts.add(vi)
    return verts


@pytest.mark.parametrize("fixture_name", FIXTURES)
def test_layer_separation_invariant(fixture_name: str) -> None:
    """A vertex in layer k MUST NOT appear in any layer m where |k - m| >= 2."""
    mesh = getattr(examples, fixture_name)()
    n = mesh.n_layers

    if n < 3:
        pytest.skip(f"{fixture_name} has only {n} layers; invariant trivially holds")

    layer_verts = [_layer_vertices(mesh, k) for k in range(n)]

    violations = []
    for k in range(n):
        for m in range(k + 2, n):
            shared = layer_verts[k] & layer_verts[m]
            if shared:
                violations.append((k, m, len(shared), sorted(shared)[:5]))

    if violations:
        msg_lines = [f"Layer separation violations in {fixture_name} ({n} layers):"]
        for k, m, count, sample in violations[:10]:
            msg_lines.append(
                f"  Layer {k} <-> Layer {m}: {count} shared vertices "
                f"(sample: {sample})"
            )
        pytest.fail("\n".join(msg_lines))


@pytest.mark.parametrize("fixture_name", FIXTURES)
def test_layer_assignment_partition(fixture_name: str) -> None:
    """Every element MUST be assigned to exactly one layer (no duplicates, no gaps)."""
    mesh = getattr(examples, fixture_name)()

    all_assignments = []
    for k in range(mesh.n_layers):
        all_assignments.extend(int(e) for e in mesh.layers["OE"][k])
        all_assignments.extend(int(e) for e in mesh.layers["IE"][k])

    assigned = set(all_assignments)
    expected = set(range(mesh.n_elems))

    missing = expected - assigned
    duplicates_count = len(all_assignments) - len(assigned)

    assert not missing, (
        f"{fixture_name}: {len(missing)} elements unassigned "
        f"(sample: {sorted(missing)[:5]})"
    )
    assert duplicates_count == 0, (
        f"{fixture_name}: {duplicates_count} elements assigned to multiple layers"
    )
