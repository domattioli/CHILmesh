"""Cross-fixture invariants required by AUDIT_REPORT.md Phase 1.

Replaces the value-pinned annulus assertions with structural invariants
that hold for any well-formed CHILmesh.
"""
from __future__ import annotations

from itertools import chain

import numpy as np
import pytest

import chilmesh

FIXTURES = ["annulus", "donut", "block_o", "structured"]


@pytest.fixture(scope="module", params=FIXTURES)
def mesh(request):
    return getattr(chilmesh.examples, request.param)()


def test_connectivity_wellformed(mesh):
    """No duplicate vertex per element, all indices in range, no zero-area."""
    conn = mesh.connectivity_list
    assert (conn >= 0).all(), (
        f"connectivity has negative vertex IDs: min={conn.min()}"
    )
    assert (conn < mesh.n_verts).all(), (
        f"connectivity references vertex IDs >= n_verts ({mesh.n_verts}): "
        f"max={conn.max()}"
    )

    if mesh.type == "Triangular":
        rows = conn[:, :3]
        # No duplicates within a row
        for row in rows:
            assert len(set(row.tolist())) == 3, f"Degenerate row {row}"
    # Positive area
    areas = mesh.signed_area()
    assert (areas > 0).all(), (
        f"min signed area {areas.min()} is non-positive; "
        f"{int((areas <= 0).sum())} of {areas.size} elements affected"
    )


def test_layers_disjoint_cover(mesh):
    """Every element appears in exactly one (OE or IE) entry across layers."""
    seen = []
    for oe in mesh.layers["OE"]:
        seen.extend(int(e) for e in oe)
    for ie in mesh.layers["IE"]:
        seen.extend(int(e) for e in ie)

    assert len(seen) == len(set(seen)), "Duplicate element across layers"
    assert set(seen) == set(range(mesh.n_elems)), (
        f"Layer cover misses {set(range(mesh.n_elems)) - set(seen)} or "
        f"includes spurious ids"
    )


def test_layers_at_least_one(mesh):
    assert mesh.n_layers >= 1, f"mesh has no layers (n_layers={mesh.n_layers})"
    assert len(mesh.layers["OE"]) == mesh.n_layers, (
        f"layers['OE'] length {len(mesh.layers['OE'])} != n_layers {mesh.n_layers}"
    )


def test_layers_outer_layer_touches_boundary(mesh):
    """The first layer's bEdgeIDs must equal mesh.boundary_edges()."""
    first_b = set(int(e) for e in mesh.layers["bEdgeIDs"][0])
    actual_b = set(int(e) for e in mesh.boundary_edges())
    assert first_b == actual_b, (
        f"First-layer bEdgeIDs {len(first_b)} != boundary_edges() {len(actual_b)}"
    )


def test_interior_angles_no_nan(mesh):
    angles = mesh.interior_angles()
    nan_count = int(np.isnan(angles).sum())
    assert nan_count == 0, (
        f"interior_angles produced {nan_count} NaN entries out of {angles.size}"
    )
