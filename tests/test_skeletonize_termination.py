"""
Regression test for infinite-loop bug in _skeletonize() (#129 fix).

The bug: broken per-element loop selected boundary edges via `edge2elem_work[:, 1] == -1`,
which only checks slot 1. After an outer ring was consumed, newly-exposed boundary edges
became `[-1, b]` (live element in slot 0), so they were never selected, while fully-
consumed `[-1, -1]` edges WERE selected but yielded zero OE. Nothing got consumed, and
the `while np.any(edge2elem_work >= 0)` loop hung forever.

The fix: replaced the broken per-element loop with the correct vectorized algorithm
that uses `np.sum(edge2elem_work >= 0, axis=1)` to count active adjacencies and
selects edges with `active_count == 1` for layers > 0.
"""

import chilmesh


def test_skeletonize_terminates_on_donut():
    """
    Test that _skeletonize() terminates quickly on the bundled donut fixture.

    Before fix: hangs forever.
    After fix: completes in <1s with n_layers > 0.
    """
    fixture_path = chilmesh.examples.fixture_path("donut_domain.fort.14")
    mesh = chilmesh.CHILmesh.read_from_fort14(fixture_path, compute_layers=True)

    # Basic checks: skeletonization completed with reasonable layer count
    assert mesh.n_layers > 0, "Donut should have multiple layers"
    assert len(mesh.layers["OV"]) == mesh.n_layers
    assert len(mesh.layers["OE"]) == mesh.n_layers
    assert len(mesh.layers["IE"]) == mesh.n_layers
    assert len(mesh.layers["IV"]) == mesh.n_layers
    assert len(mesh.layers["bEdgeIDs"]) == mesh.n_layers
