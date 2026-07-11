"""Import-safety smoke test for scripts/generate_hero_animation.py (T005).

Per plan.md Test Plan item 4 / Constitution III ("visualization is
smoke-test only"): this asserts the generator module imports headless
under MPLBACKEND=Agg WITHOUT executing a render, and that the functions
and constants the rest of the feature depends on exist and are sane.
No full render / _stage_data() call happens here — that is a manual
render-gate step (quickstart.md), not a unit test.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import generate_hero_animation as gha


def test_module_imports_headless():
    """Importing the module must not execute a render (import-only cost)."""
    assert gha.OUT_PATH.name == "readme_pipeline_annulus.gif"


def test_key_functions_exist():
    for name in (
        "_stage_data",
        "main",
        "_quality",
        "_boundary_indices",
        "_draw_hist",
        "_draw_metrics",
        "_draw_mesh",
        "_setup_axes",
        "_q_colors",
        "_layer_color",
        "_peel_facecolors",
        "_draw_peel_hist",
        "_optimize_gif",
        "_ease",
    ):
        assert callable(getattr(gha, name, None)), f"missing callable {name}"


def test_pacing_constants_sane():
    assert isinstance(gha.N_SNAP, int) and 8 <= gha.N_SNAP <= 64
    assert isinstance(gha.TRUSS_HOLD, int) and gha.TRUSS_HOLD >= 1
    assert isinstance(gha.HBINS, int) and gha.HBINS > 0


def test_metric_colors_distinct():
    # US5: Median (GOOD/green) and Mean (RED) must be visually distinct
    # colors, and RED must be the only new constant (F-08 — no GREEN/
    # NEUTRAL duplicates; TEXT is reused for the neutral Min).
    assert gha.GOOD != gha.RED
    assert gha.RED == "#ff5555"


def test_peel_facecolors_reveal_shape():
    import numpy as np
    q = np.array([0.1, 0.5, 0.9, 0.3])
    elem_layer = np.array([0, 1, 0, 1], dtype=np.int32)
    fc0 = gha._peel_facecolors(q, elem_layer, n_layers=2, k=0)
    fc1 = gha._peel_facecolors(q, elem_layer, n_layers=2, k=1)
    assert fc0.shape == (4, 4)
    # k=1 reveals layer 0 -> those rows must differ from the k=0 (unrevealed) coloring.
    assert not np.array_equal(fc0[elem_layer == 0], fc1[elem_layer == 0])
    # layer 1 elements are untouched at k=1 (still quality-colored).
    assert np.array_equal(fc0[elem_layer == 1], fc1[elem_layer == 1])


def test_draw_peel_hist_bin_totals_invariant():
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(0)
    q = rng.uniform(0.0, 1.0, size=200)
    elem_layer = rng.integers(0, 4, size=200).astype(np.int32)
    D = dict(ymax=100)

    ref_counts, _ = np.histogram(q, bins=gha.HBINS, range=(0.0, 1.0))

    for k in range(5):
        fig, (ax_mesh, ax_hist) = plt.subplots(1, 2)
        gha._draw_peel_hist(ax_hist, q, elem_layer, n_layers=4, k=k, D=D)
        # Sum the heights of all bar patches per bin position and confirm
        # the stacked total matches the plain histogram (partition property).
        totals = np.zeros(gha.HBINS)
        for patch in ax_hist.patches:
            bin_idx = int(round(patch.get_x() * gha.HBINS))
            bin_idx = min(max(bin_idx, 0), gha.HBINS - 1)
            totals[bin_idx] += patch.get_height()
        plt.close(fig)
        assert np.allclose(totals, ref_counts), f"k={k} stacked totals diverge from plain histogram"
