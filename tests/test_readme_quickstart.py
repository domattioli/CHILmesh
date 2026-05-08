"""End-to-end regeneration of the README quickstart figure.

Renders the 3-row pipeline that lives in ``README.md``:

    Row 1: Raw annulus from chilmesh.examples.annulus()
    Row 2: ADMESH warm-start truss applied to row 1
    Row 3: FEM smoother applied to row 2

All three rows use chilmesh's own optimizer — no external admesh package
required. Each row asserts geometric validity before plotting.

If the test passes, the freshly-generated PNG at ``output/annulus_quickstart.png``
is the canonical README image.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest
import warnings

import chilmesh
from chilmesh import optimize_with_admesh_truss

REPO_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = REPO_ROOT / "output"
OUTPUT_PNG = OUTPUT_DIR / "annulus_quickstart.png"

_ANNULUS_SDF = lambda p: np.maximum(
    np.linalg.norm(p, axis=1) - 1.0,
    0.3 - np.linalg.norm(p, axis=1),
)


def _assert_valid(mesh, label: str) -> None:
    areas = mesh.signed_area()
    assert (areas > 0).all(), f"{label}: {(areas <= 0).sum()} non-positive-area elements"
    conn = mesh.connectivity_list
    assert (conn >= 0).all() and (conn < mesh.n_verts).all(), f"{label}: index OOB"
    if mesh.type == "Triangular":
        for row in conn[:, :3]:
            assert len(set(row.tolist())) == 3, f"{label}: degenerate row {row}"
    assert not np.isnan(mesh.interior_angles()).any(), f"{label}: NaN angles"
    seen = set()
    for oe in mesh.layers["OE"]:
        seen.update(int(e) for e in oe)
    for ie in mesh.layers["IE"]:
        seen.update(int(e) for e in ie)
    assert seen == set(range(mesh.n_elems)), f"{label}: layer cover incomplete"


def _plot_row(axs_row, mesh, label: str) -> None:
    """Render three CHILmesh views for one figure row."""
    _, ax = mesh.plot(ax=axs_row[0])
    mesh.plot_point(1, ax=ax)
    mesh.plot_edge(1, ax=ax)
    mesh.plot_elem(1, ax=ax)
    ax.set_title(f"{label}\nmesh + highlights")

    _, ax = mesh.plot_layer(ax=axs_row[1])
    ax.set_title(f"{label}\n{mesh.n_layers} layers")

    q, _, _ = mesh.elem_quality()
    _, ax = mesh.plot_quality(ax=axs_row[2])
    ax.set_title(f"{label}\nquality med={np.median(q):.2f} σ={np.std(q):.2f}")


def test_readme_annulus_regenerates_and_validates(tmp_path):
    # Row 1 — raw annulus
    mesh = chilmesh.examples.annulus()
    _assert_valid(mesh, "annulus (raw)")

    # Row 2 — warm-start truss (chilmesh-internal, no external admesh needed)
    # deltat=0.02 + track_best_quality=False: small steps to convergence, no
    # premature exit from the degeneracy guard (outer enforce_non_degradation
    # handles quality check after the fact).
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        truss = optimize_with_admesh_truss(
            mesh, _ANNULUS_SDF, enforce_non_degradation=False,
            deltat=0.02, track_best_quality=False,
        )
    _assert_valid(truss, "annulus (truss)")

    # Row 3 — FEM smoother on top of truss result
    smoothed = truss.copy()
    smoothed.smooth_mesh(method="fem", acknowledge_change=True)
    _assert_valid(smoothed, "annulus (FEM)")

    # Roundtrip fort.14 check
    out14 = tmp_path / "annulus_smoothed.fort.14"
    assert smoothed.write_to_fort14(str(out14), grid_name="annulus_smoothed")
    reloaded = chilmesh.CHILmesh.read_from_fort14(out14)
    _assert_valid(reloaded, "annulus (reloaded)")

    # 3-row figure
    fig, axs = plt.subplots(3, 3, figsize=(15, 13.5))
    fig.suptitle(
        "CHILmesh annulus: Raw Delaunay → Warm-Start Truss → FEM Smoother",
        fontsize=14,
    )

    _plot_row(axs[0], mesh,     "Row 1 — Raw Delaunay")
    _plot_row(axs[1], truss,    "Row 2 — + Truss (warm-start)")
    _plot_row(axs[2], smoothed, "Row 3 — + FEM Smoother")

    plt.tight_layout()
    plt.subplots_adjust(top=0.95)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PNG, dpi=120, bbox_inches="tight")
    plt.close(fig)

    assert OUTPUT_PNG.exists() and OUTPUT_PNG.stat().st_size > 5_000, (
        f"figure not written or too small: {OUTPUT_PNG}"
    )
