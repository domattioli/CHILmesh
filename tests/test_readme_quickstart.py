"""End-to-end regeneration of the README quickstart figure.

Renders the 4x3 grid that lives in ``README.md``:

    Row 1: Original CHILmesh annulus
    Row 2: CHILmesh FEM smoother
    Row 3: ADMESH triangulate result (skipped if admesh not installed)
    Row 4: ADMESH + right-angle smoother (skipped if admesh not installed)

Each row asserts the mesh is geometrically valid before plotting.

If the test passes, the freshly-generated PNG at
``tests/output/annulus_quickstart.png`` is the canonical README image.
Re-run ``pytest tests/test_readme_quickstart.py`` after any change that
might affect the visualization and commit the updated PNG.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

import chilmesh

try:
    import admesh
    HAS_ADMESH = True
except ImportError:
    HAS_ADMESH = False

REPO_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = REPO_ROOT / "tests" / "output"
OUTPUT_PNG = OUTPUT_DIR / "annulus_quickstart.png"


def _assert_geometrically_valid(mesh, label: str) -> None:
    areas = mesh.signed_area()
    assert (areas > 0).all(), f"{label}: {(areas <= 0).sum()} elements with non-positive area"

    conn = mesh.connectivity_list
    assert (conn >= 0).all() and (conn < mesh.n_verts).all(), f"{label}: index out of range"
    if mesh.type == "Triangular":
        for row in conn[:, :3]:
            assert len(set(row.tolist())) == 3, f"{label}: degenerate triangle row {row}"

    angles = mesh.interior_angles()
    assert not np.isnan(angles).any(), f"{label}: NaN in interior angles"

    seen = set()
    for oe in mesh.layers["OE"]:
        seen.update(int(e) for e in oe)
    for ie in mesh.layers["IE"]:
        seen.update(int(e) for e in ie)
    assert seen == set(range(mesh.n_elems)), (
        f"{label}: layer cover misses {len(set(range(mesh.n_elems)) - seen)} elements"
    )


def _annulus_sdf(p: np.ndarray) -> np.ndarray:
    r = np.sqrt(p[:, 0] ** 2 + p[:, 1] ** 2)
    return np.maximum(r - 2.0, 1.0 - r)


def _plot_row(axs_row, mesh, label: str) -> None:
    """Render the three CHILmesh views for a single row of the figure."""
    _, ax = mesh.plot(ax=axs_row[0])
    mesh.plot_point(1, ax=ax)
    mesh.plot_edge(1, ax=ax)
    mesh.plot_elem(1, ax=ax)
    ax.set_title(f"{label}: mesh + highlighted entities")

    _, ax = mesh.plot_layer(ax=axs_row[1])
    ax.set_title(f"{label}: {mesh.n_layers} mesh layers")

    q, _, _ = mesh.elem_quality()
    _, ax = mesh.plot_quality(ax=axs_row[2])
    ax.set_title(f"{label}: quality (median {np.median(q):.2f}, std {np.std(q):.2f})")


def test_readme_annulus_regenerates_and_validates(tmp_path):
    # 1. Original + FEM-smoothed (CHILmesh-only).
    mesh = chilmesh.examples.annulus()
    _assert_geometrically_valid(mesh, "annulus (original)")

    smoothed = mesh.copy()
    smoothed.smooth_mesh(method="fem", acknowledge_change=True)
    _assert_geometrically_valid(smoothed, "annulus (FEM-smoothed)")

    # Roundtrip-write the smoothed mesh and confirm it reloads.
    out14 = tmp_path / "annulus_smoothed.fort.14"
    assert smoothed.write_to_fort14(str(out14), grid_name="annulus_smoothed") is True
    reloaded = chilmesh.CHILmesh.read_from_fort14(out14)
    _assert_geometrically_valid(reloaded, "annulus (reloaded)")

    # 2. ADMESH + ADMESH+right-smoother (skipped if admesh not installed).
    if HAS_ADMESH:
        pfix = np.array(
            [[2 * np.cos(a), 2 * np.sin(a)]
             for a in np.linspace(0, 2 * np.pi, 80, endpoint=False)] +
            [[1 * np.cos(a), 1 * np.sin(a)]
             for a in np.linspace(0, 2 * np.pi, 40, endpoint=False)]
        )
        domain = admesh.Domain(_annulus_sdf, bbox=(-2.5, -2.5, 2.5, 2.5), pfix=pfix)
        admesh_mesh = admesh.triangulate(domain, h_max=0.15)
        mesh_adm = chilmesh.from_admesh(admesh_mesh)
        _assert_geometrically_valid(mesh_adm, "annulus (admesh)")

        p_right, t_right = admesh.smooth_for_quadrangulation(
            admesh_mesh.nodes, admesh_mesh.elements,
            _annulus_sdf, n_outer=3, pair_hint=True,
        )
        mesh_right = chilmesh.from_admesh_arrays(p_right, t_right)
        _assert_geometrically_valid(mesh_right, "annulus (admesh+right)")
    else:
        mesh_adm = mesh_right = None

    # 3. Render the 4x3 README quickstart figure.
    n_rows = 4 if HAS_ADMESH else 2
    fig, axs = plt.subplots(n_rows, 3, figsize=(15, 4.5 * n_rows))
    fig.suptitle(
        "CHILmesh annulus: Original vs FEM-Smoothed"
        + (" vs ADMESH vs ADMESH + Right-Angle Smoother" if HAS_ADMESH else ""),
        fontsize=16,
    )

    _plot_row(axs[0], mesh, "Original")
    _plot_row(axs[1], smoothed, "FEM-Smoothed")
    if HAS_ADMESH:
        _plot_row(axs[2], mesh_adm, "ADMESH")
        _plot_row(axs[3], mesh_right, "ADMESH + Right-Angle Smoother")

    plt.tight_layout()
    plt.subplots_adjust(top=0.96)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PNG, dpi=120, bbox_inches="tight")
    plt.close(fig)

    assert OUTPUT_PNG.exists() and OUTPUT_PNG.stat().st_size > 5_000, (
        f"figure not written or implausibly small: {OUTPUT_PNG}"
    )
