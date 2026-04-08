"""End-to-end regeneration of the README quickstart figure.

This test loads the bundled annulus, asserts the mesh is geometrically
valid (positive signed area, well-formed connectivity, finite interior
angles, complete layer cover), runs the FEM smoother, asserts the result
is *still* valid, and re-renders the 2x3 visualization that lives in
``README.md``.

If the test passes, the freshly-generated PNG at
``tests/output/annulus_quickstart.png`` is the canonical README image.
Re-run ``pytest tests/test_readme_quickstart.py`` after any change that
might affect the visualization and commit the updated PNG.

The test does not commit anything itself; CI runs it as a smoke check.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import chilmesh

REPO_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = REPO_ROOT / "tests" / "output"
OUTPUT_PNG = OUTPUT_DIR / "annulus_quickstart.png"


def _assert_geometrically_valid(mesh, label: str) -> None:
    """The audit's working definition of "geometrically valid" 2D mesh.

    Aggregates the invariants the rest of the suite already enforces:

    - every element has strictly positive signed area (CCW orientation),
    - the connectivity table contains no out-of-range or duplicated
      vertices within an element,
    - interior angles are real and finite (no NaN slipping past
      ``elem_quality``'s degenerate-element zero-out),
    - the layer decomposition is a disjoint cover of every element.
    """
    areas = mesh.signed_area()
    assert (areas > 0).all(), f"{label}: {(areas <= 0).sum()} elements with non-positive area"

    conn = mesh.connectivity_list
    assert (conn >= 0).all() and (conn < mesh.n_verts).all(), f"{label}: index out of range"
    if mesh.type == "Triangular":
        for row in conn[:, :3]:
            assert len(set(row.tolist())) == 3, f"{label}: degenerate triangle row {row}"

    angles = mesh.interior_angles()
    assert not np.isnan(angles).any(), f"{label}: NaN in interior angles"
    if mesh.type == "Triangular":
        np.testing.assert_allclose(
            angles[:, :3].sum(axis=1), 180.0, atol=1e-6,
            err_msg=f"{label}: triangle angles do not sum to 180",
        )

    seen = set()
    for oe in mesh.layers["OE"]:
        seen.update(int(e) for e in oe)
    for ie in mesh.layers["IE"]:
        seen.update(int(e) for e in ie)
    assert seen == set(range(mesh.n_elems)), (
        f"{label}: layer cover misses {len(set(range(mesh.n_elems)) - seen)} elements"
    )


def test_readme_annulus_regenerates_and_validates(tmp_path):
    # 1. Load + validate the original.
    mesh = chilmesh.examples.annulus()
    _assert_geometrically_valid(mesh, "annulus (original)")

    # 2. Smooth + validate the smoothed copy.
    smoothed = mesh.copy()
    smoothed.smooth_mesh(method="fem", acknowledge_change=True)
    _assert_geometrically_valid(smoothed, "annulus (FEM-smoothed)")

    # 3. Roundtrip-write the smoothed mesh and confirm it reloads.
    out14 = tmp_path / "annulus_smoothed.fort.14"
    assert smoothed.write_to_fort14(str(out14), grid_name="annulus_smoothed") is True
    reloaded = chilmesh.CHILmesh.read_from_fort14(out14)
    _assert_geometrically_valid(reloaded, "annulus (reloaded)")
    assert reloaded.n_elems == smoothed.n_elems
    assert reloaded.n_verts == smoothed.n_verts

    # 4. Render the 2x3 README quickstart figure.
    fig, axs = plt.subplots(2, 3, figsize=(18, 10))
    axs = axs.flatten()
    fig.suptitle("CHILmesh quickstart: annulus (original vs FEM-smoothed)", fontsize=16)

    _, ax = mesh.plot(ax=axs[0])
    mesh.plot_point(1, ax=ax)
    mesh.plot_edge(1, ax=ax)
    mesh.plot_elem(1, ax=ax)
    ax.set_title("Original: mesh + highlighted entities")

    _, ax = mesh.plot_layer(ax=axs[1])
    ax.set_title(f"Original: {mesh.n_layers} mesh layers")

    q0, _, _ = mesh.elem_quality()
    _, ax = mesh.plot_quality(ax=axs[2])
    ax.set_title(
        f"Original: quality (median {np.median(q0):.2f}, std {np.std(q0):.2f})"
    )

    _, ax = smoothed.plot(ax=axs[3])
    smoothed.plot_point(1, ax=ax)
    smoothed.plot_edge(1, ax=ax)
    smoothed.plot_elem(1, ax=ax)
    ax.set_title("Smoothed: mesh + highlighted entities")

    _, ax = smoothed.plot_layer(ax=axs[4])
    ax.set_title(f"Smoothed: {smoothed.n_layers} mesh layers")

    q1, _, _ = smoothed.elem_quality()
    _, ax = smoothed.plot_quality(ax=axs[5])
    ax.set_title(
        f"Smoothed: quality (median {np.median(q1):.2f}, std {np.std(q1):.2f})"
    )

    plt.tight_layout()
    plt.subplots_adjust(top=0.92)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PNG, dpi=150, bbox_inches="tight")
    plt.close(fig)

    assert OUTPUT_PNG.exists() and OUTPUT_PNG.stat().st_size > 5_000, (
        f"figure not written or implausibly small: {OUTPUT_PNG}"
    )
