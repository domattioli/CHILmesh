#!/usr/bin/env python3
"""Render the size-controlled cartogram gallery (CHILmesh #219).

One figure, one row per mesh, five columns:

    plain plot_filled | noncontig | unrolled | noncontig, reference="mean"
    | unrolled, reference="mean"

The field is element skew quality (``mesh.elem_quality()``). Each row shares
one colour bar. Headless (Agg).

Usage:
    python scripts/gallery_cartograms.py [--out docs/gallery/cartograms.png]
                                         [--mesh PATH.fort.14]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402

from chilmesh import CHILmesh, examples  # noqa: E402
from chilmesh.chilplotting import cartogram, plot_filled  # noqa: E402

COLUMNS = (
    "plot_filled (plain)",
    "noncontig",
    "unrolled",
    'noncontig\nreference="mean"',
    'unrolled\nreference="mean"',
)
CMAP = "inferno"


def _meshes(extra):
    rows = [(n, getattr(examples, n)()) for n in ("annulus", "donut", "block_o")]
    if extra:
        rows.append((Path(extra).name, CHILmesh.read_from_fort14(extra)))
    return rows


def render(out: Path, extra=None) -> Path:
    rows = _meshes(extra)
    n = len(rows)
    fig = plt.figure(figsize=(21, 4.6 * n), dpi=150, layout="constrained")
    gs = fig.add_gridspec(n, 6, width_ratios=[1, 1, 1, 1, 1, 0.04])
    for r, (name, mesh) in enumerate(rows):
        q = mesh.elem_quality()[0]
        axes = [fig.add_subplot(gs[r, c]) for c in range(5)]
        plot_filled(mesh.points, mesh.connectivity_list, values=q, cmap=CMAP,
                    edgecolor="none", ax=axes[0])
        cartogram(mesh, q, kind="noncontig", colorbar=False, ax=axes[1])
        cartogram(mesh, q, kind="unrolled", colorbar=False, ax=axes[2])
        cartogram(mesh, q, kind="noncontig", reference="mean", colorbar=False,
                  ax=axes[3])
        cartogram(mesh, q, kind="unrolled", reference="mean", colorbar=False,
                  ax=axes[4])
        for ax in axes:
            ax.set_xlabel("")
            ax.set_box_aspect(1)  # equal-size square panels keep titles aligned
        for ax in (axes[2], axes[4]):
            ax.set_xticks([])  # band position, not a coordinate
        for ax in (axes[0], axes[1], axes[3]):
            ax.tick_params(labelsize=8)
        if r == 0:
            for ax, title in zip(axes, COLUMNS):
                ax.set_title(title, fontsize=13, fontweight="bold")
        axes[0].set_ylabel(f"{name}\nn_elems = {mesh.n_elems}", fontsize=12,
                           fontweight="bold")
        pc = next(c for c in axes[0].collections if isinstance(c, PolyCollection))
        cb = fig.colorbar(pc, cax=fig.add_subplot(gs[r, 5]))
        cb.set_label("element skew quality", fontsize=10)
    fig.suptitle("Size-controlled cartograms of element quality", fontsize=16)
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="docs/gallery/cartograms.png")
    ap.add_argument("--mesh", default=None, help="optional extra row (.fort.14)")
    args = ap.parse_args(argv)
    print("wrote", render(Path(args.out), args.mesh))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
