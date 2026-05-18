"""Quickstart: load a built-in mesh, print stats, save a quick plot.

Run from repo root or after `pip install chilmesh`:

    python examples/01_quickstart.py

Writes `annulus_quickstart.png` to the current directory.
"""
from __future__ import annotations

import matplotlib.pyplot as plt

import chilmesh


def main() -> None:
    mesh = chilmesh.examples.annulus()

    print(f"chilmesh {chilmesh.__version__}")
    print(f"  fixture: annulus")
    print(f"  vertices: {mesh.n_verts}")
    print(f"  elements: {mesh.n_elems}")
    print(f"  edges:    {mesh.n_edges}")
    print(f"  layers:   {mesh.n_layers}")

    quality, _, _ = mesh.elem_quality()
    print(f"  quality:  min={quality.min():.3f}  median={float(__import__('numpy').median(quality)):.3f}")

    fig, ax = plt.subplots(figsize=(6, 6))
    mesh.plot(ax=ax)
    ax.set_title(f"annulus ({mesh.n_elems} elements)")
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig("annulus_quickstart.png", dpi=120)
    print("  wrote annulus_quickstart.png")


if __name__ == "__main__":
    main()
