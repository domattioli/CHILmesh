#!/usr/bin/env python3
"""Regenerate the 4-panel mesh-concept figure.

Distinguishes four related-but-distinct constructs on a shared L-shape domain:

  1. Distance field  d(x)  -- a scalar field: distance from each interior point
     to the nearest boundary (a field, not a curve).
  2. Medial axis           -- the ridge of d(x): the continuous locus
     equidistant from >=2 boundary points (pure geometry of the domain).
  3. Skeleton              -- a 1-wide connectivity-preserving curve from
     topological thinning (a distinct algorithm, hence a distinct result).
  4. Layers (layerize)     -- actual mesh ELEMENTS grouped into concentric
     bands peeled boundary-inward, which is what CHILmesh `_layerize`
     (formerly `_skeletonize`) computes. See issue #221.

The point: these are related, not identical. distance is a field -> its ridge
is the medial axis -> the skeleton is a thinned discrete curve -> layers group
mesh elements into concentric bands.

Requires scikit-image + scipy + matplotlib + chilmesh.
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.collections import PolyCollection
from scipy.ndimage import distance_transform_edt
from scipy.spatial import Delaunay
from skimage.morphology import medial_axis, skeletonize

import chilmesh

N = 520


def _in_L(p):
    x, y = p[:, 0], p[:, 1]
    return (((x >= 60) & (x < 300) & (y >= 60) & (y < 460)) |
            ((x >= 60) & (x < 460) & (y >= 300) & (y < 460)))


def build_figure(output: Path, dpi: int) -> None:
    # shared L-shape domain (raster) for the continuous panels
    m = np.zeros((N, N), bool)
    m[60:460, 60:300] = True
    m[300:460, 60:460] = True
    d = distance_transform_edt(m)
    ma = medial_axis(m)
    sk = skeletonize(m)

    # triangulated L-shape mesh for the element-level layers panel
    gx, gy = np.meshgrid(np.linspace(62, 458, 46), np.linspace(62, 458, 46))
    P = np.c_[gx.ravel(), gy.ravel()]
    P = P[_in_L(P)]
    conn = Delaunay(P).simplices
    conn = conn[_in_L(P[conn].mean(1))]
    mesh = chilmesh.CHILmesh(connectivity=conn.astype(int), points=P, compute_layers=True)
    elay = np.zeros(mesh.n_elems, int)
    for L in range(mesh.n_layers):
        for arr in (mesh.layers["OE"][L], mesh.layers["IE"][L]):
            a = np.asarray(arr, int)
            if a.size:
                elay[a] = L

    fig, ax = plt.subplots(1, 4, figsize=(18, 5.2))
    faint = ListedColormap(["white", "#e9edf1"])

    im = ax[0].imshow(np.where(m, d, np.nan), cmap="viridis")
    ax[0].set_title("Distance field  d(x)", fontsize=13, fontweight="bold")
    ax[0].text(0.5, -0.08, "scalar field: distance from each point to the boundary",
               ha="center", va="top", transform=ax[0].transAxes, fontsize=9)
    fig.colorbar(im, ax=ax[0], fraction=0.046, pad=0.04)

    ax[1].imshow(m, cmap=faint)
    yy, xx = np.where(ma)
    ax[1].scatter(xx, yy, s=0.6, c="#ba0c2f")
    ax[1].set_title("Medial axis", fontsize=13, fontweight="bold")
    ax[1].text(0.5, -0.08, "ridge of d(x): locus equidistant from >=2 boundary points",
               ha="center", va="top", transform=ax[1].transAxes, fontsize=9)

    ax[2].imshow(m, cmap=faint)
    yy, xx = np.where(sk)
    ax[2].scatter(xx, yy, s=0.6, c="#1f6feb")
    ax[2].set_title("Skeleton", fontsize=13, fontweight="bold")
    ax[2].text(0.5, -0.08, "topological thinning: 1-wide connectivity-preserving curve",
               ha="center", va="top", transform=ax[2].transAxes, fontsize=9)

    polys = [P[c] for c in conn]
    pc = PolyCollection(polys, array=elay, cmap="Spectral", edgecolors="white", linewidths=0.25)
    ax[3].add_collection(pc)
    ax[3].set_xlim(40, 480)
    ax[3].set_ylim(40, 480)
    ax[3].invert_yaxis()
    ax[3].set_title("Layers  (layerize)", fontsize=13, fontweight="bold")
    ax[3].text(0.5, -0.08, "mesh ELEMENTS grouped into concentric bands peeled boundary-inward",
               ha="center", va="top", transform=ax[3].transAxes, fontsize=9)
    fig.colorbar(pc, ax=ax[3], fraction=0.046, pad=0.04, label="layer index")

    for a in ax:
        a.set_xticks([])
        a.set_yticks([])
        a.set_aspect("equal")
    fig.suptitle("Related, not identical:  distance is a field  ->  its ridge is the medial axis  ->  "
                 "the skeleton is a thinned curve  ->  layers group mesh elements into concentric bands",
                 fontsize=12, y=1.02)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    print(f"wrote {output}  ({mesh.n_elems} elements, {mesh.n_layers} layers)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--output", default="docs/gallery/mesh_concepts.png", help="output PNG path")
    ap.add_argument("--dpi", type=int, default=130, help="figure DPI")
    args = ap.parse_args()
    build_figure(Path(args.output), args.dpi)


if __name__ == "__main__":
    main()
