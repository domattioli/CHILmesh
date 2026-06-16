#!/usr/bin/env python3
"""Regenerate the 4-panel mesh-concept figure.

Distinguishes four related-but-distinct constructs on a shared L-shape domain:

  1. Distance field  d(x)  -- scalar: distance from each interior point to the
     nearest boundary (a field, not a curve).
  2. Medial axis           -- the ridge of d(x): the continuous locus
     equidistant from >=2 boundary points (pure geometry).
  3. Skeleton              -- a 1-px connected curve from topological thinning
     that preserves connectivity (a distinct algorithm, hence a distinct
     result from the medial axis).
  4. Layers (layerize)     -- concentric element bands peeled boundary-inward,
     which is what CHILmesh `_skeletonize` actually computes.

The figure makes the point that these are related, not identical: distance is a
field -> its ridge is the medial axis -> the skeleton is a thinned discrete
curve -> layers are concentric bands.

Requires scikit-image + scipy + matplotlib.
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.ndimage import distance_transform_edt
from skimage.morphology import medial_axis, skeletonize


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="docs/gallery/mesh_concepts.png",
        help="output PNG path (default: docs/gallery/mesh_concepts.png)",
    )
    parser.add_argument(
        "--dpi", type=int, default=130, help="output DPI (default: 130)"
    )
    args = parser.parse_args()

    OUTPUT = args.output
    DPI = args.dpi
    Path(OUTPUT).parent.mkdir(parents=True, exist_ok=True)

    N = 520
    m = np.zeros((N, N), bool)
    m[60:460, 60:300] = True
    m[300:460, 60:460] = True
    d = distance_transform_edt(m)
    ma = medial_axis(m)
    sk = skeletonize(m)
    step = d[m].max() / 8.0
    lay = np.where(m, np.ceil(d / step), 0).astype(int)

    fig, ax = plt.subplots(1, 4, figsize=(18, 5.2))
    faint = ListedColormap(['white', '#e9edf1'])
    im = ax[0].imshow(np.where(m, d, np.nan), cmap='viridis')
    ax[0].set_title("Distance field  d(x)", fontsize=13, fontweight='bold')
    ax[0].text(0.5, -0.08, "scalar: distance from each interior point\nto the nearest boundary", ha='center', va='top', transform=ax[0].transAxes, fontsize=9)
    fig.colorbar(im, ax=ax[0], fraction=0.046, pad=0.04)
    ax[1].imshow(m, cmap=faint)
    yy, xx = np.where(ma); ax[1].scatter(xx, yy, s=0.6, c='#ba0c2f')
    ax[1].set_title("Medial axis", fontsize=13, fontweight='bold')
    ax[1].text(0.5, -0.08, "ridge of d(x): continuous locus equidistant\nfrom >=2 boundary points (pure geometry)", ha='center', va='top', transform=ax[1].transAxes, fontsize=9)
    ax[2].imshow(m, cmap=faint)
    yy, xx = np.where(sk); ax[2].scatter(xx, yy, s=0.6, c='#1f6feb')
    ax[2].set_title("Skeleton", fontsize=13, fontweight='bold')
    ax[2].text(0.5, -0.08, "topological thinning: 1-px connected curve that\npreserves connectivity (distinct algorithm -> distinct result)", ha='center', va='top', transform=ax[2].transAxes, fontsize=9)
    lc = ax[3].imshow(np.where(m, lay, np.nan), cmap='Spectral')
    ax[3].set_title("Layers  (layerize)", fontsize=13, fontweight='bold')
    ax[3].text(0.5, -0.08, "concentric element bands peeled boundary-inward\n- what CHILmesh `_skeletonize` computes", ha='center', va='top', transform=ax[3].transAxes, fontsize=9)
    fig.colorbar(lc, ax=ax[3], fraction=0.046, pad=0.04, label='layer index')
    for a in ax: a.set_xticks([]); a.set_yticks([]); a.set_aspect('equal')
    fig.suptitle("Related, not identical:  distance is a field  ->  its ridge is the medial axis  ->  the skeleton is a thinned discrete curve  ->  layers are concentric bands", fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(OUTPUT, dpi=DPI, bbox_inches='tight')


if __name__ == "__main__":
    main()
