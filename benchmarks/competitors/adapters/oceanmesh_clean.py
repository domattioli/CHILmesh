"""Standalone adapter: clean/smooth a mesh with oceanmesh and print one JSON line.

Run ONLY in a subprocess inside the oceanmesh environment (see envs/). oceanmesh
is GPL-licensed; this is the only file that imports it.

Input/output are ``.npz`` files with ``points`` (n,2) and ``cells`` (m,3).
``--mode frozen`` (Track 1) smooths with the boundary pinned and topology fixed.
``--mode full`` (Track 2) runs oceanmesh's full clean, which may delete vertices.
The oceanmesh function names below are unverified until research.md (T005) maps
them; a missing function exits non-zero and is recorded as ``failed``.
"""
import argparse
import json
import sys
import time

import numpy as np


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", choices=["frozen", "full"], required=True)
    ap.add_argument("--iters", type=int, default=50)
    args = ap.parse_args()

    import oceanmesh

    data = np.load(args.inp)
    points, cells = data["points"], data["cells"]
    t0 = time.perf_counter()
    if args.mode == "frozen":
        # Pin every boundary node (edges used by exactly one cell), smooth interior only.
        edges = np.sort(np.concatenate([cells[:, [0, 1]], cells[:, [1, 2]], cells[:, [2, 0]]]), axis=1)
        uniq, counts = np.unique(edges, axis=0, return_counts=True)
        pfix = points[np.unique(uniq[counts == 1])]
        points_out, cells_out = oceanmesh.laplacian2(points, cells, max_iter=args.iters, pfix=pfix)
    else:
        points_out, cells_out = oceanmesh.mesh_clean(points, cells)
    wall = time.perf_counter() - t0
    np.savez(args.out, points=np.asarray(points_out), cells=np.asarray(cells_out))
    print(json.dumps({"tool": "oceanmesh", "version": getattr(oceanmesh, "__version__", "unknown"),
                      "op_s": wall, "mode": args.mode}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
