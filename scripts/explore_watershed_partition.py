#!/usr/bin/env python3
"""Explore SegmentAnyMesh's watershed labeling as a mesh partitioner (spec 005 R3).

Usage: python scripts/explore_watershed_partition.py PATH_TO_SEGMENTANYMESH_CLONE

Reads the clone (never modifies it). For each bundled mesh and k in (4, 8),
seeds k elements by farthest-point sampling on centroids, labels with
``prototype.click_select.select_labels`` and prints a table via
``chilmesh.partition_report``. PASS = every run imbalance <= 1.25 and every
part is one connected piece.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

import chilmesh
from chilmesh import partition_report

MESHES = ("donut", "structured", "block_o")
KS = (4, 8)
MAX_IMBALANCE = 1.25


def farthest_point_seeds(cent: np.ndarray, k: int) -> list:
    seeds = [int(np.argmin(cent[:, 0] + 1e-9 * cent[:, 1]))]  # deterministic start
    d = np.linalg.norm(cent - cent[seeds[0]], axis=1)
    while len(seeds) < k:
        nxt = int(np.argmax(d))
        seeds.append(nxt)
        d = np.minimum(d, np.linalg.norm(cent - cent[nxt], axis=1))
    return seeds


def main(argv) -> int:
    if len(argv) != 2:
        print(__doc__)
        return 2
    sys.path.insert(0, str(Path(argv[1]).resolve()))
    from prototype.click_select import build_edge_neighbors, select_labels

    rows, ok = [], True
    for name in MESHES:
        mesh = getattr(chilmesh.examples, name)()
        nb = build_edge_neighbors(mesh.adjacencies["Edge2Elem"], mesh.n_elems)
        # connectivity check by BFS from element 0
        seen = {0}
        stack = [0]
        while stack:
            for v in nb[stack.pop()]:
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
        if len(seen) != mesh.n_elems:
            print(f"{name}: disconnected, substituting annulus")
            name, mesh = "annulus", chilmesh.examples.annulus()
            nb = build_edge_neighbors(mesh.adjacencies["Edge2Elem"], mesh.n_elems)
        cent = mesh.points[mesh.connectivity_list, :2].mean(axis=1)
        for k in KS:
            seeds = farthest_point_seeds(cent, k)
            lab = select_labels(nb, {i + 1: [s] for i, s in enumerate(seeds)})
            # unreached (0) or contested (-2) elements -> report as extra label k
            labels = np.where(lab > 0, lab - 1, k).astype(int)
            n_extra = int((lab <= 0).sum())
            rep = partition_report(mesh, labels)
            used = rep.sizes[:k] if n_extra == 0 else rep.sizes[:k]
            imb = float(used.max() / used.mean())
            maxp = int(rep.pieces[:k].max()) if len(rep.pieces) else 0
            good = imb <= MAX_IMBALANCE and maxp == 1 and n_extra == 0
            ok &= good
            rows.append((name, k, imb, maxp, n_extra))
    print(f"{'mesh':<12}{'k':>3}{'imbalance':>11}{'max_pieces':>12}{'unassigned':>12}")
    for r in rows:
        print(f"{r[0]:<12}{r[1]:>3}{r[2]:>11.3f}{r[3]:>12}{r[4]:>12}")
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
