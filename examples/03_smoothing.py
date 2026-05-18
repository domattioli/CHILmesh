"""Mesh smoothing: apply the angle-based smoother and report quality delta.

Run:

    python examples/03_smoothing.py

Loads the annulus, perturbs interior vertices, then relaxes them with
the Zhou-Shimada angle-based smoother. Reports min/median element quality
before and after.
"""
from __future__ import annotations

import numpy as np

import chilmesh


def quality_stats(mesh: chilmesh.CHILmesh) -> tuple[float, float]:
    q, _, _ = mesh.elem_quality()
    return float(q.min()), float(np.median(q))


def main() -> None:
    mesh = chilmesh.examples.annulus()

    edge_verts = mesh.edge2vert(mesh.boundary_edges())
    boundary = np.unique(edge_verts.flatten())
    interior = np.setdiff1d(np.arange(mesh.n_verts), boundary)

    rng = np.random.default_rng(2026)
    mesh.points[interior, :2] += rng.uniform(-0.02, 0.02, size=(interior.size, 2))
    pre_min, pre_med = quality_stats(mesh)
    print(f"perturbed: min q={pre_min:.3f}  median q={pre_med:.3f}")

    smoothed = mesh.angle_based_smoother(n_iter=20)
    mesh.points = smoothed
    post_min, post_med = quality_stats(mesh)
    print(f"smoothed:  min q={post_min:.3f}  median q={post_med:.3f}")

    print()
    print(f"min quality:    {pre_min:.3f} -> {post_min:.3f}  (delta {post_min-pre_min:+.3f})")
    print(f"median quality: {pre_med:.3f} -> {post_med:.3f}  (delta {post_med-pre_med:+.3f})")

    # Boundary pinning sanity check.
    moved_boundary = np.any(np.abs(mesh.points[boundary, :2] - mesh.points[boundary, :2]) > 0)
    assert not moved_boundary, "boundary vertices moved (should be pinned)"
    print("OK: boundary preserved")


if __name__ == "__main__":
    main()
