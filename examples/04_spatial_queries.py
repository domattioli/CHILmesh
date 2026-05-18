"""Spatial queries: point-in-element and nearest-neighbor lookups.

Run:

    python examples/04_spatial_queries.py

Demonstrates the Phase 5 spatial-indexing API: locate the element
containing an arbitrary point, list elements within a radius, and find
the k nearest vertices to a query point.
"""
from __future__ import annotations

import numpy as np

import chilmesh


def main() -> None:
    mesh = chilmesh.examples.annulus()
    print(f"annulus: {mesh.n_verts} verts, {mesh.n_elems} elems")

    query = np.array([0.5, 0.0])
    print(f"\nquery point: {query.tolist()}")

    elem_id = mesh.find_element(query)
    if elem_id >= 0:
        verts = mesh.connectivity_list[elem_id]
        verts = verts[verts >= 0]
        print(f"  containing element: {elem_id} (verts {verts.tolist()})")
    else:
        print(f"  containing element: <outside mesh>")

    in_radius = mesh.find_elements_in_radius(query, radius=0.2)
    print(f"  elements within r=0.2: {len(in_radius)}")

    nearest = mesh.nearest_vertices(query, k=5)
    coords = mesh.points[nearest, :2]
    dists = np.linalg.norm(coords - query, axis=1)
    print(f"  5 nearest vertices:")
    for vid, (x, y), d in zip(nearest, coords, dists):
        print(f"    vert {int(vid):4d}  at ({x:+.3f}, {y:+.3f})  d={d:.3f}")


if __name__ == "__main__":
    main()
