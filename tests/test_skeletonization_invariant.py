"""
Regression tests for skeletonization layer separation invariant (issue #74).

The layer separation invariant states that for any two layers k and m with |k-m| >= 2,
the set of vertices in layer k must be disjoint from the set of vertices in layer m.
This ensures that the skeletonization produces medially-correct concentric rings.
"""

import pytest
import numpy as np
from chilmesh import examples


@pytest.mark.parametrize(
    "fixture_name",
    ["annulus", "donut", "structured", "block_o"]
)
def test_layer_separation_invariant(fixture_name):
    """
    Verify that layers separated by 2+ do not share any vertices.

    For each fixture, enumerate all layer pairs (k, m) with |k-m| >= 2
    and assert that the set of vertices in layer k is disjoint from
    the set of vertices in layer m.
    """
    factory = getattr(examples, fixture_name)
    mesh = factory()

    # Collect vertices for each layer
    layer_vertices = []
    for layer_idx in range(mesh.n_layers):
        layer = mesh.get_layer(layer_idx)

        # Get all vertices from outer and inner elements
        vertices = set()
        for elem_list in [layer["OE"], layer["IE"]]:
            for elem_id in elem_list:
                # Get vertices for this element
                elem_verts = mesh.connectivity_list[elem_id]
                # Filter out padding (-1 values)
                for v in elem_verts:
                    if v >= 0:
                        vertices.add(int(v))

        layer_vertices.append(vertices)

    # Check layer separation invariant
    for k in range(mesh.n_layers):
        for m in range(mesh.n_layers):
            if abs(k - m) >= 2:
                # Layers k and m should have disjoint vertex sets
                intersection = layer_vertices[k] & layer_vertices[m]
                assert len(intersection) == 0, (
                    f"Layer separation invariant violated: "
                    f"layer {k} and layer {m} share vertices {intersection}"
                )
