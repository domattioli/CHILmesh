"""Test layers discretization on the annulus domain and generate visualization."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from chilmesh import CHILmesh


def test_structured_grid_elem2elem():
    """Element 0 should have <= 3 neighbors (triangle), not 10+."""
    pts = np.array([[float(i), float(j)] for j in range(5) for i in range(5)])
    conn = []
    for j in range(4):
        for i in range(4):
            v0 = j * 5 + i
            conn.append([v0, v0 + 1, v0 + 6])
            conn.append([v0, v0 + 6, v0 + 5])
    conn = np.array(conn)

    mesh = CHILmesh(connectivity=conn, points=pts, grid_name='test_grid')

    e2e = mesh.adjacencies['Edge2Elem']
    elem2elem = [set() for _ in range(mesh.n_elems)]
    for e1, e2 in e2e:
        if e1 >= 0 and e2 >= 0:
            elem2elem[e1].add(e2)
            elem2elem[e2].add(e1)

    assert len(elem2elem[0]) <= 3, f"Element 0 has {len(elem2elem[0])} neighbors, expected <= 3"


def test_structured_grid_boundary_edges():
    """4x4 quad grid split into triangles should have 16 boundary edges."""
    pts = np.array([[float(i), float(j)] for j in range(5) for i in range(5)])
    conn = []
    for j in range(4):
        for i in range(4):
            v0 = j * 5 + i
            conn.append([v0, v0 + 1, v0 + 6])
            conn.append([v0, v0 + 6, v0 + 5])
    conn = np.array(conn)

    mesh = CHILmesh(connectivity=conn, points=pts, grid_name='test_grid')
    bnd = mesh.boundary_edges()
    assert len(bnd) == 16, f"Expected 16 boundary edges, got {len(bnd)}"


def test_structured_grid_layers():
    """4x4 grid should produce >= 2 layers with all elements assigned."""
    pts = np.array([[float(i), float(j)] for j in range(5) for i in range(5)])
    conn = []
    for j in range(4):
        for i in range(4):
            v0 = j * 5 + i
            conn.append([v0, v0 + 1, v0 + 6])
            conn.append([v0, v0 + 6, v0 + 5])
    conn = np.array(conn)

    mesh = CHILmesh(connectivity=conn, points=pts, grid_name='test_grid')
    assert mesh.n_layers >= 2

    all_assigned = set()
    for layer_oe in mesh.layers['OE']:
        all_assigned.update(layer_oe)
    for layer_ie in mesh.layers['IE']:
        all_assigned.update(layer_ie)
    assert len(all_assigned) == mesh.n_elems, f"Only {len(all_assigned)}/{mesh.n_elems} elements assigned"


def test_no_duplicate_layer_assignments():
    """Each element should appear in exactly one layer."""
    pts = np.array([[float(i), float(j)] for j in range(5) for i in range(5)])
    conn = []
    for j in range(4):
        for i in range(4):
            v0 = j * 5 + i
            conn.append([v0, v0 + 1, v0 + 6])
            conn.append([v0, v0 + 6, v0 + 5])
    conn = np.array(conn)

    mesh = CHILmesh(connectivity=conn, points=pts, grid_name='test_grid')

    all_elems = []
    for layer_oe in mesh.layers['OE']:
        all_elems.extend(layer_oe)
    for layer_ie in mesh.layers['IE']:
        all_elems.extend(layer_ie)
    assert len(all_elems) == len(set(all_elems)), "Duplicate element assignments across layers"


def test_annulus_layers():
    """Annulus domain should produce multiple layers with all elements assigned."""
    import chilmesh
    mesh = chilmesh.examples.annulus()

    assert mesh.n_layers >= 2, f"Annulus should have >= 2 layers, got {mesh.n_layers}"

    all_assigned = set()
    for layer_oe in mesh.layers['OE']:
        all_assigned.update(layer_oe)
    for layer_ie in mesh.layers['IE']:
        all_assigned.update(layer_ie)
    assert len(all_assigned) == mesh.n_elems, f"Only {len(all_assigned)}/{mesh.n_elems} elements assigned"


def test_edge2elem_sentinel():
    """Boundary edges should have -1 sentinel, not 0."""
    pts = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], dtype=float)
    conn = np.array([[0, 1, 2], [0, 2, 3]])
    mesh = CHILmesh(connectivity=conn, points=pts, grid_name='sentinel_test')

    e2e = mesh.adjacencies['Edge2Elem']
    bnd_edges = mesh.boundary_edges()

    for edge_id in bnd_edges:
        assert e2e[edge_id, 1] == -1, f"Boundary edge {edge_id} has e2={e2e[edge_id, 1]}, expected -1"


def generate_annulus_visualization():
    """Generate layer visualization for the annulus domain."""
    import chilmesh
    mesh = chilmesh.examples.annulus()

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: mesh with layers
    _, ax = mesh.plot_layer(ax=axes[0])
    ax.set_title(f"Annulus Mesh Layers ({mesh.n_layers} layers, {mesh.n_elems} elements)")

    # Right: mesh quality
    q, _, stats = mesh.elem_quality()
    _, ax = mesh.plot_quality(ax=axes[1])
    ax.set_title(f"Element Quality (Median: {np.median(q):.3f}, Std: {np.std(q):.3f})")

    plt.tight_layout()
    output_dir = os.path.join(os.path.dirname(__file__), 'output')
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, 'annulus_layers.png')
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved visualization to {output_path}")
    return output_path


if __name__ == '__main__':
    test_structured_grid_elem2elem()
    print("PASS: test_structured_grid_elem2elem")
    test_structured_grid_boundary_edges()
    print("PASS: test_structured_grid_boundary_edges")
    test_structured_grid_layers()
    print("PASS: test_structured_grid_layers")
    test_no_duplicate_layer_assignments()
    print("PASS: test_no_duplicate_layer_assignments")
    test_annulus_layers()
    print("PASS: test_annulus_layers")
    test_edge2elem_sentinel()
    print("PASS: test_edge2elem_sentinel")
    print("\n=== All tests passed ===\n")
    generate_annulus_visualization()
