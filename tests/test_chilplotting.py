"""Tests for the standalone chilplotting module (CHILmesh #75).

Covers the generator-agnostic array primitives — they must work on bare
``points`` + ``connectivity`` with NO CHILmesh object — and confirms the mixin
still delegates correctly so the historical API is unchanged.
"""
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

import chilmesh
from chilmesh import chilplotting as cp
from chilmesh import examples


# A bare two-triangle square — the kind of raw output any mesh generator emits.
SQUARE_PTS = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
SQUARE_TRIS = np.array([[0, 1, 2], [0, 2, 3]])
SQUARE_QUAD = np.array([[0, 1, 2, 3]])


def _artists(ax):
    return len(ax.collections) + len(ax.lines) + len(ax.patches)


def test_module_importable_from_package():
    assert chilmesh.chilplotting is cp
    for name in ("plot_mesh", "plot_filled", "plot_edges", "plot_points",
                 "plot_quality_histogram", "build_polygons", "unique_edges"):
        assert hasattr(cp, name)


# --- array primitives: no CHILmesh needed ---------------------------------

def test_unique_edges_triangles():
    edges = cp.unique_edges(SQUARE_TRIS)
    # 5 unique edges: 4 perimeter + 1 shared diagonal
    assert edges.shape == (5, 2)
    assert (edges[:, 0] < edges[:, 1]).all()


def test_unique_edges_handles_padded_triangle():
    # 4-col row with a repeated vertex is a triangle -> 3 edges, not 4
    padded = np.array([[0, 1, 2, 0]])
    assert cp.unique_edges(padded).shape == (3, 2)


def test_build_polygons_tri_quad_mixed():
    assert len(cp.build_polygons(SQUARE_PTS, SQUARE_TRIS)) == 2
    assert cp.build_polygons(SQUARE_PTS, SQUARE_TRIS)[0].shape == (3, 2)
    assert cp.build_polygons(SQUARE_PTS, SQUARE_QUAD)[0].shape == (4, 2)
    mixed = np.array([[0, 1, 2, 0], [0, 1, 2, 3]])  # padded tri + quad
    polys = cp.build_polygons(SQUARE_PTS, mixed)
    assert sorted(len(p) for p in polys) == [3, 4]


def test_plot_mesh_on_raw_arrays():
    fig, ax = cp.plot_mesh(SQUARE_PTS, SQUARE_TRIS)
    assert _artists(ax) >= 1
    assert ax.get_aspect() == 1.0
    plt.close(fig)


def test_plot_filled_scalar_field_with_colorbar():
    values = np.array([0.3, 0.9])  # one per triangle
    fig, ax = cp.plot_filled(SQUARE_PTS, SQUARE_TRIS, values=values,
                             cmap="viridis", colorbar=True, clabel="quality")
    assert len(ax.collections) == 1
    plt.close(fig)


def test_plot_filled_value_length_mismatch_raises():
    with pytest.raises(ValueError):
        cp.plot_filled(SQUARE_PTS, SQUARE_TRIS, values=np.array([0.1]))


def test_plot_points_from_coords_and_ids():
    fig, ax = cp.plot_points(SQUARE_PTS)
    assert _artists(ax) >= 1
    plt.close(fig)
    fig, ax = cp.plot_points(SQUARE_PTS, ids=[0, 2])
    plt.close(fig)
    fig, ax = cp.plot_points(SQUARE_PTS, coords=np.array([[0.5, 0.5]]))
    plt.close(fig)


def test_plot_quality_histogram_array():
    q = np.linspace(0, 1, 50)
    fig, ax = cp.plot_quality_histogram(q, bins=10)
    assert _artists(ax) >= 1
    assert ax.get_xlim() == (0.0, 1.0)
    plt.close(fig)


def test_configure_axes_equal_aspect():
    fig, ax = plt.subplots()
    cp.configure_axes(SQUARE_PTS, ax=ax)
    assert ax.get_aspect() == 1.0
    xmin, xmax = ax.get_xlim()
    assert xmin <= 0.0 and xmax >= 1.0
    plt.close(fig)


# --- mesh functions match the mixin (delegation sanity) -------------------

@pytest.fixture(params=["annulus", "donut", "quad_2x2"])
def mesh(request):
    return getattr(examples, request.param)()


def test_module_and_mixin_agree_on_edges(mesh):
    fig1, ax1 = mesh.plot_boundary()
    fig2, ax2 = cp.plot_boundary(mesh)
    assert len(ax1.collections) == len(ax2.collections) >= 1
    plt.close(fig1)
    plt.close(fig2)


def test_mesh_plot_filled_via_module(mesh):
    fig, ax = cp.plot_filled(mesh.points, mesh.connectivity_list,
                             facecolor="C1")
    assert len(ax.collections) == 1
    plt.close(fig)
