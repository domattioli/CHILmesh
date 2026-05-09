"""Tests for CHILmeshPlotMixin in plot_utils.py.

Runs headless via the Agg backend — no display required.
Parametrized over fast fixtures only (no block_o).
"""
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from chilmesh import examples


SMALL_FIXTURES = ["annulus", "donut", "structured", "quad_2x2"]


@pytest.fixture(params=SMALL_FIXTURES)
def small_mesh(request):
    return getattr(examples, request.param)()


# ---------------------------------------------------------------------------
# Basic plot methods
# ---------------------------------------------------------------------------

def test_plot_returns_fig_ax(small_mesh):
    fig, ax = small_mesh.plot()
    assert isinstance(fig, plt.Figure)
    assert isinstance(ax, plt.Axes)
    plt.close(fig)


def test_plot_with_elem_color(small_mesh):
    fig, ax = small_mesh.plot(elem_color='steelblue')
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_subset_elem_ids(small_mesh):
    ids = np.arange(min(5, small_mesh.n_elems))
    fig, ax = small_mesh.plot(elem_ids=ids)
    plt.close(fig)


def test_plot_reuses_provided_ax(small_mesh):
    fig, ax = plt.subplots()
    _, returned_ax = small_mesh.plot(ax=ax)
    assert returned_ax is ax
    plt.close(fig)


# ---------------------------------------------------------------------------
# plot_edge
# ---------------------------------------------------------------------------

def test_plot_edge_returns_fig_ax(small_mesh):
    fig, ax = small_mesh.plot_edge()
    assert isinstance(fig, plt.Figure)
    assert isinstance(ax, plt.Axes)
    plt.close(fig)


def test_plot_edge_subset(small_mesh):
    ids = np.arange(min(10, small_mesh.n_edges))
    fig, ax = small_mesh.plot_edge(edge_ids=ids)
    plt.close(fig)


# ---------------------------------------------------------------------------
# plot_boundary and plot_interior_edges (new methods)
# ---------------------------------------------------------------------------

def test_plot_boundary_returns_fig_ax(small_mesh):
    fig, ax = small_mesh.plot_boundary()
    assert isinstance(fig, plt.Figure)
    assert isinstance(ax, plt.Axes)
    plt.close(fig)


def test_plot_boundary_non_empty(small_mesh):
    boundary = small_mesh.boundary_edges()
    assert len(boundary) > 0, "expected at least one boundary edge"


def test_plot_interior_edges_returns_fig_ax(small_mesh):
    fig, ax = small_mesh.plot_interior_edges()
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_boundary_and_interior_partition_all_edges(small_mesh):
    """Boundary + interior must cover every edge exactly once."""
    boundary = small_mesh.boundary_edges()
    interior = np.setdiff1d(np.arange(small_mesh.n_edges), boundary)
    combined = np.sort(np.concatenate([boundary, interior]))
    np.testing.assert_array_equal(combined, np.arange(small_mesh.n_edges))


# ---------------------------------------------------------------------------
# plot_elem / plot_face
# ---------------------------------------------------------------------------

def test_plot_elem_returns_fig_ax(small_mesh):
    fig, ax = small_mesh.plot_elem()
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_face_is_alias(small_mesh):
    fig, ax = small_mesh.plot_face()
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


# ---------------------------------------------------------------------------
# plot_point
# ---------------------------------------------------------------------------

def test_plot_point_vertex(small_mesh):
    fig, ax = small_mesh.plot_point(point_type='vertex')
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_point_edge_midpoint(small_mesh):
    fig, ax = small_mesh.plot_point(point_type='edge')
    plt.close(fig)


def test_plot_point_centroid(small_mesh):
    fig, ax = small_mesh.plot_point(point_type='centroid')
    plt.close(fig)


def test_plot_point_invalid_type_raises(small_mesh):
    with pytest.raises(ValueError, match="Unknown point_type"):
        small_mesh.plot_point(point_type='nonsense')


# ---------------------------------------------------------------------------
# plot_label
# ---------------------------------------------------------------------------

def test_plot_label_all(small_mesh):
    fig, ax = small_mesh.plot_label(label='all')
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_label_vertex_only(small_mesh):
    fig, ax = small_mesh.plot_label(label='vertex')
    plt.close(fig)


def test_plot_label_edge_only(small_mesh):
    fig, ax = small_mesh.plot_label(label='edge')
    plt.close(fig)


def test_plot_label_element_only(small_mesh):
    fig, ax = small_mesh.plot_label(label='element')
    plt.close(fig)


# ---------------------------------------------------------------------------
# plot_layer
# ---------------------------------------------------------------------------

def test_plot_layer_returns_fig_ax(small_mesh):
    fig, ax = small_mesh.plot_layer()
    assert isinstance(fig, plt.Figure)
    assert isinstance(ax, plt.Axes)
    plt.close(fig)


def test_plot_layer_subset(small_mesh):
    layers = list(range(min(2, small_mesh.n_layers)))
    fig, ax = small_mesh.plot_layer(layers=layers)
    plt.close(fig)


# ---------------------------------------------------------------------------
# plot_quality
# ---------------------------------------------------------------------------

def test_plot_quality_returns_fig_ax(small_mesh):
    fig, ax = small_mesh.plot_quality()
    assert isinstance(fig, plt.Figure)
    assert isinstance(ax, plt.Axes)
    plt.close(fig)


def test_plot_quality_subset_elem_ids(small_mesh):
    ids = np.arange(min(5, small_mesh.n_elems))
    fig, ax = small_mesh.plot_quality(elem_ids=ids)
    plt.close(fig)


# ---------------------------------------------------------------------------
# axis_chilmesh helper
# ---------------------------------------------------------------------------

def test_axis_chilmesh_sets_equal_aspect(small_mesh):
    fig, ax = plt.subplots()
    small_mesh.axis_chilmesh(ax=ax)
    assert ax.get_aspect() == 'equal'
    plt.close(fig)


def test_axis_chilmesh_without_ax_uses_gca(small_mesh):
    fig, _ = plt.subplots()
    ax = small_mesh.axis_chilmesh()   # no ax supplied — uses gca
    assert ax.get_aspect() == 'equal'
    plt.close(fig)
