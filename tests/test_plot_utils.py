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


def _count_artists(ax):
    """Sum of artist counts across the categories CHILmesh plot methods write to."""
    return len(ax.collections) + len(ax.lines) + len(ax.patches)


def _assert_limits_enclose(ax, mesh, *, slack: float = 1e-6):
    """Axis limits must enclose every mesh vertex (with small slack)."""
    xs, ys = mesh.points[:, 0], mesh.points[:, 1]
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    assert xmin - slack <= xs.min() and xs.max() <= xmax + slack, (
        f"x-limits [{xmin}, {xmax}] do not enclose data [{xs.min()}, {xs.max()}]"
    )
    assert ymin - slack <= ys.min() and ys.max() <= ymax + slack, (
        f"y-limits [{ymin}, {ymax}] do not enclose data [{ys.min()}, {ys.max()}]"
    )


# ---------------------------------------------------------------------------
# Basic plot methods
# ---------------------------------------------------------------------------

def test_plot_returns_fig_ax(small_mesh):
    fig, ax = small_mesh.plot()
    assert isinstance(fig, plt.Figure)
    assert isinstance(ax, plt.Axes)
    assert _count_artists(ax) >= 1, "plot() drew nothing onto the axes"
    _assert_limits_enclose(ax, small_mesh)
    plt.close(fig)


def test_plot_with_elem_color(small_mesh):
    fig, ax = small_mesh.plot(elem_color='steelblue')
    assert isinstance(fig, plt.Figure)
    # Filled-element path uses PolyCollection; expect at least one collection.
    assert len(ax.collections) >= 1, (
        "plot(elem_color='steelblue') produced no PolyCollection"
    )
    _assert_limits_enclose(ax, small_mesh)
    plt.close(fig)


def test_plot_subset_elem_ids(small_mesh):
    ids = np.arange(min(5, small_mesh.n_elems))
    fig, ax = small_mesh.plot(elem_ids=ids)
    assert _count_artists(ax) >= 1, "plot(elem_ids=...) drew nothing"
    plt.close(fig)


def test_plot_reuses_provided_ax(small_mesh):
    fig, ax = plt.subplots()
    before = _count_artists(ax)
    _, returned_ax = small_mesh.plot(ax=ax)
    assert returned_ax is ax, "plot(ax=ax) returned a different axes"
    assert _count_artists(ax) > before, (
        f"plot(ax=ax) did not add any artists (before={before}, after={_count_artists(ax)})"
    )
    plt.close(fig)


# ---------------------------------------------------------------------------
# plot_edge
# ---------------------------------------------------------------------------

def test_plot_edge_returns_fig_ax(small_mesh):
    fig, ax = small_mesh.plot_edge()
    assert isinstance(fig, plt.Figure)
    assert isinstance(ax, plt.Axes)
    # Every plot_edge call uses a single LineCollection on the axes.
    assert len(ax.collections) >= 1, "plot_edge() produced no LineCollection"
    n_segments = sum(len(c.get_segments()) for c in ax.collections
                     if hasattr(c, "get_segments"))
    assert n_segments == small_mesh.n_edges, (
        f"plot_edge() drew {n_segments} segments, expected {small_mesh.n_edges}"
    )
    _assert_limits_enclose(ax, small_mesh)
    plt.close(fig)


def test_plot_edge_subset(small_mesh):
    n = min(10, small_mesh.n_edges)
    ids = np.arange(n)
    fig, ax = small_mesh.plot_edge(edge_ids=ids)
    n_segments = sum(len(c.get_segments()) for c in ax.collections
                     if hasattr(c, "get_segments"))
    assert n_segments == n, (
        f"plot_edge(edge_ids=[{n}]) drew {n_segments} segments, expected {n}"
    )
    plt.close(fig)


# ---------------------------------------------------------------------------
# plot_boundary and plot_interior_edges (new methods)
# ---------------------------------------------------------------------------

def test_plot_boundary_returns_fig_ax(small_mesh):
    fig, ax = small_mesh.plot_boundary()
    assert isinstance(fig, plt.Figure)
    assert isinstance(ax, plt.Axes)
    n_segments = sum(len(c.get_segments()) for c in ax.collections
                     if hasattr(c, "get_segments"))
    n_boundary = len(small_mesh.boundary_edges())
    assert n_segments == n_boundary, (
        f"plot_boundary() drew {n_segments} segments, expected {n_boundary}"
    )
    plt.close(fig)


def test_plot_boundary_non_empty(small_mesh):
    boundary = small_mesh.boundary_edges()
    assert len(boundary) > 0, "expected at least one boundary edge"


def test_plot_interior_edges_returns_fig_ax(small_mesh):
    fig, ax = small_mesh.plot_interior_edges()
    assert isinstance(fig, plt.Figure)
    n_segments = sum(len(c.get_segments()) for c in ax.collections
                     if hasattr(c, "get_segments"))
    expected = small_mesh.n_edges - len(small_mesh.boundary_edges())
    assert n_segments == expected, (
        f"plot_interior_edges() drew {n_segments} segments, expected {expected}"
    )
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
    # plot_elem calls plot() then adds a filled PolyCollection on top.
    assert len(ax.collections) >= 1, "plot_elem() produced no PolyCollection"
    _assert_limits_enclose(ax, small_mesh)
    plt.close(fig)


def test_plot_face_is_alias(small_mesh):
    fig_face, ax_face = small_mesh.plot_face()
    assert isinstance(fig_face, plt.Figure)
    # Mirror plot_elem; both should add a PolyCollection.
    assert len(ax_face.collections) >= 1, "plot_face() produced no PolyCollection"
    plt.close(fig_face)


# ---------------------------------------------------------------------------
# plot_point
# ---------------------------------------------------------------------------

def test_plot_point_vertex(small_mesh):
    fig, ax = small_mesh.plot_point(point_type='vertex')
    assert isinstance(fig, plt.Figure)
    # plot_point uses ax.plot(linestyle='none', marker='o'); one Line2D per call.
    point_lines = [ln for ln in ax.lines if ln.get_linestyle() == 'None']
    assert len(point_lines) >= 1, "plot_point(vertex) produced no marker Line2D"
    n_points = sum(len(ln.get_xdata()) for ln in point_lines)
    assert n_points == small_mesh.n_verts, (
        f"plot_point(vertex) drew {n_points} points, expected {small_mesh.n_verts}"
    )
    plt.close(fig)


def test_plot_point_edge_midpoint(small_mesh):
    fig, ax = small_mesh.plot_point(point_type='edge')
    point_lines = [ln for ln in ax.lines if ln.get_linestyle() == 'None']
    n_points = sum(len(ln.get_xdata()) for ln in point_lines)
    assert n_points == small_mesh.n_edges, (
        f"plot_point(edge) drew {n_points} midpoints, expected {small_mesh.n_edges}"
    )
    plt.close(fig)


def test_plot_point_centroid(small_mesh):
    fig, ax = small_mesh.plot_point(point_type='centroid')
    point_lines = [ln for ln in ax.lines if ln.get_linestyle() == 'None']
    n_points = sum(len(ln.get_xdata()) for ln in point_lines)
    assert n_points == small_mesh.n_elems, (
        f"plot_point(centroid) drew {n_points} centroids, expected {small_mesh.n_elems}"
    )
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
    # 'all' must annotate vertices, edges, and elements — at least n_verts texts.
    assert len(ax.texts) >= small_mesh.n_verts, (
        f"plot_label('all') produced {len(ax.texts)} text annotations, "
        f"expected at least {small_mesh.n_verts}"
    )
    plt.close(fig)


def test_plot_label_vertex_only(small_mesh):
    fig, ax = small_mesh.plot_label(label='vertex')
    assert len(ax.texts) == small_mesh.n_verts, (
        f"plot_label('vertex') produced {len(ax.texts)} texts, "
        f"expected {small_mesh.n_verts}"
    )
    plt.close(fig)


def test_plot_label_edge_only(small_mesh):
    fig, ax = small_mesh.plot_label(label='edge')
    assert len(ax.texts) == small_mesh.n_edges, (
        f"plot_label('edge') produced {len(ax.texts)} texts, "
        f"expected {small_mesh.n_edges}"
    )
    plt.close(fig)


def test_plot_label_element_only(small_mesh):
    fig, ax = small_mesh.plot_label(label='element')
    assert len(ax.texts) == small_mesh.n_elems, (
        f"plot_label('element') produced {len(ax.texts)} texts, "
        f"expected {small_mesh.n_elems}"
    )
    plt.close(fig)


# ---------------------------------------------------------------------------
# plot_layer
# ---------------------------------------------------------------------------

def test_plot_layer_returns_fig_ax(small_mesh):
    fig, ax = small_mesh.plot_layer()
    assert isinstance(fig, plt.Figure)
    assert isinstance(ax, plt.Axes)
    # plot_layer renders one PolyCollection per layer + adds a colorbar +
    # writes the title "Mesh Layers (N total)".
    assert len(ax.collections) >= small_mesh.n_layers, (
        f"plot_layer() produced {len(ax.collections)} collections, "
        f"expected at least {small_mesh.n_layers}"
    )
    title = ax.get_title()
    assert "Mesh Layers" in title, f"expected 'Mesh Layers' in title, got {title!r}"
    assert f"{small_mesh.n_layers}" in title, (
        f"expected layer count {small_mesh.n_layers} in title, got {title!r}"
    )
    _assert_limits_enclose(ax, small_mesh)
    plt.close(fig)


def test_plot_layer_subset(small_mesh):
    n = min(2, small_mesh.n_layers)
    layers = list(range(n))
    fig, ax = small_mesh.plot_layer(layers=layers)
    # Subset must produce at least one collection per requested layer.
    assert len(ax.collections) >= n, (
        f"plot_layer(layers={layers}) produced {len(ax.collections)} collections, "
        f"expected at least {n}"
    )
    plt.close(fig)


# ---------------------------------------------------------------------------
# plot_quality
# ---------------------------------------------------------------------------

def test_plot_quality_returns_fig_ax(small_mesh):
    fig, ax = small_mesh.plot_quality()
    assert isinstance(fig, plt.Figure)
    assert isinstance(ax, plt.Axes)
    # plot_quality writes the title "Element Quality" and adds at least one
    # filled PolyCollection per non-empty quality bin.
    assert ax.get_title() == "Element Quality", (
        f"expected title 'Element Quality', got {ax.get_title()!r}"
    )
    assert len(ax.collections) >= 1, "plot_quality() produced no PolyCollection"
    _assert_limits_enclose(ax, small_mesh)
    plt.close(fig)


def test_plot_quality_subset_elem_ids(small_mesh):
    ids = np.arange(min(5, small_mesh.n_elems))
    fig, ax = small_mesh.plot_quality(elem_ids=ids)
    assert len(ax.collections) >= 1, (
        "plot_quality(elem_ids=...) produced no PolyCollection"
    )
    plt.close(fig)


# ---------------------------------------------------------------------------
# plot_quality_histogram
# ---------------------------------------------------------------------------

def test_plot_quality_histogram_returns_fig_ax(small_mesh):
    fig, ax = small_mesh.plot_quality_histogram(bins=20)
    assert isinstance(fig, plt.Figure)
    assert isinstance(ax, plt.Axes)
    assert ax.get_title() == "Element Quality Distribution", (
        f"expected title 'Element Quality Distribution', got {ax.get_title()!r}"
    )
    assert ax.get_xlabel() == "Quality"
    assert ax.get_ylabel() == "# of Elements"
    # x-axis covers full [0, 1] quality range regardless of data extent
    xmin, xmax = ax.get_xlim()
    assert xmin == 0.0 and xmax == 1.0, f"xlim {xmin, xmax} != (0.0, 1.0)"
    # one Rectangle patch per bin
    assert len(ax.patches) == 20, (
        f"expected 20 bins, got {len(ax.patches)} Rectangle patches"
    )
    plt.close(fig)


def test_plot_quality_histogram_total_count_matches_subset(small_mesh):
    ids = np.arange(min(5, small_mesh.n_elems))
    fig, ax = small_mesh.plot_quality_histogram(elem_ids=ids, bins=10)
    total = sum(int(p.get_height()) for p in ax.patches)
    assert total == len(ids), (
        f"histogram total {total} != requested elem count {len(ids)}"
    )
    plt.close(fig)


def test_plot_quality_histogram_bin_colors_track_midpoints(small_mesh):
    """Each bar's face colour must equal the cool_r colormap value at its
    bin midpoint, matching the 2-D plot_quality colormap."""
    import matplotlib
    from matplotlib.colors import Normalize

    fig, ax = small_mesh.plot_quality_histogram(bins=5)
    cmap = matplotlib.colormaps['cool_r']
    norm = Normalize(vmin=0.0, vmax=1.0)
    expected_midpoints = np.linspace(0.1, 0.9, 5)  # bin centres for 5 bins on [0,1]
    expected_colors = cmap(norm(expected_midpoints))

    for patch, expected in zip(ax.patches, expected_colors):
        actual = np.array(patch.get_facecolor())
        np.testing.assert_allclose(
            actual, expected, atol=1e-6,
            err_msg=f"bar at x={patch.get_x():.2f} colour mismatch",
        )
    plt.close(fig)


# ---------------------------------------------------------------------------
# axis_chilmesh helper
# ---------------------------------------------------------------------------

def test_axis_chilmesh_sets_equal_aspect(small_mesh):
    fig, ax = plt.subplots()
    small_mesh.axis_chilmesh(ax=ax)
    # matplotlib normalises set_aspect('equal') to the numeric value 1.0
    assert ax.get_aspect() in ('equal', 1.0)
    plt.close(fig)


def test_axis_chilmesh_without_ax_uses_gca(small_mesh):
    fig, _ = plt.subplots()
    ax = small_mesh.axis_chilmesh()   # no ax supplied — uses gca
    assert ax.get_aspect() in ('equal', 1.0)
    plt.close(fig)
