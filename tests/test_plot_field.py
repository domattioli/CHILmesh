"""plot_field: node-aligned field on a mesh (spec 005 US6). Headless (Agg)."""
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402

from chilmesh import examples  # noqa: E402
from chilmesh.chilplotting import plot_field  # noqa: E402

# annulus/donut are triangular, quad_2x2 quad-only, block_o and structured mixed/other
CASES = ["annulus", "quad_2x2", "structured", "donut"]


@pytest.fixture(autouse=True)
def _close():
    yield
    plt.close("all")


def _mesh(name):
    return getattr(examples, name).__wrapped__()


def _pc(ax):
    return next(c for c in ax.collections if isinstance(c, PolyCollection))


@pytest.mark.parametrize("name", CASES)
def test_colorbar_present_and_one_polygon_per_element(name):
    m = _mesh(name)
    fig, ax = plot_field(m, m.points[:, 0])
    assert len(fig.axes) == 2  # main axes + colorbar axes
    assert len(_pc(ax).get_paths()) == m.n_elems


def test_no_colorbar_when_disabled():
    m = _mesh("annulus")
    fig, _ = plot_field(m, m.points[:, 0], colorbar=False)
    assert len(fig.axes) == 1


def test_mixed_mesh_covered():
    # a mixed tri/quad mesh built from the quad fixture plus one padded triangle
    m = _mesh("block_o")
    fig, ax = plot_field(m, m.points[:, 1])
    assert len(_pc(ax).get_paths()) == m.n_elems


def test_nan_nodes_are_blank():
    m = _mesh("annulus")
    vals = m.points[:, 0].astype(float).copy()
    vals[0] = np.nan
    fig, ax = plot_field(m, vals)
    fig.canvas.draw()
    pc = _pc(ax)
    fc = pc.get_facecolor()
    touching = np.array([0 in set(map(int, row)) for row in m.connectivity_list])
    assert touching.any()
    assert (fc[touching][:, 3] == 0).all()      # blank (transparent)
    assert (fc[~touching][:, 3] > 0).all()      # others drawn


def test_all_nan_does_not_crash():
    m = _mesh("quad_2x2")
    fig, ax = plot_field(m, np.full(m.n_verts, np.nan))
    assert fig is not None


def test_length_mismatch_message():
    m = _mesh("annulus")
    with pytest.raises(ValueError) as ei:
        plot_field(m, np.zeros(m.n_verts - 1))
    assert str(m.n_verts - 1) in str(ei.value) and str(m.n_verts) in str(ei.value)


def test_uses_given_axes_and_cmap():
    m = _mesh("annulus")
    fig, ax0 = plt.subplots()
    fig2, ax = plot_field(m, m.points[:, 0], ax=ax0, cmap="plasma")
    assert ax is ax0 and fig2 is fig
    assert _pc(ax).get_cmap().name == "plasma"


def test_accepts_output_series_row(tmp_path):
    m = _mesh("annulus")
    vals = np.linspace(0, 1, m.n_verts)[None, :]
    fig, ax = plot_field(m, vals[0])
    assert _pc(ax).get_array().size == m.n_elems
