"""cartogram: size-controlled mesh cartograms (CHILmesh #219). Headless (Agg)."""
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402

from chilmesh import CHILmesh, examples  # noqa: E402
from chilmesh.chilplotting import cartogram  # noqa: E402


@pytest.fixture(autouse=True)
def _close():
    yield
    plt.close("all")


def _annulus():
    return examples.annulus.__wrapped__()


def _pc(ax):
    return next(c for c in ax.collections if isinstance(c, PolyCollection))


def _area(vertices):
    v = np.asarray(vertices)
    if len(v) > 1 and np.allclose(v[0], v[-1]):
        v = v[:-1]  # closed path repeats the first vertex
    x, y = v[:, 0], v[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def _mixed_mesh():
    # 3x2 cells on a 4x3 grid; cell (1, 0) split into two padded triangles
    xs, ys = np.meshgrid(np.arange(4.0), np.arange(3.0))
    pts = np.c_[xs.ravel(), ys.ravel(), np.zeros(12)]

    def v(i, j):
        return j * 4 + i

    conn = []
    for j in range(2):
        for i in range(3):
            a, b, c, d = v(i, j), v(i + 1, j), v(i + 1, j + 1), v(i, j + 1)
            if (i, j) == (1, 0):
                conn += [[a, b, c, a], [a, c, d, a]]
            else:
                conn.append([a, b, c, d])
    return CHILmesh(np.array(conn), pts)


def _ring_mesh():
    # 2m triangles on 2m vertices: n_elems == n_verts
    n = 6
    ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
    unit = np.c_[np.cos(ang), np.sin(ang)]
    pts = np.c_[np.r_[unit, 2 * unit], np.zeros(2 * n)]
    conn = []
    for i in range(n):
        j = (i + 1) % n
        conn += [[i, j, n + i], [j, n + j, n + i]]
    return CHILmesh(np.array(conn), pts)


def _field(m):
    return m.elem_quality()[0]


@pytest.mark.parametrize("kind", ["noncontig", "unrolled"])
@pytest.mark.parametrize("reference", [None, "mean"])
def test_returns_figure_axes(kind, reference):
    m = _annulus()
    fig, ax = cartogram(m, _field(m), kind=kind, reference=reference)
    assert isinstance(fig, plt.Figure) and isinstance(ax, plt.Axes)
    assert len(fig.axes) == 2  # main axes + colorbar


def test_mixin_delegates_and_colorbar_label():
    m = _annulus()
    fig, ax = m.cartogram(_field(m), kind="unrolled", clabel="skew")
    assert fig.axes[1].get_ylabel() == "skew"
    fig, ax = m.cartogram(_field(m), colorbar=False)
    assert len(fig.axes) == 1


def test_noncontig_polygon_count_is_finite_element_count():
    m = _annulus()
    v = _field(m).copy()
    v[[3, 17, 40]] = np.nan
    _, ax = cartogram(m, v)
    assert len(_pc(ax).get_paths()) == m.n_elems - 3


def test_noncontig_draw_order_ascending_by_value():
    m = _annulus()
    v = _field(m)
    _, ax = cartogram(m, v)
    arr = _pc(ax).get_array()
    assert np.all(np.diff(arr) >= 0)


def test_noncontig_equal_area_when_reference_none():
    m = _annulus()
    _, ax = cartogram(m, _field(m))
    areas = np.array([_area(p.vertices) for p in _pc(ax).get_paths()])
    assert areas.size == m.n_elems
    np.testing.assert_allclose(areas, areas[0], rtol=1e-9)
    verts = np.asarray(m.points)[np.asarray(m.connectivity_list)[:, :3], :2]
    med = np.median([_area(t) for t in verts])
    np.testing.assert_allclose(areas[0], med, rtol=1e-9)


def test_noncontig_focus_context_changes_sizes_and_clips_focus_factor():
    m = _annulus()
    v = _field(m)
    _, ax = cartogram(m, v, reference="mean", scale_clip=(0.5, 2.0))
    areas = np.array([_area(p.vertices) for p in _pc(ax).get_paths()])
    assert areas.max() / areas.min() > 1.5  # sizes now differ
    # focus factor in [0.5, 2] -> area factor in [0.25, 4] around the median area
    verts = np.asarray(m.points)[np.asarray(m.connectivity_list)[:, :3], :2]
    base_area = np.median([_area(t) for t in verts])
    assert areas.max() <= 4.0 * base_area * (1 + 1e-9)
    assert areas.min() >= 0.25 * base_area * (1 - 1e-9)


def test_unrolled_every_element_in_exactly_one_band():
    m = _annulus()
    layers = [np.unique(np.concatenate((m.layers["OE"][i], m.layers["IE"][i])))
              for i in range(m.n_layers)]
    allids = np.concatenate(layers)
    assert len(allids) == m.n_elems and len(np.unique(allids)) == m.n_elems
    _, ax = cartogram(m, _field(m), kind="unrolled")
    paths = _pc(ax).get_paths()
    assert len(paths) == m.n_elems
    ys = np.array([p.vertices[:, 1].min() for p in paths])
    bands = np.unique(ys)
    assert len(bands) == m.n_layers
    counts = [int((ys == y).sum()) for y in bands[::-1]]  # layer 0 is at the top
    assert counts == [len(b) for b in layers]
    assert sum(counts) == m.n_elems


def test_unrolled_equal_width_without_reference_and_proportional_with():
    m = _annulus()
    v = _field(m)
    _, ax = cartogram(m, v, kind="unrolled")
    widths = np.array([np.ptp(p.vertices[:, 0]) for p in _pc(ax).get_paths()])
    ys = np.array([p.vertices[:, 1].min() for p in _pc(ax).get_paths()])
    for y in np.unique(ys):
        w = widths[ys == y]
        np.testing.assert_allclose(w, 1.0 / len(w), rtol=1e-9)
        assert w.sum() == pytest.approx(1.0)
    _, ax2 = cartogram(m, v, kind="unrolled", reference="mean")
    w2 = np.array([np.ptp(p.vertices[:, 0]) for p in _pc(ax2).get_paths()])
    assert w2.max() / w2.min() > 2
    # widths follow |v - ref| + eps within a band (cells are in band-angle order)
    ref = float(np.mean(v))
    eps = 0.004 * (v.max() - v.min())
    layer0 = np.unique(np.concatenate((m.layers["OE"][0], m.layers["IE"][0])))
    arr = _pc(ax2).get_array()
    top = np.isclose(np.array([p.vertices[:, 1].min() for p in _pc(ax2).get_paths()]),
                     m.n_layers - 1)
    expect = np.abs(arr[top] - ref) + eps
    np.testing.assert_allclose(w2[top], expect / expect.sum(), rtol=1e-9)
    assert top.sum() == len(layer0)


def test_per_vertex_input_matches_mean_of_nodes():
    m = _annulus()
    pv = m.points[:, 0].astype(float)
    _, ax = cartogram(m, pv, kind="unrolled")
    conn = np.asarray(m.connectivity_list)
    expect = pv[conn].mean(axis=1)
    got = np.sort(_pc(ax).get_array())
    np.testing.assert_allclose(got, np.sort(expect))


def test_per_vertex_nan_node_blanks_element():
    m = _annulus()
    pv = m.points[:, 0].astype(float).copy()
    pv[0] = np.nan
    touching = np.array([0 in set(map(int, r)) for r in m.connectivity_list])
    _, ax = cartogram(m, pv)
    assert len(_pc(ax).get_paths()) == m.n_elems - touching.sum()


def test_length_mismatch_raises_naming_both_counts():
    m = _annulus()
    with pytest.raises(ValueError) as ei:
        cartogram(m, np.zeros(m.n_elems + 7))
    assert str(m.n_elems) in str(ei.value) and str(m.n_verts) in str(ei.value)
    with pytest.raises(ValueError):
        cartogram(m, np.zeros(m.n_elems - 1), on="elem")
    with pytest.raises(ValueError):
        cartogram(m, np.zeros(m.n_elems), on="vertex")
    with pytest.raises(ValueError):
        cartogram(m, np.zeros(m.n_elems), on="face")


def test_equal_elem_and_vert_counts_require_on():
    m = _ring_mesh()
    assert m.n_elems == m.n_verts
    v = np.linspace(0, 1, m.n_elems)
    with pytest.raises(ValueError) as ei:
        cartogram(m, v)
    assert "on=" in str(ei.value) and str(m.n_elems) in str(ei.value)
    for on in ("elem", "vertex"):
        fig, ax = cartogram(m, v, on=on)
        assert len(_pc(ax).get_paths()) == m.n_elems


def test_mode_ties_pick_smallest():
    from chilmesh.chilplotting import _cartogram_reference
    # 5 and 1 tie on count (2 each); 1 is the smallest
    assert _cartogram_reference(np.array([5.0, 1.0, 5.0, 1.0, 3.0]), "mode") == 1.0
    assert _cartogram_reference(np.array([2.0, 2.0, 7.0]), "mode") == 2.0
    # exact values: 0.1 + 0.2 and 0.3 are different values
    assert _cartogram_reference(np.array([0.3, 0.1 + 0.2, 0.1 + 0.2]), "mode") == 0.1 + 0.2


def test_mode_reference_end_to_end_picks_smallest_tied_value():
    m = _annulus()
    n = m.n_elems
    v = np.arange(n, dtype=float) + 10.0          # all distinct
    v[:100] = 5.0                                  # tie: 100 x 5.0 ...
    v[100:200] = 1.0                               # ... and 100 x 1.0
    _, ax = cartogram(m, v, kind="unrolled", reference="mode")
    assert "|value - 1|" in ax.get_xlabel()


def test_reference_values():
    from chilmesh.chilplotting import _cartogram_reference
    f = np.array([1.0, 2.0, 2.0, 9.0])
    assert _cartogram_reference(f, "mean") == pytest.approx(3.5)
    assert _cartogram_reference(f, "median") == 2.0
    assert _cartogram_reference(f, "max") == 9.0
    assert _cartogram_reference(f, "min") == 1.0
    assert _cartogram_reference(f, 4) == 4.0
    assert _cartogram_reference(f, None) is None
    with pytest.raises(ValueError):
        _cartogram_reference(f, "modal")
    with pytest.raises(ValueError):
        _cartogram_reference(f, [1.0])


def test_nan_element_not_drawn_and_excluded_from_colour_limits():
    m = _annulus()
    v = _field(m).copy()
    hi = int(np.argmax(v))
    lo = int(np.argmin(v))
    fin = v.copy()
    fin[[hi, lo]] = np.nan
    _, ax = cartogram(m, fin)
    pc = _pc(ax)
    assert len(pc.get_paths()) == m.n_elems - 2
    assert pc.get_clim() == (np.nanmin(fin), np.nanmax(fin))


def test_unrolled_nan_keeps_empty_slot_with_eps_width_under_reference():
    m = _annulus()
    v = _field(m).copy()
    v[5] = np.nan
    fin = v[np.isfinite(v)]
    eps = 0.004 * (fin.max() - fin.min())
    _, ax = cartogram(m, v, kind="unrolled", reference="mean")
    pc = _pc(ax)
    arr = np.ma.masked_invalid(pc.get_array())
    assert len(pc.get_paths()) == m.n_elems  # slot kept
    assert arr.mask.sum() == 1               # ... and blank
    w = np.array([np.ptp(p.vertices[:, 0]) for p in pc.get_paths()])
    ys = np.array([p.vertices[:, 1].min() for p in pc.get_paths()])
    slot = int(np.flatnonzero(arr.mask)[0])
    band = ys == ys[slot]
    ref = float(fin.mean())
    weights = np.where(arr.mask[band], eps, np.abs(arr.filled(0)[band] - ref) + eps)
    np.testing.assert_allclose(w[band], weights / weights.sum(), rtol=1e-9)
    assert pc.get_clim() == (fin.min(), fin.max())


def test_eps_override_is_absolute():
    m = _annulus()
    v = _field(m)
    _, a = cartogram(m, v, kind="unrolled", reference="mean", eps=1e3)
    w = np.array([np.ptp(p.vertices[:, 0]) for p in _pc(a).get_paths()])
    ys = np.array([p.vertices[:, 1].min() for p in _pc(a).get_paths()])
    # a huge eps swamps |v - ref|, so cells in a band are almost equal width
    for y in np.unique(ys):
        sel = w[ys == y]
        assert sel.max() / sel.min() < 1.01


def test_unknown_kind_raises():
    m = _annulus()
    with pytest.raises(ValueError):
        cartogram(m, _field(m), kind="hybrid")


def test_unknown_reference_raises():
    m = _annulus()
    with pytest.raises(ValueError):
        cartogram(m, _field(m), reference="modal")


@pytest.mark.parametrize("kind", ["noncontig", "unrolled"])
@pytest.mark.parametrize("reference", [None, "mean", "mode"])
def test_mixed_tri_quad_mesh(kind, reference):
    m = _mixed_mesh()
    conn = np.asarray(m.connectivity_list)
    assert conn.shape[1] == 4 and (conn[:, 3] == conn[:, 0]).any()
    v = np.arange(m.n_elems, dtype=float)
    fig, ax = cartogram(m, v, kind=kind, reference=reference)
    pc = _pc(ax)
    assert len(pc.get_paths()) == m.n_elems
    if kind == "noncontig" and reference is None:
        areas = np.array([_area(p.vertices) for p in pc.get_paths()])
        np.testing.assert_allclose(areas, areas[0], rtol=1e-9)
    # per-vertex input on a mixed mesh
    cartogram(m, np.arange(m.n_verts, dtype=float), kind=kind, reference=reference)


@pytest.mark.parametrize("kind", ["noncontig", "unrolled"])
@pytest.mark.parametrize("reference", [None, "mean", "median", "mode", "max", "min", 0.5])
def test_constant_field_falls_back_to_equal_size(kind, reference):
    m = _annulus()
    fig, ax = cartogram(m, np.full(m.n_elems, 0.7), kind=kind, reference=reference)
    pc = _pc(ax)
    assert len(pc.get_paths()) == m.n_elems
    if kind == "noncontig":
        areas = np.array([_area(p.vertices) for p in pc.get_paths()])
        np.testing.assert_allclose(areas, areas[0], rtol=1e-9)
    else:
        w = np.array([np.ptp(p.vertices[:, 0]) for p in pc.get_paths()])
        ys = np.array([p.vertices[:, 1].min() for p in pc.get_paths()])
        for y in np.unique(ys):
            sel = w[ys == y]
            np.testing.assert_allclose(sel, 1.0 / len(sel), rtol=1e-9)


@pytest.mark.parametrize("kind", ["noncontig", "unrolled"])
def test_all_nan_does_not_crash(kind):
    m = _annulus()
    fig, ax = cartogram(m, np.full(m.n_elems, np.nan), kind=kind, reference="mean")
    assert fig is not None


def test_uses_given_axes_and_cmap_and_no_edges():
    m = _annulus()
    fig0, ax0 = plt.subplots()
    fig, ax = cartogram(m, _field(m), ax=ax0, cmap="plasma")
    assert ax is ax0 and fig is fig0
    pc = _pc(ax)
    assert pc.get_cmap().name == "plasma"
    assert pc.get_linewidths() is not None
    assert pc.get_edgecolor().size == 0 or pc.get_edgecolor()[0][3] == 0


def test_bad_scale_clip_and_eps_raise():
    m = _annulus()
    with pytest.raises(ValueError):
        cartogram(m, _field(m), scale_clip=(2.0, 1.0))
    with pytest.raises(ValueError):
        cartogram(m, _field(m), eps=-1.0)
