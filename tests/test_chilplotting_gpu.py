"""Tests for the plot render-backend selector and the GPU backend (CHILmesh #167).

Backend-selection tests run everywhere. Tests marked ``needs_gpu`` skip cleanly
when ``pygfx``/``wgpu`` or a GPU adapter is missing, so CI without a GPU stays
green; on a machine with an adapter they really render.
"""
import matplotlib
matplotlib.use("Agg")

import sys

import numpy as np
import pytest

import chilmesh
from chilmesh import chilplotting as cp

try:
    from chilmesh import chilplotting_gpu as gpu
    GPU_INFO = gpu.gpu_adapter_info()
except Exception:  # pragma: no cover - module import itself must never fail
    gpu, GPU_INFO = None, None

needs_gpu = pytest.mark.skipif(
    GPU_INFO is None, reason="pygfx/wgpu or a GPU adapter is not available")

SQUARE_PTS = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
MIXED_PTS = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [2, 0], [2, 1.0]])
# one quad plus one padded triangle (repeated last vertex)
MIXED_CONN = np.array([[0, 1, 2, 3], [1, 4, 5, 5]])


def _block_o():
    return chilmesh.examples.block_o()


# -- backend selection (no GPU needed) ------------------------------------

def test_default_backend_is_matplotlib(monkeypatch):
    monkeypatch.delenv("CHILMESH_PLOT_BACKEND", raising=False)
    info = cp.plot_backend_info()
    assert info["selected"] == "mpl"
    assert "mpl" in info["available"]
    assert set(info) == {"available", "selected", "versions"}


def test_env_mpl_override(monkeypatch):
    monkeypatch.setenv("CHILMESH_PLOT_BACKEND", "mpl")
    assert cp.plot_backend_info()["selected"] == "mpl"


def test_env_gpu_ignored_when_unavailable(monkeypatch):
    monkeypatch.setenv("CHILMESH_PLOT_BACKEND", "gpu")
    monkeypatch.setattr(cp, "find_spec", lambda name: None)
    info = cp.plot_backend_info()
    assert info["available"] == ["mpl"]
    assert info["selected"] == "mpl"


def test_import_chilmesh_does_not_import_gpu_stack():
    import subprocess
    out = subprocess.run(
        [sys.executable, "-c",
         "import sys, chilmesh; from chilmesh import chilplotting;"
         "print(any(m.split('.')[0] in ('pygfx','wgpu') for m in sys.modules))"],
        capture_output=True, text=True, check=True).stdout.strip()
    assert out == "False"


def test_render_image_mpl_shape_and_content():
    tris = np.array([[0, 1, 2], [0, 2, 3]])
    img = cp.render_image(SQUARE_PTS, tris, values=np.array([0.0, 1.0]),
                          backend="mpl", size=(200, 100))
    assert img.shape == (100, 200, 4) and img.dtype == np.uint8
    assert (img[..., :3] < 250).any()


def test_render_image_rejects_unknown_backend():
    with pytest.raises(ValueError):
        cp.render_image(SQUARE_PTS, np.array([[0, 1, 2]]), backend="vulkan")


def test_gpu_backend_falls_back_with_warning_when_unavailable(monkeypatch):
    if gpu is None:
        pytest.skip("chilplotting_gpu failed to import")
    monkeypatch.setattr(gpu, "gpu_available", lambda: False)
    with pytest.warns(RuntimeWarning, match="GPU plot backend unavailable"):
        img = cp.render_image(SQUARE_PTS, np.array([[0, 1, 2]]), backend="gpu",
                              size=(64, 64))
    assert img.shape == (64, 64, 4)


# -- real GPU rendering ----------------------------------------------------

def _mask(img, bg=250):
    return (img[..., :3] < bg).any(axis=2)


@needs_gpu
def test_gpu_triangulation_preserves_element_ids():
    tris, elem = gpu._triangulate(MIXED_CONN)
    assert elem.tolist() == [0, 0, 1]          # quad -> 2 tris, padded tri -> 1
    assert tris[2].tolist() == [1, 4, 5]
    assert (np.diff(elem) >= 0).all()


@needs_gpu
def test_gpu_edges_match_cpu_unique_edges():
    conn = np.vstack([MIXED_CONN])
    np.testing.assert_array_equal(gpu._unique_edges_vec(conn), cp.unique_edges(conn))


@needs_gpu
def test_gpu_render_nonempty_and_inside_bbox():
    m = _block_o()
    view = gpu.GPUMeshView(m.points, m.connectivity_list, size=(640, 480))
    img = view.render()
    assert img.shape == (480, 640, 4) and img.dtype == np.uint8
    mask = _mask(img)
    assert mask.any()
    # fit() applies equal aspect: expected pixel bbox from the data bounds
    x0, x1, y0, y1 = view._bounds
    scale = max((x1 - x0) * 1.02 / 640, (y1 - y0) * 1.02 / 480)  # data per pixel
    w_px, h_px = (x1 - x0) / scale, (y1 - y0) / scale
    cols, rows = np.where(mask.any(axis=0))[0], np.where(mask.any(axis=1))[0]
    assert abs((cols.min() + cols.max()) / 2 - 320) < 4
    assert abs((rows.min() + rows.max()) / 2 - 240) < 4
    assert abs((cols.max() - cols.min()) - w_px) < 6
    assert abs((rows.max() - rows.min()) - h_px) < 6


@needs_gpu
def test_gpu_element_colors_follow_scalar():
    # two disjoint quads; low scalar -> viridis dark end, high -> bright end
    pts = np.array([[0, 0], [1, 0], [1, 1], [0, 1],
                    [2, 0], [3, 0], [3, 1], [2, 1.0]])
    conn = np.array([[0, 1, 2, 3], [4, 5, 6, 7]])
    view = gpu.GPUMeshView(pts, conn, values=np.array([0.0, 1.0]),
                           edge_color=None, size=(300, 100))
    img = view.render().astype(float)
    left, right = img[50, 50, :3], img[50, 250, :3]
    # matplotlib viridis(0) ~ (68,1,84), viridis(1) ~ (253,231,37)
    assert np.abs(left - np.array([68, 1, 84])).max() < 12
    assert np.abs(right - np.array([253, 231, 37])).max() < 12
    # both triangles of a quad share one colour (element id preserved per triangle)
    assert np.abs(img[30, 70, :3] - img[70, 30, :3]).max() < 3


@needs_gpu
def test_gpu_set_scalar_updates_colors_without_geometry_change():
    pts = np.array([[0, 0], [1, 0], [1, 1], [0, 1.0]])
    view = gpu.GPUMeshView(pts, np.array([[0, 1, 2, 3]]), values=np.array([0.0]),
                           vmin=0.0, vmax=1.0, edge_color=None, size=(100, 100))
    before = view.render()[50, 50, :3].astype(float)
    pos_id = id(view.mesh.geometry.positions)
    view.set_scalar(np.array([1.0]))
    after = view.render()[50, 50, :3].astype(float)
    assert np.abs(before - np.array([68, 1, 84])).max() < 12
    assert np.abs(after - np.array([253, 231, 37])).max() < 12
    assert id(view.mesh.geometry.positions) == pos_id
    with pytest.raises(ValueError):
        view.set_scalar(np.array([1.0, 2.0]))


@needs_gpu
def test_gpu_vs_matplotlib_coverage_iou():
    m = _block_o()
    kw = dict(values=np.arange(len(m.connectivity_list), dtype=float),
              edge_color=None, size=(640, 480))
    g = cp.render_image(m.points, m.connectivity_list, backend="gpu", **kw)
    c = cp.render_image(m.points, m.connectivity_list, backend="mpl", **kw)
    mg, mc = _mask(g), _mask(c)
    iou = (mg & mc).sum() / (mg | mc).sum()
    assert iou >= 0.9, f"coverage IoU {iou:.3f}"


@needs_gpu
def test_gpu_mixed_connectivity_renders():
    img = gpu.render_offscreen(MIXED_PTS, MIXED_CONN, size=(200, 100))
    assert _mask(img).any()


@needs_gpu
def test_gpu_rejects_bad_input():
    with pytest.raises(ValueError):
        gpu.GPUMeshView(SQUARE_PTS, np.array([[0, 1, 9]]))
    with pytest.raises(ValueError):
        gpu.GPUMeshView(SQUARE_PTS, np.empty((0, 3), dtype=int))
