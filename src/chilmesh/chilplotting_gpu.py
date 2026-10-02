"""GPU render backend for :mod:`chilmesh.chilplotting` (CHILmesh #167, Phase A).

Renders a mesh from bare ``points`` + ``connectivity`` arrays through
``pygfx`` / ``wgpu`` (Metal, Vulkan, D3D12 or GL, whichever ``wgpu`` finds).
This module is imported lazily by ``chilplotting``; importing ``chilmesh``
never imports ``pygfx`` or ``wgpu``. Install the dependencies with
``pip install chilmesh[gpu]``.

Design notes
------------
* Quads are split into two triangles ``(0, 1, 2)`` and ``(0, 2, 3)``. Padded
  triangles in 4-column connectivity keep their first three vertices, the same
  convention ``chilplotting.build_polygons`` uses. Triangles stay grouped by
  element id and ``elem_of_tri`` maps each triangle back to its element.
* Every triangle owns three private vertices, so a per-element scalar becomes a
  per-vertex texture coordinate sampled through a 1-D colormap texture on the
  GPU. This costs 3 vertices per triangle but makes
  :meth:`GPUMeshView.set_scalar` an O(n_elements) buffer write with no geometry
  rebuild. Position and index buffers are uploaded once.
* Colormaps come from matplotlib (already a hard dependency), so colours match
  the matplotlib backend for any registered colormap name.
"""
from __future__ import annotations

from typing import Optional, Sequence, Tuple

import numpy as np

try:  # the optional [gpu] extra
    import wgpu
    import pygfx as gfx
    from rendercanvas.offscreen import OffscreenRenderCanvas
    _IMPORT_ERROR: Optional[Exception] = None
except Exception as exc:  # noqa: BLE001 - any failure means "GPU unavailable"
    wgpu = gfx = OffscreenRenderCanvas = None  # type: ignore[assignment]
    _IMPORT_ERROR = exc

__all__ = ["GPUMeshView", "gpu_adapter_info", "gpu_available", "render_offscreen"]


def gpu_adapter_info() -> Optional[dict]:
    """Return a summary of the wgpu adapter, or ``None`` if none is available."""
    if _IMPORT_ERROR is not None:
        return None
    try:
        adapter = wgpu.gpu.request_adapter_sync(power_preference="high-performance")
    except Exception:  # noqa: BLE001 - no adapter on headless CI
        return None
    if adapter is None:
        return None
    info = dict(adapter.info)
    info["summary"] = adapter.summary
    info["wgpu_version"] = wgpu.__version__
    info["pygfx_version"] = gfx.__version__
    return info


def gpu_available() -> bool:
    """Return True when pygfx/wgpu import and a GPU adapter can be obtained."""
    return gpu_adapter_info() is not None


def _require() -> None:
    if _IMPORT_ERROR is not None:
        raise ImportError(
            "The GPU plot backend needs pygfx and wgpu: pip install chilmesh[gpu] "
            f"(import failed with: {_IMPORT_ERROR})")


def _triangulate(connectivity: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Split elements into triangles; return ``(tris, elem_of_tri)``."""
    from .chilplotting import _is_padded_or_degenerate

    conn = np.asarray(connectivity)
    if conn.ndim != 2 or conn.shape[1] not in (3, 4):
        raise ValueError(f"connectivity must have shape (n, 3|4), got {conn.shape}")
    if conn.shape[1] == 3:
        return conn.astype(np.int64), np.arange(len(conn))
    tri_mask = _is_padded_or_degenerate(conn)
    ids = np.arange(len(conn))
    quad_ids = ids[~tri_mask]
    q = conn[quad_ids]
    tris = np.concatenate([conn[tri_mask, :3], q[:, [0, 1, 2]], q[:, [0, 2, 3]]])
    elem = np.concatenate([ids[tri_mask], quad_ids, quad_ids])
    order = np.argsort(elem, kind="stable")
    return tris[order].astype(np.int64), elem[order]


def _unique_edges_vec(connectivity: np.ndarray) -> np.ndarray:
    """Vectorised unique undirected edges (same set as ``chilplotting.unique_edges``).

    A padded triangle repeats one vertex; its zero-length self-edge is dropped,
    leaving exactly the triangle's three edges.
    """
    conn = np.asarray(connectivity)
    e = np.stack([conn.ravel(), np.roll(conn, -1, axis=1).ravel()], axis=1)
    e = e[e[:, 0] != e[:, 1]]
    return np.unique(np.sort(e, axis=1), axis=0)


def _lut(cmap: str, n: int = 256) -> np.ndarray:
    import matplotlib
    return np.asarray(matplotlib.colormaps[cmap](np.linspace(0, 1, n)), dtype=np.float32)


class GPUMeshView:
    """A mesh uploaded to the GPU once, drawable many times.

    Parameters
    ----------
    points : ndarray, shape (n, 2) or (n, 3)
        Vertex coordinates (only x, y are used).
    connectivity : ndarray, shape (m, 3) or (m, 4)
        Element vertex ids; mixed meshes use padded triangles.
    values : ndarray, shape (m,), optional
        One scalar per element, mapped through ``cmap``. Without it elements
        use ``facecolor``.
    cmap : str
        matplotlib colormap name.
    vmin, vmax : float, optional
        Colour limits; default to the data range at construction/``set_scalar``.
    facecolor : str
        Fill colour when ``values`` is None.
    edge_color : str or None
        Edge colour; ``None`` hides edges.
    linewidth : float
        Edge width in screen pixels.
    size : (int, int)
        Offscreen canvas size ``(width, height)`` in pixels.
    canvas : optional
        An existing ``rendercanvas`` canvas; default is an offscreen canvas.
    background : str
        Background colour.
    """

    def __init__(self, points: np.ndarray, connectivity: np.ndarray, *,
                 values: Optional[np.ndarray] = None, cmap: str = "viridis",
                 vmin: Optional[float] = None, vmax: Optional[float] = None,
                 facecolor: str = "#1f77b4", edge_color: Optional[str] = "k",
                 linewidth: float = 1.0, size: Tuple[int, int] = (800, 600),
                 canvas=None, background: str = "white") -> None:
        _require()
        pts = np.asarray(points, dtype=np.float64)
        conn = np.asarray(connectivity)
        if pts.ndim != 2 or pts.shape[1] < 2:
            raise ValueError("points must have shape (n, >=2)")
        if len(conn) == 0:
            raise ValueError("connectivity is empty")
        if conn.min() < 0 or conn.max() >= len(pts):
            raise ValueError("connectivity references vertices outside points")

        self._n_elem = len(conn)
        self._vmin, self._vmax = vmin, vmax
        self._bounds = (pts[:, 0].min(), pts[:, 0].max(),
                        pts[:, 1].min(), pts[:, 1].max())
        # Shift to the bbox centre in float64 so float32 GPU buffers keep their
        # precision on large-coordinate meshes (UTM, state plane).
        self._origin = np.array([(self._bounds[0] + self._bounds[1]) / 2,
                                 (self._bounds[2] + self._bounds[3]) / 2])

        tris, self._elem_of_tri = _triangulate(conn)
        xy = (pts[:, :2] - self._origin).astype(np.float32)
        positions = np.zeros((len(tris) * 3, 3), dtype=np.float32)
        positions[:, :2] = xy[tris.ravel()]
        indices = np.arange(len(tris) * 3, dtype=np.uint32).reshape(-1, 3)

        self.canvas = canvas if canvas is not None else OffscreenRenderCanvas(
            size=size, pixel_ratio=1)
        self.renderer = gfx.renderers.WgpuRenderer(self.canvas)
        self.scene = gfx.Scene()
        self.scene.add(gfx.Background.from_color(background))

        geo = gfx.Geometry(positions=positions, indices=indices)
        self._scalar_mode = values is not None
        if self._scalar_mode:
            geo.texcoords = gfx.Buffer(np.zeros(len(tris) * 3, dtype=np.float32))
            mat = gfx.MeshBasicMaterial(
                map=gfx.TextureMap(gfx.Texture(_lut(cmap), dim=1),
                                   filter="nearest", wrap="clamp"))
        else:
            mat = gfx.MeshBasicMaterial(color=facecolor)
        self.mesh = gfx.Mesh(geo, mat)
        self.mesh.render_order = 0
        self.scene.add(self.mesh)
        if self._scalar_mode:
            self.set_scalar(values)

        self.edges = None
        if edge_color is not None:
            e = _unique_edges_vec(conn)
            seg = np.zeros((len(e) * 2, 3), dtype=np.float32)
            seg[:, :2] = xy[e.ravel()]
            self.edges = gfx.Line(
                gfx.Geometry(positions=seg),
                gfx.LineSegmentMaterial(thickness=linewidth, color=edge_color,
                                        thickness_space="screen", aa=True))
            self.edges.render_order = 1
            self.scene.add(self.edges)

        self.camera = gfx.OrthographicCamera()
        self.camera.maintain_aspect = False
        self.fit()

    def set_scalar(self, values: np.ndarray) -> None:
        """Replace the per-element scalar field without touching geometry."""
        if not self._scalar_mode:
            raise RuntimeError("view was built without values; pass values= at construction")
        v = np.asarray(values, dtype=np.float64).ravel()
        if len(v) != self._n_elem:
            raise ValueError(f"values length {len(v)} != number of elements {self._n_elem}")
        lo = float(v.min()) if self._vmin is None else self._vmin
        hi = float(v.max()) if self._vmax is None else self._vmax
        if hi <= lo:
            hi = lo + 1.0
        tc = np.clip((v[self._elem_of_tri] - lo) / (hi - lo), 0.0, 1.0)
        buf = self.mesh.geometry.texcoords
        buf.data[:] = np.repeat(tc, 3).astype(np.float32)
        buf.update_full()

    def set_view(self, xlim: Sequence[float], ylim: Sequence[float]) -> None:
        """Show exactly ``xlim`` by ``ylim`` in data units."""
        ox, oy = self._origin
        self.camera.show_rect(xlim[0] - ox, xlim[1] - ox, ylim[0] - oy, ylim[1] - oy)
        self.camera.maintain_aspect = False

    def fit(self, pad_frac: float = 0.01) -> None:
        """Fit the whole mesh with equal aspect, like ``configure_axes``."""
        x0, x1, y0, y1 = self._bounds
        w, h = self.canvas.get_logical_size()
        pad = pad_frac * max(x1 - x0, y1 - y0)
        scale = max((x1 - x0 + 2 * pad) / w, (y1 - y0 + 2 * pad) / h)  # data per pixel
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        self.set_view((cx - scale * w / 2, cx + scale * w / 2),
                      (cy - scale * h / 2, cy + scale * h / 2))

    def pan(self, dx: float, dy: float) -> None:
        """Translate the view by ``(dx, dy)`` data units."""
        self.camera.local.x += dx
        self.camera.local.y += dy

    def zoom(self, factor: float) -> None:
        """Zoom about the view centre; ``factor`` > 1 zooms in."""
        self.camera.zoom = self.camera.zoom * factor

    def render(self) -> np.ndarray:
        """Draw one frame; return an ``(h, w, 4)`` uint8 RGBA array."""
        self.renderer.render(self.scene, self.camera)
        return np.asarray(self.canvas.draw())

    def show(self) -> None:  # pragma: no cover - needs a display
        """Open an interactive window (needs a windowed rendercanvas backend)."""
        from rendercanvas.auto import RenderCanvas, loop
        canvas = RenderCanvas(size=self.canvas.get_logical_size(), title="CHILmesh")
        self.canvas = canvas
        self.renderer = gfx.renderers.WgpuRenderer(canvas)
        gfx.PanZoomController(self.camera, register_events=self.renderer)
        canvas.request_draw(lambda: self.renderer.render(self.scene, self.camera))
        loop.run()


def render_offscreen(points: np.ndarray, connectivity: np.ndarray, *,
                     size: Tuple[int, int] = (800, 600), **kwargs) -> np.ndarray:
    """Render a mesh once on the GPU; return an ``(h, w, 4)`` uint8 RGBA array."""
    return GPUMeshView(points, connectivity, size=size, **kwargs).render()
