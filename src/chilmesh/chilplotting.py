"""Standalone mesh plotting for the CHILmesh ecosystem.

This module decouples the rendering logic from the :class:`~chilmesh.CHILmesh`
class so it can be reused by *any* mesh generator in the stack (ADMESH,
QuADMesh, MADMESHing) — they emit ``points`` + ``connectivity`` arrays and can
plot directly, without constructing a CHILmesh or building its adjacency
structures::

    from chilmesh import chilplotting as cp
    fig, ax = cp.plot_mesh(points, triangles)          # edges only
    cp.plot_filled(points, triangles, values=quality)  # scalar colormap

Two layers are provided:

* **Array primitives** (``plot_mesh``, ``plot_filled``, ``plot_edges``,
  ``plot_points``, ``plot_quality_histogram``, ``build_polygons``,
  ``unique_edges``, ``configure_axes``) operate on bare numpy arrays. They are
  generator-agnostic and need no mesh topology beyond a connectivity table.
* **Mesh functions** (``plot``, ``plot_edge``, ``plot_boundary``,
  ``plot_layer``, ``plot_quality``, …) take a CHILmesh-like object and use its
  topology (edges, boundary, layers, quality). :class:`CHILmeshPlotMixin` in
  ``plot_utils`` is a thin set of wrappers over these.

All public functions return ``(Figure, Axes)`` (or ``Axes`` for
``configure_axes``) and accept an optional ``ax`` to compose subplots.
Rendering is vectorised — one ``LineCollection`` / ``PolyCollection`` per call —
so it scales to large meshes.

Render backends (CHILmesh #167): every function above draws with matplotlib,
which stays the default and the headless/CI fallback. :func:`render_image`
additionally offers an opt-in GPU path (``pygfx``/``wgpu``, ``pip install
chilmesh[gpu]``) implemented in :mod:`chilmesh.chilplotting_gpu` and imported
lazily. :func:`plot_backend_info` mirrors ``chilmesh.backend_info()`` (the
compute-backend selector in ``chilmesh/__init__.py``): same ``available`` /
``selected`` / ``versions`` keys and an environment override, here
``CHILMESH_PLOT_BACKEND={mpl,gpu}``.
"""
from __future__ import annotations

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import BoundaryNorm, Normalize
import numbers
import os
import warnings
from importlib.util import find_spec
from typing import Optional, Sequence, Tuple

__all__ = [
    # render backend selection (GPU backend is opt-in; matplotlib is the default)
    "plot_backend_info",
    "render_image",
    # array primitives (generator-agnostic)
    "configure_axes",
    "build_polygons",
    "unique_edges",
    "plot_mesh",
    "plot_filled",
    "plot_edges",
    "plot_points",
    "plot_quality_histogram",
    # mesh-object functions (CHILmesh-like)
    "axis_chilmesh",
    "plot",
    "plot_edge",
    "plot_boundary",
    "plot_interior_edges",
    "plot_elem",
    "plot_face",
    "plot_point",
    "plot_label",
    "plot_layer",
    "plot_quality",
    "plot_quality_from_mesh_histogram",
    "cartogram",
]


# ---------------------------------------------------------------------------
# Array primitives — operate on raw points + connectivity. No mesh required.
# ---------------------------------------------------------------------------

def plot_backend_info() -> dict:
    """Return the available and selected plot render backends.

    Mirrors ``chilmesh.backend_info()``. ``matplotlib`` is always ``"mpl"``;
    ``"gpu"`` is listed when ``pygfx`` and ``wgpu`` are installed (checked
    without importing them, so this call is cheap and never touches a GPU).
    Whether a GPU adapter exists is only known at render time;
    :func:`render_image` falls back to matplotlib with a warning if it does not.

    Returns
    -------
    dict
        ``available`` (list of names), ``selected`` (name used when
        ``render_image`` gets no ``backend=``), ``versions`` (name -> version).
        ``selected`` is ``"mpl"`` unless ``CHILMESH_PLOT_BACKEND=gpu`` is set
        and the GPU backend is available.
    """
    from importlib.metadata import PackageNotFoundError, version

    available = ["mpl"]
    versions = {"mpl": matplotlib.__version__}
    if find_spec("pygfx") is not None and find_spec("wgpu") is not None:
        available.append("gpu")
        try:
            versions["gpu"] = f"pygfx {version('pygfx')}, wgpu {version('wgpu')}"
        except PackageNotFoundError:  # pragma: no cover - spec found, metadata missing
            versions["gpu"] = "unknown"
    env_override = os.environ.get("CHILMESH_PLOT_BACKEND", "").strip().lower()
    selected = env_override if env_override in available else "mpl"
    return {"available": available, "selected": selected, "versions": versions}


def render_image(points: np.ndarray, connectivity: np.ndarray, *,
                 values: Optional[np.ndarray] = None, cmap: str = "viridis",
                 vmin: Optional[float] = None, vmax: Optional[float] = None,
                 edge_color: Optional[str] = "k", linewidth: float = 1.0,
                 size: Tuple[int, int] = (800, 600),
                 backend: Optional[str] = None) -> np.ndarray:
    """Render a mesh to an ``(height, width, 4)`` uint8 RGBA array.

    Parameters
    ----------
    points, connectivity : ndarray
        Bare mesh arrays (triangles, quads, or padded mixed connectivity).
    values : ndarray, optional
        One scalar per element, drawn through ``cmap``.
    backend : {"mpl", "gpu"}, optional
        Defaults to ``plot_backend_info()["selected"]``. ``"gpu"`` falls back
        to matplotlib with a ``RuntimeWarning`` when no GPU adapter is found.
    """
    name = (backend or plot_backend_info()["selected"]).lower()
    if name not in ("mpl", "gpu"):
        raise ValueError(f"backend must be 'mpl' or 'gpu', got {backend!r}")
    if name == "gpu":
        try:
            from . import chilplotting_gpu as _gpu
            if _gpu.gpu_available():
                return _gpu.render_offscreen(
                    points, connectivity, values=values, cmap=cmap, vmin=vmin,
                    vmax=vmax, edge_color=edge_color, linewidth=linewidth,
                    size=size)
            reason = "no GPU adapter found"
        except ImportError as exc:
            reason = str(exc)
        warnings.warn(f"GPU plot backend unavailable ({reason}); using matplotlib",
                      RuntimeWarning, stacklevel=2)

    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    fig = Figure(figsize=(size[0] / 100, size[1] / 100), dpi=100)
    FigureCanvasAgg(fig)
    ax = fig.add_axes([0, 0, 1, 1])
    plot_filled(points, connectivity, values=values, cmap=cmap, vmin=vmin,
                vmax=vmax, edgecolor=edge_color if edge_color else "none",
                linewidth=linewidth, ax=ax)
    ax.set_axis_off()
    fig.canvas.draw()
    return np.asarray(fig.canvas.buffer_rgba()).copy()


def configure_axes(points: np.ndarray, ax: Optional[plt.Axes] = None,
                   pad_frac: float = 0.01) -> plt.Axes:
    """Equal-aspect axes fitted to the point cloud (white background)."""
    if ax is None:
        ax = plt.gca()
    pts = np.asarray(points)
    ax.set_aspect("equal")
    ax.set_facecolor("white")
    ax.tick_params(labelsize=12, width=1.5)
    x, y = pts[:, 0], pts[:, 1]
    pad = pad_frac * max(x.max() - x.min(), y.max() - y.min())
    ax.set_xlim([x.min() - pad, x.max() + pad])
    ax.set_ylim([y.min() - pad, y.max() + pad])
    return ax


def _is_padded_or_degenerate(rows: np.ndarray) -> np.ndarray:
    """Boolean mask: 4-col rows that are really triangles (a repeated vertex)."""
    return (
        (rows[:, 0] == rows[:, 1])
        | (rows[:, 1] == rows[:, 2])
        | (rows[:, 2] == rows[:, 3])
        | (rows[:, 3] == rows[:, 0])
        | (rows[:, 0] == rows[:, 2])
        | (rows[:, 1] == rows[:, 3])
    )


def build_polygons(points: np.ndarray, connectivity: np.ndarray,
                   elem_ids: Optional[Sequence[int]] = None) -> list:
    """Return a list of ``(k, 2)`` polygon vertex arrays for the elements.

    Handles all-triangle, all-quad, and mixed / padded-triangle connectivity.
    Gathers are vectorised; only mixed meshes fall back to a Python list build
    (still O(1) matplotlib draw calls downstream). Polygons are returned in
    input row order, so ``polys[i]`` is element ``i`` (or ``elem_ids[i]``).
    """
    points = np.asarray(points)
    connectivity = np.asarray(connectivity)
    if elem_ids is None:
        rows = connectivity
    else:
        elem_ids = np.atleast_1d(np.asarray(elem_ids))
        rows = connectivity[elem_ids]
    if len(rows) == 0:
        return []

    n_cols = rows.shape[1]
    if n_cols == 3:
        return list(points[rows, :2])

    tri_mask = _is_padded_or_degenerate(rows)
    quad_mask = ~tri_mask
    if not quad_mask.any():
        return list(points[rows[:, :3], :2])
    if not tri_mask.any():
        return list(points[rows, :2])

    # Mixed: scatter triangles and quads back to their input row positions so
    # polygon i always belongs to element i (values align in plot_filled).
    polys = [None] * len(rows)
    for i, poly in zip(np.flatnonzero(tri_mask), points[rows[tri_mask, :3], :2]):
        polys[i] = poly
    for i, poly in zip(np.flatnonzero(quad_mask), points[rows[quad_mask], :2]):
        polys[i] = poly
    return polys


def unique_edges(connectivity: np.ndarray) -> np.ndarray:
    """Derive the unique undirected edges of a mesh from its connectivity.

    Adjacency-free: lets ``plot_mesh`` draw any generator's output without
    building CHILmesh edge structures. Returns an ``(m, 2)`` int array of
    vertex-index pairs (``v0 < v1``).
    """
    connectivity = np.asarray(connectivity)
    pairs = []
    for row in connectivity:
        verts = row.tolist()
        # drop a repeated padding vertex so a padded triangle yields 3 edges
        seen = []
        for v in verts:
            if v not in seen:
                seen.append(int(v))
        k = len(seen)
        for i in range(k):
            a, b = seen[i], seen[(i + 1) % k]
            pairs.append((a, b) if a < b else (b, a))
    if not pairs:
        return np.empty((0, 2), dtype=int)
    return np.unique(np.asarray(pairs, dtype=int), axis=0)


def _new_ax(ax: Optional[plt.Axes], figsize=(8, 6)) -> Tuple[plt.Figure, plt.Axes]:
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    return fig, ax


def plot_edges(points: np.ndarray, edge_verts: np.ndarray, *, color="g",
               linewidth: float = 2.5, linestyle: str = "-",
               ax: Optional[plt.Axes] = None) -> Tuple[plt.Figure, plt.Axes]:
    """Draw edges given an ``(m, 2)`` array of endpoint vertex indices."""
    fig, ax = _new_ax(ax)
    points = np.asarray(points)
    edge_verts = np.asarray(edge_verts)
    if len(edge_verts) == 0:
        return fig, ax
    p1 = points[edge_verts[:, 0], :2]
    p2 = points[edge_verts[:, 1], :2]
    segments = np.stack([p1, p2], axis=1)  # (m, 2, 2)
    lc = LineCollection(segments, colors=color, linewidths=linewidth,
                        linestyles=linestyle)
    ax.add_collection(lc)
    ax.autoscale_view()
    return fig, ax


def plot_mesh(points: np.ndarray, connectivity: np.ndarray, *,
              elem_ids: Optional[Sequence[int]] = None, edge_color="k",
              linewidth: float = 1.0, linestyle: str = "-",
              ax: Optional[plt.Axes] = None) -> Tuple[plt.Figure, plt.Axes]:
    """Draw a mesh wireframe (edges only) from raw arrays — no adjacency needed."""
    points = np.asarray(points)
    connectivity = np.asarray(connectivity)
    if elem_ids is not None:
        connectivity = connectivity[np.atleast_1d(np.asarray(elem_ids))]
    fig, ax = _new_ax(ax)
    configure_axes(points, ax=ax)
    plot_edges(points, unique_edges(connectivity), color=edge_color,
               linewidth=linewidth, linestyle=linestyle, ax=ax)
    return fig, ax


def plot_filled(points: np.ndarray, connectivity: np.ndarray, *,
                elem_ids: Optional[Sequence[int]] = None, facecolor="C0",
                values: Optional[np.ndarray] = None, cmap: str = "viridis",
                vmin: Optional[float] = None, vmax: Optional[float] = None,
                edgecolor="k", linewidth: float = 1.0, linestyle: str = "-",
                alpha: float = 1.0, colorbar: bool = False,
                clabel: Optional[str] = None,
                ax: Optional[plt.Axes] = None) -> Tuple[plt.Figure, plt.Axes]:
    """Fill mesh elements as polygons.

    With ``values`` (one scalar per drawn element) the polygons are coloured by
    a colormap — a generator-agnostic scalar field plot (quality, bathymetry,
    size field, …). Without ``values`` every polygon uses ``facecolor``.
    """
    fig, ax = _new_ax(ax)
    configure_axes(points, ax=ax)
    polys = build_polygons(points, connectivity, elem_ids)
    if not polys:
        return fig, ax

    if values is not None:
        values = np.asarray(values, dtype=float).ravel()
        if len(values) != len(polys):
            raise ValueError(
                f"values length {len(values)} != number of elements {len(polys)}")
        norm = Normalize(vmin=np.min(values) if vmin is None else vmin,
                         vmax=np.max(values) if vmax is None else vmax)
        cmap_obj = matplotlib.colormaps[cmap]
        pc = PolyCollection(polys, array=values, cmap=cmap_obj, norm=norm,
                            edgecolors=edgecolor, linewidths=linewidth,
                            linestyles=linestyle, alpha=alpha)
        ax.add_collection(pc)
        ax.autoscale_view()
        if colorbar:
            plt.colorbar(pc, ax=ax, label=clabel or "value")
    else:
        pc = PolyCollection(polys, facecolors=facecolor, edgecolors=edgecolor,
                            linewidths=linewidth, linestyles=linestyle,
                            alpha=alpha)
        ax.add_collection(pc)
        ax.autoscale_view()
    return fig, ax


def plot_points(points: np.ndarray, *, coords: Optional[np.ndarray] = None,
                ids: Optional[Sequence[int]] = None, marker: str = "o",
                color: str = "r", size: float = 5,
                ax: Optional[plt.Axes] = None, **kwargs
                ) -> Tuple[plt.Figure, plt.Axes]:
    """Scatter points. Pass ``coords`` (n,2) directly, or ``ids`` into ``points``."""
    fig, ax = _new_ax(ax)
    if coords is None:
        pts = np.asarray(points)
        coords = pts[np.arange(len(pts)) if ids is None
                     else np.atleast_1d(np.asarray(ids)), :2]
    coords = np.asarray(coords)
    ax.plot(coords[:, 0], coords[:, 1], linestyle="none", marker=marker,
            color=color, markersize=size, **kwargs)
    return fig, ax


def plot_quality_histogram(quality: np.ndarray, *, bins: int = 100,
                           cmap: str = "cool", auto_norm: bool = False,
                           ax: Optional[plt.Axes] = None
                           ) -> Tuple[plt.Figure, plt.Axes]:
    """Histogram of a per-element scalar (quality), bars coloured by ``cmap``.

    The reversed colormap (``cmap + "_r"``) matches :func:`plot_quality` so the
    1-D histogram reads as a companion to the 2-D map.
    """
    fig, ax = _new_ax(ax, figsize=(10, 4))
    q = np.asarray(quality, dtype=float).ravel()
    counts, edges = np.histogram(q, bins=bins, range=(0.0, 1.0))
    widths = np.diff(edges)
    midpoints = edges[:-1] + widths / 2.0
    if auto_norm and q.size > 0 and float(q.max() - q.min()) > 0:
        norm = Normalize(vmin=float(q.min()), vmax=float(q.max()))
    else:
        norm = Normalize(vmin=0.0, vmax=1.0)
    colors = matplotlib.colormaps[cmap + "_r"](norm(midpoints))
    ax.bar(edges[:-1], counts, width=widths, align="edge", color=colors,
           edgecolor="k", linewidth=0.3)
    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel("Quality")
    ax.set_ylabel("# of Elements")
    ax.set_title("Element Quality Distribution")
    return fig, ax


# ---------------------------------------------------------------------------
# Mesh-object functions — take a CHILmesh-like object and use its topology.
# CHILmeshPlotMixin delegates to these; they reuse the array primitives above.
# ---------------------------------------------------------------------------

def axis_chilmesh(mesh, ax: Optional[plt.Axes] = None) -> plt.Axes:
    """Configure axes with equal aspect and mesh-fitted limits."""
    return configure_axes(mesh.points, ax=ax)


def plot(mesh, elem_ids=None, elem_color="none", edge_color="k",
         linewidth=1.0, linestyle="-", ax=None) -> Tuple[plt.Figure, plt.Axes]:
    """Plot the mesh. ``elem_color='none'`` → fast edge-only wireframe."""
    if elem_ids is None:
        elem_ids = np.arange(mesh.n_elems)
    elif np.isscalar(elem_ids):
        elem_ids = np.array([elem_ids])

    fig, ax = _new_ax(ax)
    axis_chilmesh(mesh, ax=ax)

    if elem_color == "none":
        edges = mesh.elem2edge(elem_ids).flatten()
        edges = edges[edges >= 0]
        plot_edge(mesh, edges, color=edge_color, linewidth=linewidth,
                  linestyle=linestyle, ax=ax)
    else:
        _fill_mesh(mesh, elem_ids, facecolor=elem_color, edgecolor=edge_color,
                   linewidth=linewidth, linestyle=linestyle, ax=ax)
    return fig, ax


def plot_edge(mesh, edge_ids=None, color="g", linewidth=2.5, linestyle="-",
              ax=None) -> Tuple[plt.Figure, plt.Axes]:
    """Plot mesh edges (by edge ID) using a single vectorized LineCollection."""
    if edge_ids is None:
        edge_ids = np.arange(mesh.n_edges)
    if ax is None:
        ax = axis_chilmesh(mesh)
    return plot_edges(mesh.points, mesh.edge2vert(edge_ids), color=color,
                      linewidth=linewidth, linestyle=linestyle, ax=ax)


def plot_boundary(mesh, color="r", linewidth=2.5, linestyle="-",
                  ax=None) -> Tuple[plt.Figure, plt.Axes]:
    """Plot mesh exterior (boundary) edges highlighted."""
    return plot_edge(mesh, mesh.boundary_edges(), color=color,
                     linewidth=linewidth, linestyle=linestyle, ax=ax)


def plot_interior_edges(mesh, color="b", linewidth=1.0, linestyle="-",
                        ax=None) -> Tuple[plt.Figure, plt.Axes]:
    """Plot interior (shared by two elements) edges."""
    interior = np.setdiff1d(np.arange(mesh.n_edges), mesh.boundary_edges())
    return plot_edge(mesh, interior, color=color, linewidth=linewidth,
                     linestyle=linestyle, ax=ax)


def plot_elem(mesh, elem_ids=None, color="b", edge_color="k", linewidth=1.0,
              linestyle="-", ax=None) -> Tuple[plt.Figure, plt.Axes]:
    """Plot elements filled with a single color using PolyCollection."""
    if np.isscalar(elem_ids):
        elem_ids = np.array([elem_ids])
    if elem_ids is None:
        elem_ids = np.arange(mesh.n_elems)
    if ax is None:
        fig, ax = plot(mesh)
    else:
        fig = ax.figure
    _fill_mesh(mesh, elem_ids, facecolor=color, edgecolor=edge_color,
               linewidth=linewidth, linestyle=linestyle, ax=ax)
    return fig, ax


def plot_face(mesh, face_ids=None, color="b", edge_color="k", linewidth=1.0,
              linestyle="-", ax=None) -> Tuple[plt.Figure, plt.Axes]:
    """Plot faces (alias for :func:`plot_elem`)."""
    return plot_elem(mesh, elem_ids=face_ids, color=color, edge_color=edge_color,
                     linewidth=linewidth, linestyle=linestyle, ax=ax)


def plot_point(mesh, ids=None, point_type: str = "vertex", color: str = "r",
               marker: str = "o", size: float = 5, ax=None, **kwargs
               ) -> Tuple[plt.Figure, plt.Axes]:
    """Plot mesh entities as scatter points (vertices, edge midpoints, centroids)."""
    if ax is None:
        fig, ax = plot(mesh)
    else:
        fig = ax.figure

    pt = point_type.lower()
    if pt in {"vertex", "vert"}:
        if ids is None:
            ids = np.arange(mesh.n_verts)
        x, y = mesh.points[ids, 0], mesh.points[ids, 1]
    elif pt in {"edge", "midpoint"}:
        if ids is None:
            ids = np.arange(mesh.n_edges)
        edges = mesh.edge2vert(ids)
        x = np.mean(mesh.points[edges, 0], axis=1)
        y = np.mean(mesh.points[edges, 1], axis=1)
    elif pt in {"element", "centroid"}:
        if ids is None:
            ids = np.arange(mesh.n_elems)
        ids = np.atleast_1d(np.asarray(ids))
        conn = mesh.connectivity_list[ids]
        coords = mesh.points[conn, :2]
        x = coords[:, :, 0].mean(axis=1)
        y = coords[:, :, 1].mean(axis=1)
    else:
        raise ValueError(f"Unknown point_type: {point_type!r}")

    ax.plot(x, y, linestyle="none", marker=marker, color=color,
            markersize=size, **kwargs)
    return fig, ax


def plot_label(mesh, ids=None, label: str = "all", ax=None
               ) -> Tuple[plt.Figure, plt.Axes]:
    """Annotate mesh entities with their IDs."""
    if ax is None:
        fig, ax = plot(mesh)
    else:
        fig = ax.figure

    if label in {"vertex", "point", "all"}:
        vert_ids = ids if ids is not None else np.arange(mesh.n_verts)
        for i in vert_ids:
            ax.text(mesh.points[i, 0], mesh.points[i, 1], f"V{i}",
                    color="red", ha="center", fontsize=7)
    if label in {"edge", "all"}:
        edge_ids = ids if ids is not None else np.arange(mesh.n_edges)
        for i in edge_ids:
            ev = mesh.edge2vert([i])
            x = np.mean(mesh.points[ev, 0])
            y = np.mean(mesh.points[ev, 1])
            ax.text(x, y, f"E{i}", color="green", ha="center", fontsize=7)
    if label in {"element", "centroid", "all"}:
        elem_ids = ids if ids is not None else np.arange(mesh.n_elems)
        for i in elem_ids:
            verts = mesh.connectivity_list[i]
            verts = verts[verts >= 0]
            x = np.mean(mesh.points[verts, 0])
            y = np.mean(mesh.points[verts, 1])
            ax.text(x, y, f"El{i}", color="blue", ha="center", fontsize=7)
    return fig, ax


def plot_layer(mesh, layers=None, cmap="viridis", ax=None
               ) -> Tuple[plt.Figure, plt.Axes]:
    """Plot peel layers as filled elements colored by layer index."""
    if layers is None:
        layers = range(mesh.n_layers)
    fig, ax = _new_ax(ax, figsize=(10, 8))
    axis_chilmesh(mesh, ax=ax)

    cmap_obj = matplotlib.colormaps[cmap].resampled(mesh.n_layers)
    norm = BoundaryNorm(boundaries=np.arange(mesh.n_layers + 1),
                        ncolors=mesh.n_layers)
    for layer_idx in layers:
        if (layer_idx >= len(mesh.layers["OE"])
                or layer_idx >= len(mesh.layers["IE"])):
            continue
        elem_ids = np.concatenate(
            (mesh.layers["OE"][layer_idx], mesh.layers["IE"][layer_idx])
        ).astype(int)
        color = cmap_obj(norm(layer_idx))
        _fill_mesh(mesh, elem_ids, facecolor=color, edgecolor="k",
                   linewidth=0.5, alpha=0.7, ax=ax)

    sm = cm.ScalarMappable(norm=norm, cmap=cmap_obj)
    sm.set_array([])
    plt.colorbar(sm, ax=ax, label="Layer",
                 ticks=np.arange(0, mesh.n_layers) + 0.5, format="%d")
    ax.set_title(f"Mesh Layers ({mesh.n_layers} total)")
    return fig, ax


def plot_quality(mesh, elem_ids=None, cmap="cool", ax=None
                 ) -> Tuple[plt.Figure, plt.Axes]:
    """Plot element quality as a binned filled colormap."""
    if elem_ids is None:
        elem_ids = np.arange(mesh.n_elems)
    elif np.isscalar(elem_ids):
        elem_ids = np.array([elem_ids])
    fig, ax = _new_ax(ax, figsize=(10, 8))
    axis_chilmesh(mesh, ax=ax)

    q, _, _ = mesh.elem_quality(elem_ids=elem_ids)
    bins = np.linspace(0, 1, 21)
    norm = Normalize(vmin=0, vmax=1)
    colors = matplotlib.colormaps[cmap + "_r"](norm(bins[:-1]))
    for i, bin_lo in enumerate(bins[:-1]):
        bin_hi = bins[i + 1]
        if i == len(bins) - 2:
            mask = (q >= bin_lo) & (q <= bin_hi)
        else:
            mask = (q >= bin_lo) & (q < bin_hi)
        if not np.any(mask):
            continue
        _fill_mesh(mesh, elem_ids[mask], facecolor=colors[i], edgecolor="k",
                   linewidth=0.5, ax=ax)

    sm = cm.ScalarMappable(norm=norm, cmap=matplotlib.colormaps[cmap + "_r"])
    sm.set_array([])
    plt.colorbar(sm, ax=ax, label="Element Quality")
    ax.set_title("Element Quality")
    return fig, ax


def plot_quality_from_mesh_histogram(mesh, elem_ids=None, bins: int = 100,
                                     cmap: str = "cool", auto_norm: bool = False,
                                     ax=None) -> Tuple[plt.Figure, plt.Axes]:
    """Quality histogram for a mesh — computes quality then delegates to the
    array primitive :func:`plot_quality_histogram`."""
    if elem_ids is None:
        elem_ids = np.arange(mesh.n_elems)
    elif np.isscalar(elem_ids):
        elem_ids = np.array([elem_ids])
    q, _, _ = mesh.elem_quality(elem_ids=elem_ids)
    return plot_quality_histogram(q, bins=bins, cmap=cmap, auto_norm=auto_norm,
                                  ax=ax)


def _fill_mesh(mesh, elem_ids, facecolor, edgecolor="k", linewidth=1.0,
               linestyle="-", alpha=1.0, ax=None) -> None:
    """Fill specific mesh elements (by ID) on ``ax`` via a PolyCollection."""
    elem_ids = np.atleast_1d(np.asarray(elem_ids))
    if len(elem_ids) == 0:
        return
    polys = build_polygons(mesh.points, mesh.connectivity_list, elem_ids)
    if not polys:
        return
    pc = PolyCollection(polys, facecolors=facecolor, edgecolors=edgecolor,
                        linewidths=linewidth, linestyles=linestyle, alpha=alpha)
    ax.add_collection(pc)
    ax.autoscale_view()


def plot_field(mesh, values, *, ax: Optional[plt.Axes] = None, cmap: str = "viridis",
               colorbar: bool = True, clabel: Optional[str] = None,
               vmin: Optional[float] = None, vmax: Optional[float] = None,
               edgecolor="none") -> Tuple[plt.Figure, plt.Axes]:
    """Draw a node-aligned field (e.g. one ADCIRC output record) on the mesh.

    Each element takes the mean of its node values. An element with any NaN
    node (dry or fill) is left blank. ``values`` must have ``mesh.n_verts``
    entries or ``ValueError`` names both lengths.

    Example:
        >>> import matplotlib; matplotlib.use("Agg")
        >>> from chilmesh import examples
        >>> m = examples.annulus()
        >>> fig, ax = plot_field(m, m.points[:, 0])
        >>> len(fig.axes)  # main axes and the color bar
        2

    Tested by: tests/test_plot_field.py::test_colorbar_present_and_one_polygon_per_element.
    """
    vals = np.asarray(values, dtype=float).ravel()
    if len(vals) != mesh.n_verts:
        raise ValueError(
            f"values length {len(vals)} != number of mesh nodes {mesh.n_verts}")
    conn = np.asarray(mesh.connectivity_list)
    pad = (conn.shape[1] == 4) & (conn[:, 3] == conn[:, 0]) if conn.shape[1] == 4 \
        else np.zeros(len(conn), dtype=bool)
    vsum = vals[conn].sum(axis=1)
    cnt = np.where(pad, 3, conn.shape[1])
    elem_vals = vsum / cnt
    if conn.shape[1] == 4:  # padded triangles repeat node 0; do not count it twice
        elem_vals = np.where(pad, (vsum - vals[conn[:, 3]]) / 3, elem_vals)
    finite = elem_vals[np.isfinite(elem_vals)]
    lo = (float(finite.min()) if finite.size else 0.0) if vmin is None else vmin
    hi = (float(finite.max()) if finite.size else 1.0) if vmax is None else vmax
    return plot_filled(np.asarray(mesh.points), conn, values=elem_vals, cmap=cmap,
                       vmin=lo, vmax=hi, edgecolor=edgecolor, colorbar=colorbar,
                       clabel=clabel, ax=ax)


# ---------------------------------------------------------------------------
# Size-controlled cartograms (CHILmesh #219)
# ---------------------------------------------------------------------------

_CARTOGRAM_KINDS = ("noncontig", "unrolled")
_CARTOGRAM_REFERENCES = ("mean", "median", "mode", "max", "min")
_BAND_FILL = 0.92  # fraction of a unit-height band covered by its cells


def _cartogram_elem_values(mesh, values, on) -> np.ndarray:
    """Return one float per element from per-element or per-vertex ``values``."""
    if on not in (None, "elem", "vertex"):
        raise ValueError(f"on must be 'elem', 'vertex' or None, got {on!r}")
    vals = np.asarray(values, dtype=float).ravel()
    conn = np.asarray(mesh.connectivity_list)
    n_elems, n_verts = len(conn), len(mesh.points)
    n = len(vals)
    if n_elems == n_verts and n == n_elems and on is None:
        raise ValueError(
            f"values length {n} matches both the number of elements ({n_elems}) "
            f"and the number of vertices ({n_verts}); pass on='elem' or "
            f"on='vertex'")
    if on is None:
        if n == n_elems:
            on = "elem"
        elif n == n_verts:
            on = "vertex"
        else:
            raise ValueError(
                f"values length {n} matches neither the number of elements "
                f"({n_elems}) nor the number of vertices ({n_verts})")
    expected = n_elems if on == "elem" else n_verts
    if n != expected:
        raise ValueError(
            f"values length {n} != number of {'elements' if on == 'elem' else 'vertices'}"
            f" ({expected}) for on={on!r} (elements: {n_elems}, vertices: {n_verts})")
    if on == "elem":
        return vals
    # per-vertex -> per-element mean, exactly as plot_field (padded triangles
    # repeat node 0 and must not count it twice; any NaN node blanks the element)
    pad = (conn.shape[1] == 4) & (conn[:, 3] == conn[:, 0]) if conn.shape[1] == 4 \
        else np.zeros(len(conn), dtype=bool)
    vsum = vals[conn].sum(axis=1)
    cnt = np.where(pad, 3, conn.shape[1])
    elem_vals = vsum / cnt
    if conn.shape[1] == 4:
        elem_vals = np.where(pad, (vsum - vals[conn[:, 3]]) / 3, elem_vals)
    return elem_vals


def _cartogram_reference(finite: np.ndarray, reference) -> Optional[float]:
    """Resolve ``reference`` against the finite values (None stays None)."""
    if reference is None:
        return None
    if isinstance(reference, str):
        if reference not in _CARTOGRAM_REFERENCES:
            raise ValueError(
                f"reference must be None, a number or one of "
                f"{_CARTOGRAM_REFERENCES}, got {reference!r}")
        if finite.size == 0:
            return None
        if reference == "mean":
            return float(np.mean(finite))
        if reference == "median":
            return float(np.median(finite))
        if reference == "max":
            return float(np.max(finite))
        if reference == "min":
            return float(np.min(finite))
        uniq, counts = np.unique(finite, return_counts=True)  # exact values
        return float(uniq[int(np.argmax(counts))])  # argmax -> first -> smallest
    if isinstance(reference, (bool, np.bool_)) or not isinstance(
            reference, numbers.Real):
        raise ValueError(
            f"reference must be None, a number or one of "
            f"{_CARTOGRAM_REFERENCES}, got {reference!r}")
    return float(reference)


def _cartogram_geometry(mesh):
    """Per-element ``(verts (n,4,2), is_tri, centroid (n,2), area (n,))``.

    Triangles are padded to four vertices by repeating vertex 0, so the shoelace
    area and every vectorised gather work for tri, quad and mixed meshes.
    """
    pts = np.asarray(mesh.points)[:, :2]
    conn = np.asarray(mesh.connectivity_list)
    if conn.shape[1] == 3:
        is_tri = np.ones(len(conn), dtype=bool)
        conn4 = np.concatenate([conn, conn[:, :1]], axis=1)
    else:
        is_tri = _is_padded_or_degenerate(conn[:, :4])
        conn4 = conn[:, :4].copy()
        conn4[is_tri, 3] = conn4[is_tri, 0]
    verts = pts[conn4]
    nv = np.where(is_tri, 3, 4)
    cent = (verts[:, :3].sum(axis=1) + np.where(is_tri, 0, 1)[:, None] * verts[:, 3]) \
        / nv[:, None]
    x, y = verts[..., 0], verts[..., 1]
    area = 0.5 * np.abs((x * np.roll(y, -1, axis=1) - np.roll(x, -1, axis=1) * y).sum(axis=1))
    return verts, is_tri, cent, area


def _cartogram_noncontig(mesh, elem_vals, finite_mask, ref, eps, scale_clip,
                         cmap, vmin, vmax, ax):
    verts, is_tri, cent, area = _cartogram_geometry(mesh)
    target = float(np.median(area))  # A: median element area
    idx = np.flatnonzero(finite_mask)
    a = area[idx]
    if target > 0:
        base = np.where(a > 0, np.sqrt(target / np.where(a > 0, a, 1.0)), 1.0)
    else:
        base = np.ones(len(idx))
    scale = base
    if ref is not None and len(idx):
        w = np.abs(elem_vals[idx] - ref) + eps
        denom = float(np.median(w))
        if denom <= 0:
            denom = float(np.mean(w))
        wn = w / denom if denom > 0 else np.ones(len(idx))
        # clip the focus factor only, not the final scale
        scale = base * np.clip(np.sqrt(wn), scale_clip[0], scale_clip[1])

    order = np.argsort(elem_vals[idx], kind="stable")  # high values drawn last
    sel = idx[order]
    c = cent[sel][:, None, :]
    scaled = c + (verts[sel] - c) * scale[order][:, None, None]
    tri_sel = is_tri[sel]
    if tri_sel.all():
        polys = scaled[:, :3]
    elif tri_sel.any():
        polys = [s[:3] if t else s for s, t in zip(scaled, tri_sel)]
    else:
        polys = scaled
    pc = PolyCollection(polys, array=elem_vals[sel], cmap=matplotlib.colormaps[cmap],
                        norm=Normalize(vmin=vmin, vmax=vmax), edgecolors="none")
    extent = scaled.reshape(-1, 2) if len(sel) else np.asarray(mesh.points)[:, :2]
    configure_axes(extent, ax=ax)
    ax.add_collection(pc)
    return pc


def _cartogram_unrolled(mesh, elem_vals, finite_mask, ref, eps, cmap, vmin, vmax,
                        ax):
    _, _, cent, _ = _cartogram_geometry(mesh)
    mid = cent.mean(axis=0)
    angle = np.arctan2(cent[:, 1] - mid[1], cent[:, 0] - mid[0])
    oe, ie = mesh.layers["OE"], mesh.layers["IE"]
    n_bands = min(len(oe), len(ie))
    rects, cell_vals, ticks, labels = [], [], [], []
    for i in range(n_bands):
        members = np.unique(np.concatenate((oe[i], ie[i])).astype(int))
        members = members[np.argsort(angle[members], kind="stable")]
        n = len(members)
        if n == 0:
            continue
        if ref is None:
            w = np.ones(n)
        else:
            w = np.where(finite_mask[members],
                         np.abs(elem_vals[members] - ref) + eps, eps)
        total = float(w.sum())
        w = w / total if total > 0 else np.full(n, 1.0 / n)
        x = np.concatenate(([0.0], np.cumsum(w)))
        x[-1] = 1.0
        y0 = float(n_bands - 1 - i)  # layer 0 at the top
        y1 = y0 + _BAND_FILL
        rects.append(np.stack([
            np.stack([x[:-1], np.full(n, y0)], axis=1),
            np.stack([x[1:], np.full(n, y0)], axis=1),
            np.stack([x[1:], np.full(n, y1)], axis=1),
            np.stack([x[:-1], np.full(n, y1)], axis=1)], axis=1))
        cell_vals.append(elem_vals[members])  # NaN stays: slot kept, drawn blank
        ticks.append(y0 + _BAND_FILL / 2)
        labels.append(f"L{i} (n={n})")
    polys = np.concatenate(rects) if rects else np.empty((0, 4, 2))
    vals = np.concatenate(cell_vals) if cell_vals else np.empty(0)
    pc = PolyCollection(polys, array=vals, cmap=matplotlib.colormaps[cmap],
                        norm=Normalize(vmin=vmin, vmax=vmax), edgecolors="none")
    ax.set_facecolor("white")
    ax.add_collection(pc)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, max(n_bands, 1))
    step = max(1, len(ticks) // 25)  # thin labels on very deep meshes
    ax.set_yticks(ticks[::step])
    ax.set_yticklabels(labels[::step], fontsize=8)
    if ref is None:
        ax.set_xlabel("angular position about the mesh centre (-pi to pi), "
                      "equal-width cells")
    else:
        ax.set_xlabel(f"cell width ~ |value - {ref:.4g}| + eps "
                      f"(cells near the reference are compressed)")
    return pc


def cartogram(mesh, values, *, kind: str = "noncontig", reference=None,
              on: Optional[str] = None, eps: Optional[float] = None,
              scale_clip: Tuple[float, float] = (0.05, 10), cmap: str = "inferno",
              colorbar: bool = True, clabel: Optional[str] = None,
              ax: Optional[plt.Axes] = None) -> Tuple[plt.Figure, plt.Axes]:
    """Draw a size-controlled cartogram of a scalar field on a mesh.

    A plain filled mesh lets large elements (coarse regions) dominate the
    canvas while small refined elements all but vanish. A **cartogram** redraws
    the mesh so each element's visual size is set by the plot, not by its
    geometric size. Four terms are used below.

    * **layout** (``kind``): ``"noncontig"`` keeps every element in its true
      position and redraws it in its own shape scaled about its centroid to the
      median element area (a non-contiguous cartogram, so elements may overlap
      or leave gaps). ``"unrolled"`` discards geography: the mesh is unrolled
      into horizontal bands and each element becomes a rectangular cell.
    * **band**: one skeleton layer ``i``, that is the elements in
      ``mesh.layers['OE'][i]`` and ``mesh.layers['IE'][i]``. Bands have equal
      height, layer 0 at the top; inside a band cells are ordered by the angle
      ``atan2`` of the element centroid about the mean of all element
      centroids. Every element sits in exactly one band.
    * **reference**: a value the field is compared with. ``None`` gives equal
      size for every element (noncontig: median area; unrolled: equal cell
      width). Otherwise the plot is focus+context: element weight is
      ``|value - reference| + eps``, so elements near the reference shrink and
      outliers grow. In noncontig the linear size is multiplied by
      ``clip(sqrt(weight / median(weight)), *scale_clip)``; the clip limits
      only that focus factor, not the final size. In unrolled the cell width
      is proportional to the weight within its band.
    * **colour**: ``cmap`` maps the value; noncontig draws ascending by value,
      so high values lie on top.

    Args:
        mesh: CHILmesh-like object (``points``, ``connectivity_list``,
            ``layers``). Tri, quad and mixed (padded-triangle) meshes work.
        values: one scalar per element, or one per vertex (reduced to the
            element mean exactly as :func:`plot_field`, NaN node blanks the
            element). Which one is inferred from the length.
        kind: ``"noncontig"`` (default) or ``"unrolled"``.
        reference: ``None`` (equal size), ``"mean"``, ``"median"``, ``"mode"``,
            ``"max"``, ``"min"`` or a number. ``"mode"`` is the exact most
            frequent value (``np.unique`` counts, ties go to the smallest
            value) and is intended for discrete fields.
        on: ``"elem"`` or ``"vertex"``. Required only when the mesh has as many
            elements as vertices and ``values`` has that length.
        eps: absolute floor added to every weight. Default
            ``0.004 * (max - min)`` of the finite values.
        scale_clip: ``(low, high)`` limits on the noncontig focus factor.
        cmap: Matplotlib colormap name.
        colorbar: draw a colour bar labelled ``clabel``.
        clabel: colour-bar label (default ``"value"``).
        ax: axes to draw on; a new figure when ``None``.

    Returns:
        ``(Figure, Axes)``.

    Raises:
        ValueError: unknown ``kind`` or ``reference``; ``on`` not ``"elem"`` /
            ``"vertex"``; ``values`` length matches neither the element nor the
            vertex count (or the one ``on`` selects); ``on`` missing when the
            two counts are equal; negative ``eps``; bad ``scale_clip``.

    NaN (or inf) elements are blank (noncontig: not drawn; unrolled: an empty
    slot, of width ``eps`` under a reference) and are excluded from the
    reference, the colour limits and the default ``eps``. A constant field
    falls back to equal size.

    Example:
        >>> import matplotlib; matplotlib.use("Agg")
        >>> from chilmesh import examples
        >>> m = examples.annulus()
        >>> q, _, _ = m.elem_quality()
        >>> fig, ax = cartogram(m, q, kind="unrolled", reference="mean")
        >>> len(fig.axes)  # main axes and the color bar
        2

    Tested by: tests/test_cartogram.py::test_returns_figure_axes.
    """
    if kind not in _CARTOGRAM_KINDS:
        raise ValueError(f"kind must be one of {_CARTOGRAM_KINDS}, got {kind!r}")
    if eps is not None and not eps >= 0:
        raise ValueError(f"eps must be >= 0, got {eps!r}")
    if len(scale_clip) != 2 or not 0 < scale_clip[0] <= scale_clip[1]:
        raise ValueError(
            f"scale_clip must be (low, high) with 0 < low <= high, got {scale_clip!r}")
    elem_vals = _cartogram_elem_values(mesh, values, on)
    finite_mask = np.isfinite(elem_vals)
    finite = elem_vals[finite_mask]
    ref = _cartogram_reference(finite, reference)
    if finite.size:
        vmin, vmax = float(finite.min()), float(finite.max())
    else:
        vmin, vmax = 0.0, 1.0
    if finite.size == 0 or vmax == vmin:
        ref = None  # constant (or empty) field: equal size
    eps_val = float(eps) if eps is not None else 0.004 * (vmax - vmin)

    fig, ax = _new_ax(ax, figsize=(10, 8))
    if kind == "noncontig":
        pc = _cartogram_noncontig(mesh, elem_vals, finite_mask, ref, eps_val,
                                  scale_clip, cmap, vmin, vmax, ax)
    else:
        pc = _cartogram_unrolled(mesh, elem_vals, finite_mask, ref, eps_val,
                                 cmap, vmin, vmax, ax)
    if colorbar:
        fig.colorbar(pc, ax=ax, label=clabel or "value")
    return fig, ax
