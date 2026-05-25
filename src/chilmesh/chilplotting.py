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

GPU-accelerated live animation of very large meshes (CHILmesh #75) is a
separate, future effort and intentionally out of scope here.
"""
from __future__ import annotations

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import BoundaryNorm, Normalize
from typing import Optional, Sequence, Tuple

__all__ = [
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
]


# ---------------------------------------------------------------------------
# Array primitives — operate on raw points + connectivity. No mesh required.
# ---------------------------------------------------------------------------

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
    (still O(1) matplotlib draw calls downstream).
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

    polys = []
    if tri_mask.any():
        polys.extend(points[rows[tri_mask, :3], :2])
    if quad_mask.any():
        polys.extend(points[rows[quad_mask], :2])
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
    """Plot skeletonization layers as filled elements colored by layer index."""
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
