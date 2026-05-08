# plot_utils.py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import BoundaryNorm, Normalize
from typing import Optional, Tuple


class CHILmeshPlotMixin:
    """Mixin providing fast, vectorized plotting for CHILmesh objects."""

    @property
    def grid_name(self) -> str:
        return getattr(self, '_grid_name', 'CHILmesh')

    def axis_chilmesh(self, ax: Optional[plt.Axes] = None) -> plt.Axes:
        """Configure axes with equal aspect and mesh-fitted limits."""
        if ax is None:
            ax = plt.gca()
        ax.set_aspect('equal')
        ax.set_facecolor('white')
        ax.tick_params(labelsize=12, width=1.5)
        x = self.points[:, 0]
        y = self.points[:, 1]
        pad = 0.01 * max(x.max() - x.min(), y.max() - y.min())
        ax.set_xlim([x.min() - pad, x.max() + pad])
        ax.set_ylim([y.min() - pad, y.max() + pad])
        return ax

    def plot(self, elem_ids=None, elem_color='none', edge_color='k',
             linewidth=1.0, linestyle='-', ax=None) -> Tuple[plt.Figure, plt.Axes]:
        """
        Plot the mesh.

        When elem_color is 'none', only edges are drawn (fast LineCollection path).
        Otherwise elements are filled using PolyCollection.
        """
        if elem_ids is None:
            elem_ids = np.arange(self.n_elems)
        elif np.isscalar(elem_ids):
            elem_ids = np.array([elem_ids])

        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 6))
        else:
            fig = ax.figure

        self.axis_chilmesh(ax=ax)

        if elem_color == 'none':
            edges = self.elem2edge(elem_ids).flatten()
            edges = edges[edges >= 0]
            self.plot_edge(edges, color=edge_color, linewidth=linewidth,
                           linestyle=linestyle, ax=ax)
        else:
            self._plot_polys(elem_ids, facecolor=elem_color,
                             edgecolor=edge_color, linewidth=linewidth,
                             linestyle=linestyle, ax=ax)

        return fig, ax

    def plot_edge(self, edge_ids=None, color='g', linewidth=2.5,
                  linestyle='-', ax=None) -> Tuple[plt.Figure, plt.Axes]:
        """Plot edges using a single vectorized LineCollection call."""
        if edge_ids is None:
            edge_ids = np.arange(self.n_edges)

        if ax is None:
            ax = self.axis_chilmesh()
        fig = ax.figure

        v = self.edge2vert(edge_ids)
        p1 = self.points[v[:, 0], :2]
        p2 = self.points[v[:, 1], :2]

        # One call instead of n_edges calls — critical for large meshes
        segments = np.stack([p1, p2], axis=1)  # (n_edges, 2, 2)
        lc = LineCollection(segments, colors=color, linewidths=linewidth,
                            linestyles=linestyle)
        ax.add_collection(lc)

        return fig, ax

    def plot_elem(self, elem_ids=None, color='b', edge_color='k',
                  linewidth=1.0, linestyle='-',
                  ax=None) -> Tuple[plt.Figure, plt.Axes]:
        """Plot elements filled with a single color using PolyCollection."""
        elem_ids = self._ensure_array(elem_ids)
        if elem_ids is None:
            elem_ids = np.arange(self.n_elems)

        if ax is None:
            fig, ax = self.plot()
        else:
            fig = ax.figure

        self._plot_polys(elem_ids, facecolor=color, edgecolor=edge_color,
                         linewidth=linewidth, linestyle=linestyle, ax=ax)
        return fig, ax

    def plot_face(self, face_ids=None, color='b', edge_color='k',
                  linewidth=1.0, linestyle='-',
                  ax=None) -> Tuple[plt.Figure, plt.Axes]:
        """Plot faces (alias for plot_elem)."""
        return self.plot_elem(elem_ids=face_ids, color=color,
                              edge_color=edge_color, linewidth=linewidth,
                              linestyle=linestyle, ax=ax)

    def plot_point(self, ids: Optional[np.ndarray] = None,
                   point_type: str = 'vertex', color: str = 'r',
                   marker: str = 'o', size: float = 5,
                   ax: plt.Axes = None, **kwargs) -> Tuple[plt.Figure, plt.Axes]:
        """Plot mesh entities as scatter points (vertices, edge midpoints, or centroids)."""
        if ax is None:
            fig, ax = self.plot()
        else:
            fig = ax.figure

        if point_type.lower() in {'vertex', 'vert'}:
            if ids is None:
                ids = np.arange(self.n_verts)
            x, y = self.points[ids, 0], self.points[ids, 1]

        elif point_type.lower() in {'edge', 'midpoint'}:
            if ids is None:
                ids = np.arange(self.n_edges)
            edges = self.edge2vert(ids)
            x = np.mean(self.points[edges, 0], axis=1)
            y = np.mean(self.points[edges, 1], axis=1)

        elif point_type.lower() in {'element', 'centroid'}:
            if ids is None:
                ids = np.arange(self.n_elems)
            centroids = np.zeros((len(ids), 2))
            for i, elem_id in enumerate(ids):
                verts = self.connectivity_list[elem_id]
                verts = verts[verts >= 0]
                centroids[i] = np.mean(self.points[verts, :2], axis=0)
            x, y = centroids[:, 0], centroids[:, 1]

        else:
            raise ValueError(f"Unknown point_type: {point_type!r}")

        ax.plot(x, y, linestyle='none', marker=marker, color=color,
                markersize=size, **kwargs)
        return fig, ax

    def plot_label(self, ids: Optional[np.ndarray] = None, label: str = 'all',
                   ax: Optional[plt.Axes] = None) -> Tuple[plt.Figure, plt.Axes]:
        """Annotate mesh entities with their IDs."""
        if ax is None:
            fig, ax = self.plot()
        else:
            fig = ax.figure

        if label in {'vertex', 'point', 'all'}:
            vert_ids = ids if ids is not None else np.arange(self.n_verts)
            for i in vert_ids:
                ax.text(self.points[i, 0], self.points[i, 1], f'V{i}',
                        color='red', ha='center', fontsize=7)

        if label in {'edge', 'all'}:
            edge_ids = ids if ids is not None else np.arange(self.n_edges)
            for i in edge_ids:
                ev = self.edge2vert([i])
                x = np.mean(self.points[ev, 0])
                y = np.mean(self.points[ev, 1])
                ax.text(x, y, f'E{i}', color='green', ha='center', fontsize=7)

        if label in {'element', 'centroid', 'all'}:
            elem_ids = ids if ids is not None else np.arange(self.n_elems)
            for i in elem_ids:
                verts = self.connectivity_list[i]
                verts = verts[verts >= 0]
                x = np.mean(self.points[verts, 0])
                y = np.mean(self.points[verts, 1])
                ax.text(x, y, f'El{i}', color='blue', ha='center', fontsize=7)

        return fig, ax

    def plot_layer(self, layers=None, cmap='viridis',
                   ax=None) -> Tuple[plt.Figure, plt.Axes]:
        """Plot skeletonization layers as filled elements colored by layer index."""
        if layers is None:
            layers = range(self.n_layers)

        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 8))
        else:
            fig = ax.figure
        ax.set_aspect('equal')

        cmap_obj = matplotlib.colormaps[cmap].resampled(self.n_layers)
        norm = BoundaryNorm(boundaries=np.arange(self.n_layers + 1),
                            ncolors=self.n_layers)

        for layer_idx in layers:
            if layer_idx >= len(self.layers["OE"]) or layer_idx >= len(self.layers["IE"]):
                continue
            elem_ids = np.concatenate(
                (self.layers["OE"][layer_idx], self.layers["IE"][layer_idx])
            ).astype(int)
            color = cmap_obj(norm(layer_idx))
            self._plot_polys(elem_ids, facecolor=color, edgecolor='k',
                             linewidth=0.5, alpha=0.7, ax=ax)

        sm = cm.ScalarMappable(norm=norm, cmap=cmap_obj)
        sm.set_array([])
        plt.colorbar(sm, ax=ax, label='Layer',
                     ticks=np.arange(0, self.n_layers) + 0.5,
                     format='%d')
        ax.set_title(f"Mesh Layers ({self.n_layers} total)")
        return fig, ax

    def plot_quality(self, elem_ids=None, cmap='cool',
                     ax=None) -> Tuple[plt.Figure, plt.Axes]:
        """
        Plot element quality as a filled contour map.

        Uses a single PolyCollection per colormap bin, so rendering scales
        with O(n_bins) matplotlib calls instead of O(n_elements).
        """
        if elem_ids is None:
            elem_ids = np.arange(self.n_elems)
        elif np.isscalar(elem_ids):
            elem_ids = np.array([elem_ids])

        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 8))
        else:
            fig = ax.figure
        self.axis_chilmesh(ax=ax)

        q, _, _ = self.elem_quality(elem_ids=elem_ids)

        bins = np.linspace(0, 1, 21)
        norm = Normalize(vmin=0, vmax=1)
        colors = matplotlib.colormaps[cmap + "_r"](norm(bins[:-1]))

        for i, bin_lo in enumerate(bins[:-1]):
            bin_hi = bins[i + 1]
            mask = (q >= bin_lo) & (q < bin_hi)
            if not np.any(mask):
                continue
            self._plot_polys(elem_ids[mask], facecolor=colors[i],
                             edgecolor='k', linewidth=0.5, ax=ax)

        sm = cm.ScalarMappable(norm=norm,
                               cmap=matplotlib.colormaps[cmap + "_r"])
        sm.set_array([])
        plt.colorbar(sm, ax=ax, label='Element Quality')
        ax.set_title("Element Quality")
        return fig, ax

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _plot_polys(self, elem_ids: np.ndarray, facecolor, edgecolor='k',
                    linewidth=1.0, linestyle='-', alpha=1.0,
                    ax: plt.Axes = None) -> None:
        """
        Render elements as filled polygons via a single PolyCollection.

        Handles mixed tri/quad meshes by collecting each element's vertex
        coords as a variable-length polygon (PolyCollection supports this).
        """
        if len(elem_ids) == 0:
            return

        tri_ids, quad_ids = self._elem_type(elem_ids)
        polys = []

        for eid in tri_ids:
            verts = self.connectivity_list[eid]
            if self.type != "Triangular":
                verts = verts[:3]
            polys.append(self.points[verts, :2])

        for eid in quad_ids:
            verts = self.connectivity_list[eid, :4]
            polys.append(self.points[verts, :2])

        if not polys:
            return

        pc = PolyCollection(polys, facecolors=facecolor, edgecolors=edgecolor,
                            linewidths=linewidth, linestyles=linestyle,
                            alpha=alpha)
        ax.add_collection(pc)

    def _ensure_array(self, maybe_scalar):
        return np.array([maybe_scalar]) if np.isscalar(maybe_scalar) else maybe_scalar
