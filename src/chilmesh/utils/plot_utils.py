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
        ax.autoscale_view()  # add_collection doesn't auto-update limits

        return fig, ax

    def plot_boundary(self, color='r', linewidth=2.5, linestyle='-',
                      ax=None) -> Tuple[plt.Figure, plt.Axes]:
        """Plot mesh exterior (boundary) edges highlighted."""
        return self.plot_edge(self.boundary_edges(), color=color,
                              linewidth=linewidth, linestyle=linestyle, ax=ax)

    def plot_interior_edges(self, color='b', linewidth=1.0, linestyle='-',
                            ax=None) -> Tuple[plt.Figure, plt.Axes]:
        """Plot interior (shared by two elements) edges."""
        interior = np.setdiff1d(np.arange(self.n_edges), self.boundary_edges())
        return self.plot_edge(interior, color=color, linewidth=linewidth,
                              linestyle=linestyle, ax=ax)

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
            ids = np.atleast_1d(np.asarray(ids))
            conn = self.connectivity_list[ids]          # (n, 3|4)
            # Vectorised gather: for mixed meshes the padded vertex (repeated)
            # contributes twice to the mean, giving a very slight bias toward
            # v0, but for visualisation this is negligible and avoids a loop.
            coords = self.points[conn, :2]              # (n, n_cols, 2)
            x = coords[:, :, 0].mean(axis=1)
            y = coords[:, :, 1].mean(axis=1)

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
        self.axis_chilmesh(ax=ax)  # sets equal aspect + mesh-fitted limits

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
            if i == len(bins) - 2:
                mask = (q >= bin_lo) & (q <= bin_hi)
            else:
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
        Render elements as filled polygons via PolyCollection.

        Pure-triangle and pure-quad meshes take a single vectorised numpy
        gather (no Python loop).  Mixed meshes do two batched gathers and
        combine the results into a list — still O(1) matplotlib draw calls.
        """
        if len(elem_ids) == 0:
            return

        n_cols = self.connectivity_list.shape[1]

        if n_cols == 3:
            # All-triangle mesh: single vectorised gather -> 3-D array
            verts = self.connectivity_list[elem_ids]       # (n, 3)
            polys = self.points[verts, :2]                 # (n, 3, 2)
            pc = PolyCollection(polys, facecolors=facecolor, edgecolors=edgecolor,
                                linewidths=linewidth, linestyles=linestyle,
                                alpha=alpha)
            ax.add_collection(pc)
            ax.autoscale_view()
            return

        # 4-column connectivity: detect padded triangles vs real quads
        rows = self.connectivity_list[elem_ids]            # (n, 4)
        tri_mask = (
            (rows[:, 0] == rows[:, 1])
            | (rows[:, 1] == rows[:, 2])
            | (rows[:, 2] == rows[:, 3])
            | (rows[:, 3] == rows[:, 0])
            | (rows[:, 0] == rows[:, 2])
            | (rows[:, 1] == rows[:, 3])
        )
        quad_mask = ~tri_mask

        if not quad_mask.any():
            # All triangles in a 4-col connectivity table
            polys = self.points[rows[:, :3], :2]           # (n, 3, 2)
            pc = PolyCollection(polys, facecolors=facecolor, edgecolors=edgecolor,
                                linewidths=linewidth, linestyles=linestyle,
                                alpha=alpha)
            ax.add_collection(pc)
            ax.autoscale_view()
            return

        if not tri_mask.any():
            # All quads
            polys = self.points[rows, :2]                  # (n, 4, 2)
            pc = PolyCollection(polys, facecolors=facecolor, edgecolors=edgecolor,
                                linewidths=linewidth, linestyles=linestyle,
                                alpha=alpha)
            ax.add_collection(pc)
            ax.autoscale_view()
            return

        # Mixed tri+quad: two batched gathers, polygon sizes differ so we
        # must build a Python list, but gathering is vectorised.
        all_polys = []
        if tri_mask.any():
            all_polys.extend(self.points[rows[tri_mask, :3], :2])   # (n_tri, 3, 2)
        if quad_mask.any():
            all_polys.extend(self.points[rows[quad_mask], :2])      # (n_quad, 4, 2)

        pc = PolyCollection(all_polys, facecolors=facecolor, edgecolors=edgecolor,
                            linewidths=linewidth, linestyles=linestyle,
                            alpha=alpha)
        ax.add_collection(pc)
        ax.autoscale_view()

    def _ensure_array(self, maybe_scalar):
        return np.array([maybe_scalar]) if np.isscalar(maybe_scalar) else maybe_scalar
