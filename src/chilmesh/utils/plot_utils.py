from __future__ import annotations

# plot_utils.py
#
# Thin mesh-bound wrappers over the standalone ``chilmesh.chilplotting`` module.
# The rendering logic lives in ``chilplotting`` so it can be reused by any mesh
# generator (ADMESH / QuADMesh / MADMESHing) on raw point + connectivity arrays
# (CHILmesh #75). These methods just forward ``self`` to the module functions,
# preserving the historical ``CHILmesh.plot*`` API.
import matplotlib.pyplot as plt
from typing import Optional, Tuple

from .. import chilplotting as _cp


class CHILmeshPlotMixin:
    """Mixin providing fast, vectorized plotting for CHILmesh objects.

    All methods delegate to :mod:`chilmesh.chilplotting`; see that module for
    the generator-agnostic array primitives (``plot_mesh``, ``plot_filled``, …).
    """

    def axis_chilmesh(self, ax: Optional[plt.Axes] = None) -> plt.Axes:
        """Configure axes with equal aspect and mesh-fitted limits."""
        return _cp.axis_chilmesh(self, ax=ax)

    def plot(self, elem_ids=None, elem_color='none', edge_color='k',
             linewidth=1.0, linestyle='-', ax=None) -> Tuple[plt.Figure, plt.Axes]:
        """Plot the mesh (edge-only when ``elem_color='none'``, else filled)."""
        return _cp.plot(self, elem_ids=elem_ids, elem_color=elem_color,
                        edge_color=edge_color, linewidth=linewidth,
                        linestyle=linestyle, ax=ax)

    def plot_edge(self, edge_ids=None, color='g', linewidth=2.5,
                  linestyle='-', ax=None) -> Tuple[plt.Figure, plt.Axes]:
        """Plot edges using a single vectorized LineCollection call."""
        return _cp.plot_edge(self, edge_ids=edge_ids, color=color,
                             linewidth=linewidth, linestyle=linestyle, ax=ax)

    def plot_boundary(self, color='r', linewidth=2.5, linestyle='-',
                      ax=None) -> Tuple[plt.Figure, plt.Axes]:
        """Plot mesh exterior (boundary) edges highlighted."""
        return _cp.plot_boundary(self, color=color, linewidth=linewidth,
                                 linestyle=linestyle, ax=ax)

    def plot_interior_edges(self, color='b', linewidth=1.0, linestyle='-',
                            ax=None) -> Tuple[plt.Figure, plt.Axes]:
        """Plot interior (shared by two elements) edges."""
        return _cp.plot_interior_edges(self, color=color, linewidth=linewidth,
                                       linestyle=linestyle, ax=ax)

    def plot_elem(self, elem_ids=None, color='b', edge_color='k',
                  linewidth=1.0, linestyle='-',
                  ax=None) -> Tuple[plt.Figure, plt.Axes]:
        """Plot elements filled with a single color using PolyCollection."""
        return _cp.plot_elem(self, elem_ids=elem_ids, color=color,
                             edge_color=edge_color, linewidth=linewidth,
                             linestyle=linestyle, ax=ax)

    def plot_face(self, face_ids=None, color='b', edge_color='k',
                  linewidth=1.0, linestyle='-',
                  ax=None) -> Tuple[plt.Figure, plt.Axes]:
        """Plot faces (alias for plot_elem)."""
        return _cp.plot_face(self, face_ids=face_ids, color=color,
                             edge_color=edge_color, linewidth=linewidth,
                             linestyle=linestyle, ax=ax)

    def plot_point(self, ids=None, point_type: str = 'vertex', color: str = 'r',
                   marker: str = 'o', size: float = 5,
                   ax: plt.Axes = None, **kwargs) -> Tuple[plt.Figure, plt.Axes]:
        """Plot mesh entities as scatter points (vertices, edge midpoints, centroids)."""
        return _cp.plot_point(self, ids=ids, point_type=point_type, color=color,
                              marker=marker, size=size, ax=ax, **kwargs)

    def plot_label(self, ids=None, label: str = 'all',
                   ax: Optional[plt.Axes] = None) -> Tuple[plt.Figure, plt.Axes]:
        """Annotate mesh entities with their IDs."""
        return _cp.plot_label(self, ids=ids, label=label, ax=ax)

    def plot_layer(self, layers=None, cmap='viridis',
                   ax=None) -> Tuple[plt.Figure, plt.Axes]:
        """Plot skeletonization layers as filled elements colored by layer index."""
        return _cp.plot_layer(self, layers=layers, cmap=cmap, ax=ax)

    def plot_quality(self, elem_ids=None, cmap='cool',
                     ax=None) -> Tuple[plt.Figure, plt.Axes]:
        """Plot element quality as a filled colormap."""
        return _cp.plot_quality(self, elem_ids=elem_ids, cmap=cmap, ax=ax)

    def plot_quality_histogram(self, elem_ids=None, bins: int = 100,
                               cmap: str = 'cool', auto_norm: bool = False,
                               ax: Optional[plt.Axes] = None
                               ) -> Tuple[plt.Figure, plt.Axes]:
        """Plot a histogram of element quality coloured by ``cmap``."""
        return _cp.plot_quality_from_mesh_histogram(
            self, elem_ids=elem_ids, bins=bins, cmap=cmap,
            auto_norm=auto_norm, ax=ax)
