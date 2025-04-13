# plot_utils.py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap, BoundaryNorm
from abc import abstractproperty
from typing import Optional, List, Union, Tuple, Any, Dict

class CHILmeshPlotMixin:
    """
    A mixin class that provides plotting functionality for CHILmesh objects.
    """
    @property
    def grid_name( self ) -> str:
        """
        Return the name of the grid.
        """
        return getattr( self, '_grid_name', 'CHILmesh' )
    
    def axis_chilmesh( self, ax: Optional[plt.Axes] = None ) -> plt.Axes:
        """
        Set up the axes for plotting the mesh.
        
        Parameters:
            ax: Matplotlib axes to set up. If None, the current axes are used.
        
        Returns:
            The configured axes
        """
        if ax is None:
            ax = plt.gca()
        
        ax.set_aspect( 'equal' )
        ax.set_facecolor( 'white' )
        ax.tick_params( labelsize=12, width=1.5 )
        
        x = self.points[:, 0]
        y = self.points[:, 1]
        offset = 0.01 * max( x.max() - x.min(), y.max() - y.min() )
        
        ax.set_xlim( [x.min() - offset, x.max() + offset] )
        ax.set_ylim( [y.min() - offset, y.max() + offset] )
        
        return ax


    def plot(self, elem_ids=None, elem_color='none', edge_color='k',
            linewidth=1.0, linestyle='-', ax=None) -> Tuple[plt.Figure, plt.Axes]:
        """
        Plot the mesh.

        Parameters:
            elem_ids: IDs of elements to plot. If None, all elements are plotted.
            elem_color: Color of elements (set to 'none' for edge-only).
            edge_color: Color of element or edge outlines.
            linewidth: Line width.
            linestyle: Line style.
            ax: Optional matplotlib axis. Creates a new one if not provided.

        Returns:
            (fig, ax): matplotlib Figure and Axes
        """
        if elem_ids is None:
            elem_ids = np.arange(self.n_elems)
        elif np.isscalar(elem_ids):
            elem_ids = np.array([elem_ids])

        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 6))
        else:
            fig = ax.figure

        ax.set_aspect('equal')
        ax.set_facecolor('white')
        ax.tick_params(labelsize=12, width=1.5)

        x = self.points[:, 0]
        y = self.points[:, 1]
        offset = 0.01 * max(x.max() - x.min(), y.max() - y.min())
        ax.set_xlim([x.min() - offset, x.max() + offset])
        ax.set_ylim([y.min() - offset, y.max() + offset])

        if elem_color == 'none':
            edges = self.elem2edge(elem_ids).flatten()
            edges = edges[edges > 0]
            self.plot_edge(edges, color=edge_color, linewidth=linewidth, linestyle=linestyle, ax=ax)
        else:
            tri_elems, quad_elems = self._elem_type(elem_ids)

            for elem_id in tri_elems:
                vertices = self.connectivity_list[elem_id]
                if self.type != "Triangular":
                    vertices = vertices[:3]
                coords = self.points[vertices, :2]
                tri = plt.Polygon(coords, facecolor=elem_color, edgecolor=edge_color,
                                linewidth=linewidth, linestyle=linestyle)
                ax.add_patch(tri)

            for elem_id in quad_elems:
                vertices = self.connectivity_list[elem_id]
                coords = self.points[vertices, :2]
                quad = plt.Polygon(coords, facecolor=elem_color, edgecolor=edge_color,
                                linewidth=linewidth, linestyle=linestyle)
                ax.add_patch(quad)

        return fig, ax

    def plot_edge(self, edge_ids=None, color='g', linewidth=2.5, linestyle='-', ax=None):
        if edge_ids is None:
            edge_ids = np.arange(self.n_edges)

        if ax is None:
            ax = self.axis_chilmesh()
        fig = ax.figure


        v = self.edge2vert(edge_ids)
        p1 = self.points[v[:, 0], :2]
        p2 = self.points[v[:, 1], :2]

        for (x1, y1), (x2, y2) in zip(p1, p2):
            ax.plot([x1, x2], [y1, y2], color=color, linewidth=linewidth, linestyle=linestyle)

        return fig, ax

    def plot_elem(self, elem_ids=None, color='b', edge_color='k', linewidth=1.0, linestyle='-', ax=None):
        elem_ids = self._ensure_array(elem_ids)
        if elem_ids is None:
            elem_ids = np.arange(self.n_elems)

        if ax is None:
            fig, ax = self.plot()
        else:
            fig = ax.figure

        tri_elems, quad_elems = self._elem_type(elem_ids)

        for elem_id in tri_elems:
            vertices = self.connectivity_list[elem_id]
            if self.type != "Triangular":
                vertices = vertices[:3]
            coords = self.points[vertices, :2]
            tri = plt.Polygon(coords, facecolor=color, edgecolor=edge_color, linewidth=linewidth, linestyle=linestyle)
            ax.add_patch(tri)

        for elem_id in quad_elems:
            vertices = self.connectivity_list[elem_id]
            coords = self.points[vertices, :2]
            quad = plt.Polygon(coords, facecolor=color, edgecolor=edge_color, linewidth=linewidth, linestyle=linestyle)
            ax.add_patch(quad)

        return fig, ax
    
    def plot_face(self, face_ids=None, color='b', edge_color='k', linewidth=1.0, linestyle='-', ax=None):
        """
        Plot faces of the mesh (alias for plot_elem).
        
        Parameters:
            face_ids: IDs of faces to plot. If None, all faces are plotted.
            color: Color of faces.
            edge_color: Color of edges.
            linewidth: Width of edge lines.
            linestyle: Style of edge lines.
            ax: Optional matplotlib axis to plot on.
            
        Returns:
            A tuple containing the figure and axis objects.
        """
        return self.plot_elem(elem_ids=face_ids, color=color, edge_color=edge_color, 
                             linewidth=linewidth, linestyle=linestyle, ax=ax)

    def plot_point( self, ids: Optional[np.ndarray] = None, point_type: str = 'vertex',
                   color: str = 'r', marker: str = 'o', size: float = 5, ax: plt.Axes=None, **kwargs ) -> Tuple[plt.Figure, plt.Axes]:
        """
        Plot points of the mesh.
        
        Parameters:
            ids: IDs of points to plot. If None, all points of the specified type are plotted.
            point_type: Type of points to plot ('vertex', 'edge', or 'element').
            color: Color of points.
            marker: Marker style.
            size: Marker size.
        """
        # ax = self.axis_chilmesh()
        if ax is None:
            fig, ax = self.plot()  # fallback
        else:
            fig = ax.figure
        
        if point_type.lower() in {'vertex', 'vert'}:
            if ids is None:
                ids = np.arange( self.n_verts )
            x, y = self.points[ids, 0], self.points[ids, 1]
            
        elif point_type.lower() in {'edge', 'midpoint'}:
            if ids is None:
                ids = np.arange( self.n_edges )
            
            # Get vertices for each edge
            edges = self.edge2vert( ids )
            
            # Calculate midpoints
            x = np.mean( self.points[edges, 0], axis=1 )
            y = np.mean( self.points[edges, 1], axis=1 )
            
        elif point_type.lower() in {'element', 'centroid'}:
            if ids is None:
                ids = np.arange( self.n_elems )
            
            # Calculate centroids
            centroids = np.zeros( ( len( ids ), 2 ) )
            for i, elem_id in enumerate( ids ):
                vertices = self.connectivity_list[elem_id]
                vertices = vertices[vertices > 0]  # Remove zero vertices
                centroids[i, 0] = np.mean( self.points[vertices, 0] )
                centroids[i, 1] = np.mean( self.points[vertices, 1] )
            
            x, y = centroids[:, 0], centroids[:, 1]
            
        else:
            raise ValueError( f"Unknown point type: {point_type}" )
        
        ax.plot( x, y, linestyle='none', marker=marker, color=color, markersize=size, **kwargs )
        return fig, ax
    
    def plot_label(self, ids: Optional[np.ndarray] = None, label: str = 'all', ax: Optional[plt.Axes] = None) -> Tuple[plt.Figure, plt.Axes]:
        """
        Label mesh entities.

        Parameters:
            ids: IDs of entities to label. If None, all entities of the specified type are labeled.
            label: Type of entities to label ('vertex', 'edge', 'element', or 'all').
            ax: Optional matplotlib axis to plot on.
        """
        if ax is None:
            fig, ax = self.plot()
        else:
            fig = ax.figure

        if ids is None:
            ids = np.arange(self.n_elems)

        if label in {'vertex', 'point', 'all'}:
            for i in ids:
                x, y = self.points[i, 0], self.points[i, 1]
                ax.text(x, y, f'V{i}', color='red', ha='center')

        if label in {'edge', 'all'}:
            for i in ids:
                if i < self.n_edges:
                    edges = self.edge2vert([i])
                    x = np.mean(self.points[edges, 0])
                    y = np.mean(self.points[edges, 1])
                    ax.text(x, y, f'E{i}', color='green', ha='center')

        if label in {'element', 'centroid', 'all'}:
            for i in ids:
                if i < self.n_elems:
                    vertices = self.connectivity_list[i]
                    vertices = vertices[vertices > 0]
                    x = np.mean(self.points[vertices, 0])
                    y = np.mean(self.points[vertices, 1])
                    ax.text(x, y, f'E{i}', color='blue', ha='center')
        return fig, ax
    
    def plot_layer(self, layers=None, cmap='viridis', ax=None) -> Tuple[plt.Figure, plt.Axes]:
        """
        Plot the specified mesh layers with different colors.

        Parameters:
            layers: Indices of layers to plot. If None, all layers are plotted.
            cmap: Colormap to use for the layers.
            ax: Optional matplotlib axis to plot on.
        
        Returns:
            (fig, ax): Matplotlib Figure and Axes
        """
        if layers is None:
            layers = range(self.n_layers)

        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 8))
        else:
            fig = ax.figure

        ax.set_aspect('equal')

        cmap_obj = cm.get_cmap(cmap, self.n_layers)
        norm = BoundaryNorm(boundaries=np.arange(self.n_layers + 1), ncolors=self.n_layers)

        for layer_idx in layers:
            if layer_idx >= len(self.layers["OE"]) or layer_idx >= len(self.layers["IE"]):
                continue

            outer_elems = self.layers["OE"][layer_idx]
            inner_elems = self.layers["IE"][layer_idx]
            elem_ids = np.concatenate((outer_elems, inner_elems)).astype(int)

            color = cmap_obj(norm(layer_idx))

            for elem_id in elem_ids:
                if elem_id < 0 or elem_id >= len(self.connectivity_list):
                    continue
                vertices = self.connectivity_list[elem_id]
                valid_indices = [v for v in vertices if v > 0 and v < len(self.points)]
                if len(valid_indices) >= 3:
                    polygon = self.points[valid_indices, :2]
                    ax.fill(polygon[:, 0], polygon[:, 1], color=color, edgecolor='k',
                            linewidth=0.5, alpha=0.7)

        sm = cm.ScalarMappable(norm=norm, cmap=cmap_obj)
        sm.set_array([])
        plt.colorbar(sm, ax=ax, label='Layer',
                    ticks=np.arange(0, self.n_layers) + 0.5,
                    format='%d')

        ax.set_title(f"Mesh Layers ({self.n_layers} total)")
        return fig, ax

    def plot_quality(self, elem_ids=None, cmap='cool', ax=None) -> Tuple[plt.Figure, plt.Axes]:
        """
        Plot element quality as a contour map.

        Parameters:
            elem_ids: IDs of elements to plot. If None, all elements are plotted.
            cmap: Colormap to use.
            ax: Optional matplotlib axis to plot on.

        Returns:
            (fig, ax): Matplotlib Figure and Axes
        """
        if elem_ids is None:
            elem_ids = np.arange(self.n_elems)
        elif np.isscalar(elem_ids):
            elem_ids = np.array([elem_ids])

        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 8))
        else:
            fig = ax.figure
        self.axis_chilmesh(ax=ax)  # ensures consistent scaling and view

        # Compute element quality (returns q, angles, stats)
        q, _, _ = self.elem_quality(elem_ids=elem_ids)

        # Bin and color by quality
        bins = np.linspace(0, 1, 21)
        norm = plt.Normalize(vmin=0, vmax=1)
        colors = cm.get_cmap(cmap + "_r")(norm(bins[:-1]))
        for i in range(len(bins) - 1):
            bin_ids = elem_ids[(q >= bins[i]) & (q < bins[i + 1])]
            if len(bin_ids):
                self.plot_elem(bin_ids, color=colors[i], edge_color='k', linewidth=0.5, ax=ax)

        sm = cm.ScalarMappable(norm=norm, cmap=cmap + "_r")
        sm.set_array([])
        plt.colorbar(sm, ax=ax, label='Element Quality')
        ax.set_title("Element Quality")
        return fig, ax


    def _ensure_array(self, maybe_scalar):
        return np.array([maybe_scalar]) if np.isscalar(maybe_scalar) else maybe_scalar
