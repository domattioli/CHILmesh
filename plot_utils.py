import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

class CHILmeshPlotMixin:
    def axis_chilmesh(self, ax=None):
        if ax is None:
            ax = plt.gca()
        ax.set_aspect('equal')
        ax.set_facecolor('white')
        ax.tick_params(labelsize=12, width=1.5)
        x = self.points[:, 0]
        y = self.points[:, 1]
        offset = 0.01 * max(x.max() - x.min(), y.max() - y.min())
        ax.set_xlim([x.min() - offset, x.max() + offset])
        ax.set_ylim([y.min() - offset, y.max() + offset])
        return ax

    def plot(self, elem_ids=None, elem_color='none', edge_color='k', linewidth=1.0, linestyle='-'):
        if elem_ids is None:
            elem_ids = np.arange(self.nElems)
        ax = self.axis_chilmesh()
        if elem_color == 'none':
            self.plot_edge(self.elem2edge(elem_ids).flatten(), color=edge_color, linewidth=linewidth, linestyle=linestyle)
        else:
            tris, quads = self.elem_type(elem_ids)
            for tri in tris:
                coords = self.points[self.connectivity[tri, 1:4]]
                ax.fill(coords[:, 0], coords[:, 1], facecolor=elem_color, edgecolor=edge_color, linewidth=linewidth, linestyle=linestyle)
            for quad in quads:
                coords = self.points[self.connectivity[quad]]
                ax.fill(coords[:, 0], coords[:, 1], facecolor=elem_color, edgecolor=edge_color, linewidth=linewidth, linestyle=linestyle)
        plt.show()

    def plot_edge(self, edge_ids=None, color='k', linewidth=1.0, linestyle='-'):
        if edge_ids is None:
            edge_ids = np.arange(self.nEdges)
        v = self.edge2vert(edge_ids)
        p1 = self.points[v[:, 0]]
        p2 = self.points[v[:, 1]]
        for x1, x2 in zip(p1, p2):
            plt.plot([x1[0], x2[0]], [x1[1], x2[1]], color=color, linewidth=linewidth, linestyle=linestyle)

    def plot_elem(self, elem_ids=None, color='b', edge_color='k', linewidth=1.0, linestyle='-'):
        if elem_ids is None:
            elem_ids = np.arange(self.nElems)
        ax = self.axis_chilmesh()
        tris, quads = self.elem_type(elem_ids)
        for tri in tris:
            coords = self.points[self.connectivity[tri, 1:4]]
            ax.fill(coords[:, 0], coords[:, 1], facecolor=color, edgecolor=edge_color, linewidth=linewidth, linestyle=linestyle)
        for quad in quads:
            coords = self.points[self.connectivity[quad]]
            ax.fill(coords[:, 0], coords[:, 1], facecolor=color, edgecolor=edge_color, linewidth=linewidth, linestyle=linestyle)
        plt.show()

    def plot_point(self, ids=None, point_type='vertex', color='r', marker='o', size=5):
        if point_type.lower() in {'vertex', 'vert'}:
            if ids is None:
                ids = np.arange(self.nVerts)
            x, y, _ = self.vert_coordinates(ids)
        elif point_type.lower() in {'edge', 'midpoint'}:
            if ids is None:
                ids = np.arange(self.nEdges)
            x, y, _ = self.edge_midpoint(ids)
        elif point_type.lower() in {'element', 'centroid'}:
            if ids is None:
                ids = np.arange(self.nElems)
            x, y, _ = self.centroid(ids)
        else:
            raise ValueError(f"Unknown point type: {point_type}")
        plt.plot(x, y, linestyle='none', marker=marker, color=color, markersize=size)

    def plot_label(self, ids=None, label='all'):
        ax = self.axis_chilmesh()
        if ids is None:
            ids = np.arange(self.nElems)
        if label in {'vertex', 'point', 'all'}:
            for i in ids:
                x, y, _ = self.vert_coordinates([i])
                ax.text(x[0], y[0], f'V{i}', color='red', ha='center')
        if label in {'edge', 'all'}:
            for i in ids:
                x, y, _ = self.edge_midpoint([i])
                ax.text(x[0], y[0], f'E{i}', color='green', ha='center')
        if label in {'element', 'centroid', 'all'}:
            x, y, _ = self.centroid(ids)
            for i, (xi, yi) in zip(ids, zip(x, y)):
                ax.text(xi, yi, f'E{i}', color='blue', ha='center')

    def plot_layer(self, layers=None, cmap='viridis'):
        if layers is None:
            layers = np.arange(self.nLayers)
        norm = plt.Normalize(vmin=0, vmax=len(layers))
        colors = cm.get_cmap(cmap)(norm(np.arange(len(layers))))
        for i, idx in enumerate(layers):
            ids = self.Layers['OE'][idx] + self.Layers['IE'][idx]
            self.plot_elem(ids, color=colors[i], edge_color='k', linewidth=0.5)
        plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), label='Layer')

    def plot_quality(self, elem_ids=None, cmap='cool'):
        if elem_ids is None:
            elem_ids = np.arange(self.nElems)
        q, _ = self.elem_quality(elem_ids)
        bins = np.linspace(0, 1, 21)
        norm = plt.Normalize(vmin=0, vmax=1)
        colors = cm.get_cmap(cmap)(norm(bins[:-1]))
        for i in range(len(bins) - 1):
            bin_ids = elem_ids[(q >= bins[i]) & (q < bins[i + 1])]
            if len(bin_ids):
                self.plot_elem(bin_ids, color=colors[i], edge_color='k', linewidth=0.5)
        plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), label='Element Quality')
