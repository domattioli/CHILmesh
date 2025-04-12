import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import os
import re
from typing import List, Dict, Tuple, Union, Optional, Any, Set


class CHILmesh:
    """
    CHILmesh: Triangulation, Quadrangulation, or Mixed-Elements in 2-D.
    
    CHILmesh supports topological and geometric queries for 2D
    triangular, quadrangular, and mixed finite element meshes. CHILmesh
    provides methods for analyzing the element nodal connectivity list.
    Various methods exist that allow for queries involving the geometry and
    topology of vertices, edges, and elements. Additional methods allow for
    plotting the mesh.
    
    Note: "vertex", "node", and "point" are used interchangeably throughout 
    this class's documentation.
    """

    def __init__(self, *args):
        """
        Initialize a CHILmesh object from different input types.
        
        Parameters
        ----------
        *args : varies
            Different initialization methods are supported:
            - filename (str): Load mesh from a file (.msh or .14 formats)
            - connectivity_list, points: Create mesh from connectivity and points
            - connectivity_list, points, boundary_condition: Create mesh with BC
        """
        # Initialize properties
        self.grid_name = None
        self.points = None
        self.connectivity_list = None
        self.boundary_condition = None
        self.mesh_type = None
        
        # Hidden properties
        self.adjacencies = None
        self.n_verts = 0
        self.n_elems = 0
        self.n_edges = 0
        self.layers = None
        self.n_layers = 0
        self.filename = None
        
        # Process input arguments
        self._process_inputs(args)
        
        # Setup mesh
        if self.connectivity_list is not None and self.points is not None:
            self.n_verts = self.points.shape[0]
            self.n_elems = self.connectivity_list.shape[0]
            
            # Compute signed area using the shoelace formula
            areas[quad_mask] = 0.5 * (
                (x1 * y2 - x2 * y1) +
                (x2 * y3 - x3 * y2) +
                (x3 * y4 - x4 * y3) +
                (x4 * y1 - x1 * y4)
            )
        
        # Return unsigned area if requested
        if sign.lower() == 'unsigned':
            return np.abs(areas)
        return areas

    def elem_type(self, elem_ids=None):
        """
        Identify triangular and quadrilateral elements.
        
        Parameters
        ----------
        elem_ids : array_like, optional
            Element IDs to check. If None, check all elements.
            
        Returns
        -------
        tri_ids : ndarray
            Indices of triangular elements
        quad_ids : ndarray
            Indices of quadrilateral elements
        """
        if elem_ids is None:
            elem_ids = np.arange(self.n_elems)
        else:
            elem_ids = np.array(elem_ids)
        
        if self.connectivity_list.shape[1] == 3:  # All triangles
            return elem_ids, np.array([], dtype=np.int32)
        
        # Mixed mesh - identify triangles
        tri_mask = np.zeros(len(elem_ids), dtype=bool)
        
        for i, elem_id in enumerate(elem_ids):
            conn = self.connectivity_list[elem_id]
            # Check for triangles (identified by repeated vertices or 0 value)
            if (conn[0] == conn[1] or 
                conn[1] == conn[2] or 
                conn[2] == conn[3] or 
                conn[3] == conn[0] or
                conn[3] == 0):
                tri_mask[i] = True
        
        tri_ids = elem_ids[tri_mask]
        quad_ids = elem_ids[~tri_mask]
        
        return tri_ids, quad_ids

    def boundary_edges(self, bc_ids=None):
        """
        Find edges that define the boundary of the mesh.
        
        Parameters
        ----------
        bc_ids : array_like, optional
            Boundary condition IDs to filter by
            
        Returns
        -------
        edge_ids : ndarray
            Indices of boundary edges
        i_bc : ndarray, optional
            Boundary condition indices for each edge
        """
        # Edges with only one neighboring element are boundary edges
        edge_ids = np.where(np.any(self.adjacencies['edge2elem'] == -1, axis=1))[0]
        
        # For the second output (boundary condition indices)
        if len(edge_ids) > 0:
            i_bc = np.ones(len(edge_ids), dtype=np.int32)  # Placeholder
            return edge_ids, i_bc
        else:
            return edge_ids

    def centroid(self, elem_ids=None):
        """
        Calculate the centroid of elements.
        
        Parameters
        ----------
        elem_ids : array_like, optional
            Element IDs to calculate centroids for. If None, calculate for all elements.
            
        Returns
        -------
        x, y, z : ndarray
            Coordinates of centroids
        """
        if elem_ids is None:
            elem_ids = np.arange(self.n_elems)
        else:
            elem_ids = np.array(elem_ids)
        
        # Get triangles and quads
        tri_ids, quad_ids = self.elem_type(elem_ids)
        
        # Initialize output arrays
        x = np.zeros(len(elem_ids))
        y = np.zeros(len(elem_ids))
        z = np.zeros(len(elem_ids))
        
        # Process triangles
        if len(tri_ids) > 0:
            tri_mask = np.isin(elem_ids, tri_ids)
            tri_indices = elem_ids[tri_mask]
            
            # Get vertices
            conn = self.connectivity_list[tri_indices]
            
            # For triangles, use first 3 vertices (in case of mixed mesh)
            if conn.shape[1] > 3:
                conn = conn[:, :3]
            
            # Calculate centroids by averaging vertices
            v_coords = np.array([self.points[v] for v in conn.flat]).reshape(len(tri_indices), 3, 3)
            centroids = np.mean(v_coords, axis=1)
            
            x[tri_mask] = centroids[:, 0]
            y[tri_mask] = centroids[:, 1]
            z[tri_mask] = centroids[:, 2]
        
        # Process quads
        if len(quad_ids) > 0:
            quad_mask = np.isin(elem_ids, quad_ids)
            quad_indices = elem_ids[quad_mask]
            
            # Get vertices
            conn = self.connectivity_list[quad_indices]
            
            # Calculate centroids by averaging vertices
            v_coords = np.array([self.points[v] for v in conn.flat]).reshape(len(quad_indices), 4, 3)
            centroids = np.mean(v_coords, axis=1)
            
            x[quad_mask] = centroids[:, 0]
            y[quad_mask] = centroids[:, 1]
            z[quad_mask] = centroids[:, 2]
        
        if len(elem_ids) == 1:
            return x[0], y[0], z[0]
        
        return x, y, z
    
    def edge2vert(self, edge_ids=None):
        """
        Get vertices defining the specified edges.
        
        Parameters
        ----------
        edge_ids : array_like, optional
            Edge IDs to get vertices for. If None, get for all edges.
            
        Returns
        -------
        vert_ids : ndarray
            Array of vertex ID pairs defining each edge
        """
        if edge_ids is None:
            return self.adjacencies['edge2vert']
        else:
            edge_ids = np.array(edge_ids)
            return self.adjacencies['edge2vert'][edge_ids]
    
    def edge2elem(self, edge_ids=None):
        """
        Get elements attached to the specified edges.
        
        Parameters
        ----------
        edge_ids : array_like, optional
            Edge IDs to get elements for. If None, get for all edges.
            
        Returns
        -------
        elem_ids : ndarray
            Array of element ID pairs attached to each edge
        """
        if edge_ids is None:
            return self.adjacencies['edge2elem']
        else:
            edge_ids = np.array(edge_ids)
            return self.adjacencies['edge2elem'][edge_ids]
    
    def elem2edge(self, elem_ids=None):
        """
        Get edges defining the specified elements.
        
        Parameters
        ----------
        elem_ids : array_like, optional
            Element IDs to get edges for. If None, get for all elements.
            
        Returns
        -------
        edge_ids : ndarray
            Array of edge IDs defining each element
        """
        if elem_ids is None:
            return self.adjacencies['elem2edge']
        else:
            elem_ids = np.array(elem_ids)
            return self.adjacencies['elem2edge'][elem_ids]
    
    def elem2vert(self, elem_ids=None):
        """
        Get vertices defining the specified elements.
        
        Parameters
        ----------
        elem_ids : array_like, optional
            Element IDs to get vertices for. If None, get for all elements.
            
        Returns
        -------
        vert_ids : ndarray
            Array of vertex IDs defining each element
        """
        if elem_ids is None:
            return self.adjacencies['elem2vert']
        else:
            elem_ids = np.array(elem_ids)
            return self.connectivity_list[elem_ids]
    
    def vert2edge(self, vert_ids=None, store='matrix'):
        """
        Get edges attached to the specified vertices.
        
        Parameters
        ----------
        vert_ids : array_like, optional
            Vertex IDs to get edges for. If None, get for all vertices.
        store : str, optional
            Output format: 'matrix' or 'cell' (list of lists)
            
        Returns
        -------
        edge_ids : ndarray or list
            Edge IDs attached to each vertex
        """
        if vert_ids is None:
            vert_ids = np.arange(self.n_verts)
        else:
            vert_ids = np.array(vert_ids)
        
        if store.lower() == 'matrix':
            return self.adjacencies['vert2edge'][vert_ids]
        else:  # Return as list of lists
            result = []
            for vid in vert_ids:
                edges = self.adjacencies['vert2edge'][vid]
                edges = edges[edges >= 0]  # Filter out -1 values
                result.append(edges.tolist())
            return result
    
    def vert2elem(self, vert_ids=None, store='matrix'):
        """
        Get elements attached to the specified vertices.
        
        Parameters
        ----------
        vert_ids : array_like, optional
            Vertex IDs to get elements for. If None, get for all vertices.
        store : str, optional
            Output format: 'matrix' or 'cell' (list of lists)
            
        Returns
        -------
        elem_ids : ndarray or list
            Element IDs attached to each vertex
        """
        if vert_ids is None:
            vert_ids = np.arange(self.n_verts)
        else:
            vert_ids = np.array(vert_ids)
        
        if store.lower() == 'matrix':
            return self.adjacencies['vert2elem'][vert_ids]
        else:  # Return as list of lists
            result = []
            for vid in vert_ids:
                elems = self.adjacencies['vert2elem'][vid]
                elems = elems[elems >= 0]  # Filter out -1 values
                result.append(elems.tolist())
            return result
    
    def edge_coordinates(self, edge_ids=None):
        """
        Get coordinates of the vertices defining edges.
        
        Parameters
        ----------
        edge_ids : array_like, optional
            Edge IDs to get coordinates for. If None, get for all edges.
            
        Returns
        -------
        x, y, z : ndarray
            Coordinates of edge vertices
        """
        if edge_ids is None:
            edge_ids = np.arange(self.n_edges)
        else:
            edge_ids = np.array(edge_ids)
        
        # Get vertices defining each edge
        vert_ids = self.edge2vert(edge_ids)
        
        # Get coordinates
        x = np.zeros((len(edge_ids), 2))
        y = np.zeros((len(edge_ids), 2))
        z = np.zeros((len(edge_ids), 2))
        
        for i, (v1, v2) in enumerate(vert_ids):
            x[i, 0] = self.points[v1, 0]
            x[i, 1] = self.points[v2, 0]
            y[i, 0] = self.points[v1, 1]
            y[i, 1] = self.points[v2, 1]
            z[i, 0] = self.points[v1, 2]
            z[i, 1] = self.points[v2, 2]
        
        return x, y, z
    
    def edge_length(self, edge_ids=None):
        """
        Calculate the length of edges.
        
        Parameters
        ----------
        edge_ids : array_like, optional
            Edge IDs to calculate lengths for. If None, calculate for all edges.
            
        Returns
        -------
        lengths : ndarray
            Lengths of the specified edges
        """
        if edge_ids is None:
            edge_ids = np.arange(self.n_edges)
        else:
            edge_ids = np.array(edge_ids)
        
        # Get edge coordinates
        x, y, z = self.edge_coordinates(edge_ids)
        
        # Calculate distances
        dx = x[:, 1] - x[:, 0]
        dy = y[:, 1] - y[:, 0]
        dz = z[:, 1] - z[:, 0]
        
        return np.sqrt(dx**2 + dy**2 + dz**2)
    
    def edge_midpoint(self, edge_ids=None):
        """
        Calculate the midpoint of edges.
        
        Parameters
        ----------
        edge_ids : array_like, optional
            Edge IDs to calculate midpoints for. If None, calculate for all edges.
            
        Returns
        -------
        x, y, z : ndarray
            Coordinates of edge midpoints
        """
        if edge_ids is None:
            edge_ids = np.arange(self.n_edges)
        else:
            edge_ids = np.array(edge_ids)
        
        # Get edge coordinates
        x, y, z = self.edge_coordinates(edge_ids)
        
        # Calculate midpoints
        x_mid = np.mean(x, axis=1)
        y_mid = np.mean(y, axis=1)
        z_mid = np.mean(z, axis=1)
        
        return x_mid, y_mid, z_mid
    
    def elem_coordinates(self, elem_ids=None):
        """
        Get coordinates of the vertices defining elements.
        
        Parameters
        ----------
        elem_ids : array_like, optional
            Element IDs to get coordinates for. If None, get for all elements.
            
        Returns
        -------
        x, y, z : ndarray
            Coordinates of element vertices
        """
        if elem_ids is None:
            elem_ids = np.arange(self.n_elems)
        else:
            elem_ids = np.array(elem_ids)
        
        # Get element connectivity
        conn = self.elem2vert(elem_ids)
        
        # Initialize output arrays
        n_cols = conn.shape[1]
        x = np.zeros((len(elem_ids), n_cols))
        y = np.zeros((len(elem_ids), n_cols))
        z = np.zeros((len(elem_ids), n_cols))
        
        # Get coordinates for each vertex
        for i, elem_conn in enumerate(conn):
            for j, v in enumerate(elem_conn):
                if v >= 0 and v < self.n_verts:  # Check for valid vertex index
                    x[i, j] = self.points[v, 0]
                    y[i, j] = self.points[v, 1]
                    z[i, j] = self.points[v, 2]
        
        return x, y, z
    
    def vert_coordinates(self, vert_ids=None):
        """
        Get coordinates of vertices.
        
        Parameters
        ----------
        vert_ids : array_like, optional
            Vertex IDs to get coordinates for. If None, get for all vertices.
            
        Returns
        -------
        x, y, z : ndarray
            Coordinates of vertices
        """
        if vert_ids is None:
            vert_ids = np.arange(self.n_verts)
        else:
            vert_ids = np.array(vert_ids, dtype=np.int32).flatten()
        
        # Handle 2D array input (e.g., from edge2vert)
        if vert_ids.ndim > 1:
            flat_ids = vert_ids.flatten()
            coords = self.points[flat_ids]
            x = coords[:, 0].reshape(vert_ids.shape)
            y = coords[:, 1].reshape(vert_ids.shape)
            z = coords[:, 2].reshape(vert_ids.shape)
        else:
            coords = self.points[vert_ids]
            x = coords[:, 0]
            y = coords[:, 1]
            z = coords[:, 2]
        
        return x, y, z
    
    def elem2elem(self, elem_ids=None):
        """
        Get neighboring elements for the specified elements.
        
        Parameters
        ----------
        elem_ids : array_like, optional
            Element IDs to get neighbors for. If None, get for all elements.
            
        Returns
        -------
        neighbor_ids : ndarray
            Element IDs of neighbors
        """
        if elem_ids is None:
            elem_ids = np.arange(self.n_elems)
        else:
            elem_ids = np.array(elem_ids)
        
        # Get edges of each element
        elem_edges = self.elem2edge(elem_ids)
        
        # Maximum number of edges per element
        max_edges = elem_edges.shape[1]
        
        # Initialize output array
        neighbors = np.zeros((len(elem_ids), max_edges), dtype=np.int32) - 1
        
        # For each element, find its neighbors through shared edges
        for i, elem_id in enumerate(elem_ids):
            edge_ids = elem_edges[i]
            
            for j, edge_id in enumerate(edge_ids):
                if edge_id >= 0:  # Valid edge
                    # Get elements attached to this edge
                    edge_elems = self.edge2elem(edge_id)
                    
                    # Add neighbors (excluding the element itself)
                    for k, neighbor_id in enumerate(edge_elems):
                        if neighbor_id >= 0 and neighbor_id != elem_id:
                            neighbors[i, j] = neighbor_id
        
        return neighbors
    
    def interior_angles(self, elem_ids=None):
        """
        Calculate interior angles of elements.
        
        Parameters
        ----------
        elem_ids : array_like, optional
            Element IDs to calculate angles for. If None, calculate for all elements.
            
        Returns
        -------
        angles : ndarray
            Interior angles of elements (in degrees)
        """
        if elem_ids is None:
            elem_ids = np.arange(self.n_elems)
        else:
            elem_ids = np.array(elem_ids)
        
        # Get triangles and quads
        tri_ids, quad_ids = self.elem_type(elem_ids)
        
        # Maximum number of vertices per element
        max_verts = self.connectivity_list.shape[1]
        
        # Initialize output array
        angles = np.zeros((len(elem_ids), max_verts))
        
        # Process triangles
        if len(tri_ids) > 0:
            tri_mask = np.isin(elem_ids, tri_ids)
            tri_indices = elem_ids[tri_mask]
            
            # Get edge lengths
            edge_ids = self.elem2edge(tri_indices).flatten()
            edge_ids = edge_ids[edge_ids >= 0]
            
            # Get edge lengths and organize by element
            edge_lengths = self.edge_length(edge_ids)
            elem_edges = {}
            for i, elem_id in enumerate(tri_indices):
                edges = self.elem2edge(elem_id)[0]
                edges = edges[edges >= 0]
                
                # Map edge IDs to their indices in edge_ids
                edge_indices = [np.where(edge_ids == edge_id)[0][0] for edge_id in edges]
                elem_edges[elem_id] = edge_lengths[edge_indices]
            
            # Calculate angles using Law of Cosines
            for i, elem_id in enumerate(tri_indices):
                if elem_id in elem_edges:
                    a, b, c = elem_edges[elem_id]
                    
                    # Calculate angles
                    angle_a = np.arccos((b**2 + c**2 - a**2) / (2 * b * c))
                    angle_b = np.arccos((a**2 + c**2 - b**2) / (2 * a * c))
                    angle_c = np.arccos((a**2 + b**2 - c**2) / (2 * a * b))
                    
                    # Convert to degrees
                    angles[tri_mask][i, :3] = np.degrees([angle_a, angle_b, angle_c])
        
        # Process quads
        if len(quad_ids) > 0:
            quad_mask = np.isin(elem_ids, quad_ids)
            quad_indices = elem_ids[quad_mask]
            
            # For each quad, calculate angles between edges
            for i, elem_id in enumerate(quad_indices):
                x, y, _ = self.elem_coordinates([elem_id])
                x, y = x[0], y[0]
                
                # For each vertex, calculate the angle
                for j in range(4):
                    # Get adjacent vertices (with wraparound)
                    prev_j = (j - 1) % 4
                    next_j = (j + 1) % 4
                    
                    # Calculate vectors
                    v1 = np.array([x[prev_j] - x[j], y[prev_j] - y[j]])
                    v2 = np.array([x[next_j] - x[j], y[next_j] - y[j]])
                    
                    # Normalize vectors
                    v1_norm = np.linalg.norm(v1)
                    v2_norm = np.linalg.norm(v2)
                    
                    if v1_norm > 0 and v2_norm > 0:
                        v1 = v1 / v1_norm
                        v2 = v2 / v2_norm
                        
                        # Calculate angle (dot product)
                        cos_angle = np.clip(np.dot(v1, v2), -1.0, 1.0)
                        angle = np.arccos(cos_angle)
                        
                        # Convert to degrees
                        angles[quad_mask][i, j] = np.degrees(angle)
        
        return angles
    
    def elem_quality(self, elem_ids=None, quality_measure='skew'):
        """
        Calculate the quality of elements.
        
        Parameters
        ----------
        elem_ids : array_like, optional
            Element IDs to calculate quality for. If None, calculate for all elements.
        quality_measure : str, optional
            Method to measure quality: 'skew' for angular skewness
            
        Returns
        -------
        quality : ndarray
            Quality measure for each element
        angles : ndarray, optional
            Interior angles if requested
        """
        if elem_ids is None:
            elem_ids = np.arange(self.n_elems)
        else:
            elem_ids = np.array(elem_ids)
        
        # Get triangles and quads
        tri_ids, quad_ids = self.elem_type(elem_ids)
        
        # Calculate interior angles
        angles = self.interior_angles(elem_ids)
        
        # Initialize quality array
        quality = np.zeros(len(elem_ids))
        
        if quality_measure.lower() in ['skew', 'skewness', 'angular skewness']:
            # Process triangles
            if len(tri_ids) > 0:
                tri_mask = np.isin(elem_ids, tri_ids)
                
                # For triangles, ideal angle is 60 degrees
                tri_angles = angles[tri_mask, :3]
                
                angle_max = np.max(tri_angles, axis=1)
                angle_min = np.min(tri_angles, axis=1)
                
                # Calculate equiangular skew
                quality[tri_mask] = 1 - np.maximum(
                    (angle_max - 60) / (180 - 60),
                    (60 - angle_min) / 60
                )
            
            # Process quads
            if len(quad_ids) > 0:
                quad_mask = np.isin(elem_ids, quad_ids)
                
                # For quads, ideal angle is 90 degrees
                quad_angles = angles[quad_mask]
                
                angle_max = np.max(quad_angles, axis=1)
                angle_min = np.min(quad_angles, axis=1)
                
                # Calculate equiangular skew
                quality[quad_mask] = 1 - np.maximum(
                    (angle_max - 90) / (180 - 90),
                    (90 - angle_min) / 90
                )
            
            # Check for invalid angles (e.g., concave elements)
            sum_angles_tri = np.sum(angles[np.isin(elem_ids, tri_ids), :3], axis=1)
            sum_angles_quad = np.sum(angles[np.isin(elem_ids, quad_ids)], axis=1)
            
            # Flag elements with invalid angle sums
            invalid_tri = np.isin(elem_ids, tri_ids) & (sum_angles_tri <= 179.99)
            invalid_quad = np.isin(elem_ids, quad_ids) & (sum_angles_quad <= 359.99)
            
            quality[invalid_tri | invalid_quad] = 0
        
        return quality, angles if 'angles' in locals() else quality
    
    def diagonals(self, elem_ids=None, store='matrix'):
        """
        Get the diagonals of quadrilateral elements.
        
        Parameters
        ----------
        elem_ids : array_like, optional
            Element IDs to get diagonals for. If None, get for all elements.
        store : str, optional
            Output format: 'matrix' or 'cell' (list of lists)
            
        Returns
        -------
        vert_ids : ndarray or list
            Vertex pairs forming diagonals
        elem_ids : ndarray, optional
            Element IDs corresponding to each diagonal
        """
        if elem_ids is None:
            elem_ids = np.arange(self.n_elems)
        else:
            elem_ids = np.array(elem_ids)
        
        # Get triangles and quads
        tri_ids, quad_ids = self.elem_type(elem_ids)
        
        # Initialize output arrays
        diag1 = np.zeros((len(elem_ids), 2), dtype=np.int32)
        diag2 = np.zeros((len(elem_ids), 2), dtype=np.int32)
        
        # Process quads
        if len(quad_ids) > 0:
            quad_mask = np.isin(elem_ids, quad_ids)
            quad_indices = elem_ids[quad_mask]
            
            for i, elem_id in enumerate(quad_indices):
                conn = self.connectivity_list[elem_id]
                
                # First diagonal: vertices 0 and 2
                diag1[quad_mask][i] = [conn[0], conn[2]]
                
                # Second diagonal: vertices 1 and 3
                diag2[quad_mask][i] = [conn[1], conn[3]]
        
        # Format output
        if store.lower() == 'matrix':
            return np.array([diag1, diag2]).transpose(1, 0, 2)
        else:  # cell format
            result = []
            for i in range(len(elem_ids)):
                if elem_ids[i] in quad_ids:
                    result.append([diag1[i].tolist(), diag2[i].tolist()])
                else:  # Triangle - no diagonals
                    result.append([[0, 0], [0, 0]])
            
            if len(result) == 1:
                return result[0]
            return result

    def medians(self, edge_ids=None):
        """
        Get the median vertices for edges (vertices opposite to each edge).
        
        Parameters
        ----------
        edge_ids : array_like, optional
            Edge IDs to get medians for. If None, get for all edges.
            
        Returns
        -------
        vert_ids : list of lists
            Opposite vertices for each edge's neighboring elements
        """
        if edge_ids is None:
            edge_ids = np.arange(self.n_edges)
        else:
            edge_ids = np.array(edge_ids)
        
        result = []
        
        for edge_id in edge_ids:
            # Get vertices of the edge
            edge_verts = self.edge2vert([edge_id])[0]
            
            # Get elements attached to the edge
            edge_elems = self.edge2elem([edge_id])[0]
            edge_elems = edge_elems[edge_elems >= 0]  # Filter out -1 values
            
            # Find opposite vertices in each element
            opposite_verts = [[], []]
            
            for i, elem_id in enumerate(edge_elems):
                if i < 2:  # We only track up to 2 neighboring elements
                    # Get all vertices of the element
                    elem_verts = self.connectivity_list[elem_id]
                    
                    # Find vertices not in the edge
                    for v in elem_verts:
                        if v not in edge_verts and v > 0:  # Skip invalid vertices
                            opposite_verts[i].append(v)
            
            result.append(opposite_verts)
        
        return result

    def mesh_layers(self):
        """
        Discretize the mesh into layers starting from the boundary.
        
        Returns
        -------
        self : CHILmesh
            Updated object with layer information
        """
        # Get edges and their neighbors
        edge2vert = self.adjacencies['edge2vert']
        edge2elem = self.adjacencies['edge2elem']
        
        # Initialize layers
        self.layers = {
            'OE': [],  # Outer elements
            'IE': [],  # Inner elements
            'OV': [],  # Outer vertices
            'IV': [],  # Inner vertices
            'bEdgeIDs': []  # Boundary edge IDs
        }
        
        # Track vertices and elements
        vert_ids = np.arange(self.n_verts)
        elem_ids = np.arange(self.n_elems)
        
        layer = 0
        while np.any(edge2elem > -1):
            # Identify outer (boundary) edges of current layer
            if layer == 0:
                # First layer - use mesh boundary
                boundary_edges = self.boundary_edges()[0]
                bv = np.unique(edge2vert[boundary_edges].flatten())
                self.layers['OV'].append(bv)
            else:
                # Subsequent layers - use edges with only one neighbor
                boundary_edges = np.where(np.sum(edge2elem > -1, axis=1) == 1)[0]
                self.layers['OV'].append(np.unique(edge2vert[boundary_edges].flatten()))
            
            self.layers['bEdgeIDs'].append(boundary_edges)
            
            # Find elements connected to boundary edges (outer elements)
            outer_elems = np.unique(edge2elem[boundary_edges].flatten())
            outer_elems = outer_elems[outer_elems >= 0]
            self.layers['OE'].append(outer_elems)
            
            # Mark edges used by Ensure points are 3D (add z=0 if needed)
            if self.points.shape[1] == 2:
                self.points = np.column_stack((self.points, np.zeros(self.n_verts)))
            
            # Ensure counter-clockwise orientation
            self.is_poly_ccw()
            
            # Build mesh adjacencies
            self.build_adjacencies()
            
            # Identify mesh layers
            self.mesh_layers()
    
    def _process_inputs(self, args):
        """Process the input arguments to set up the mesh."""
        if len(args) == 0:
            # Interactive mode would go here, but skipping for now
            return
            
        elif len(args) == 1:
            # Single argument: either a filename or a triangulation
            if isinstance(args[0], str):
                # It's a filename
                self.read_grid_file(args[0])
            else:
                # Assume it's some kind of triangulation object
                # (would need adapting for specific Python mesh libraries)
                self.grid_name = "Imported Triangulation"
                if hasattr(args[0], 'points'):
                    self.points = np.array(args[0].points, dtype=float)
                if hasattr(args[0], 'connectivity_list') or hasattr(args[0], 'elements'):
                    conn_attr = 'connectivity_list' if hasattr(args[0], 'connectivity_list') else 'elements'
                    self.connectivity_list = np.array(getattr(args[0], conn_attr), dtype=np.int32)
                self.boundary_condition = {'id': []}
                
        elif len(args) >= 2:
            # Two or more arguments: connectivity_list, points, [bc], [grid_name]
            self.connectivity_list = np.array(args[0], dtype=np.int32)
            self.points = np.array(args[1], dtype=float)
            
            # Handle optional arguments
            if len(args) >= 3:
                grid_name_given = False
                for i, arg in enumerate(args[2:], start=2):
                    if isinstance(arg, str):
                        self.grid_name = arg
                        grid_name_given = True
                    elif i == 2 and not grid_name_given:
                        # Third argument and not a grid name, assume boundary condition
                        self.boundary_condition = arg
                
                if not grid_name_given:
                    self.grid_name = "Custom Input"
            else:
                self.grid_name = "Custom Input"
                self.boundary_condition = {'id': []}

    def read_grid_file(self, filename: str):
        """
        Read in a mesh grid file.
        
        Parameters
        ----------
        filename : str
            Path to the mesh file
            
        Returns
        -------
        self : CHILmesh
            Updated CHILmesh object
        """
        self.filename = filename
        _, grid_name, ext = os.path.basename(filename).rpartition('.')
        
        # Read file based on extension
        if ext == '14':
            self, bc = self.read_fort14(filename)
        elif ext == 'msh':
            self, bc = self.read_gmsh(filename)
        else:
            raise ValueError(f"Mesh file format '{ext}' not supported")
            
        self.grid_name = grid_name
        self.boundary_condition = bc
        return self

    def read_fort14(self, filename: str) -> Tuple['CHILmesh', Dict]:
        """
        Read in FE mesh from the ADCIRC text file format.
        
        Parameters
        ----------
        filename : str
            Path to the FORT14 file
            
        Returns
        -------
        mesh : CHILmesh
            Updated mesh object
        bc : dict
            Boundary condition information
        """
        with open(filename, 'r') as f:
            # Read grid name
            self.grid_name = f.readline().strip()
            
            # Read number of elements and points
            info_line = f.readline().split()
            n_elems = int(info_line[0])
            n_points = int(info_line[1])
            
            # Read points
            points = np.zeros((n_points, 3))
            for i in range(n_points):
                line = f.readline().split()
                # Skip the first column (index) and read x, y, z
                points[i, 0] = float(line[1])
                points[i, 1] = float(line[2])
                if len(line) > 3:
                    points[i, 2] = float(line[3])
            
            # Read connectivity list
            conn_list = np.zeros((n_elems, 5), dtype=np.int32)
            for i in range(n_elems):
                line = f.readline().split()
                # First column is element number, skip it
                n_nodes = int(line[1])
                for j in range(n_nodes):
                    if j+2 < len(line):
                        conn_list[i, j] = int(line[j+2])
            
            # Process connectivity list to get the right format
            # Determine max number of nodes per element
            max_nodes = max(np.count_nonzero(conn_list[i, 1:]) for i in range(n_elems)) + 1
            
            if max_nodes == 3:  # Triangle mesh
                self.connectivity_list = conn_list[:, 1:4]
                self.mesh_type = 'Triangular'
            elif max_nodes == 4:  # Quad or mixed mesh
                self.connectivity_list = conn_list[:, 1:5]
                # Check for triangles (will be updated in build_adjacencies)
                if np.any(self.connectivity_list[:, 3] == 0):
                    self.mesh_type = 'Mixed-Element'
                else:
                    self.mesh_type = 'Quadrangular'
            
            # Read boundary conditions
            n_open_ocean = int(f.readline().strip())
            
            # Skip total number of elevation specified boundary nodes
            f.readline()
            
            # Process boundary conditions
            boundary_condition = []
            
            # Open ocean segments
            for k in range(n_open_ocean):
                bc = {
                    'id': -1,
                    'des': 'Open Ocean',
                    'nodes': [],
                    'att': 0
                }
                
                # Read number of nodes in this segment
                n_nodes = int(f.readline().strip())
                
                # Read node IDs
                nodes = []
                for _ in range(n_nodes):
                    nodes.append(int(f.readline().strip()))
                
                bc['nodes'] = np.array(nodes)
                boundary_condition.append(bc)
            
            # Normal flow boundary segments
            n_normal_flow = int(f.readline().strip())
            
            # Skip the number of normal flow specified boundary nodes
            f.readline()
            
            # Process normal flow boundary segments
            for k in range(n_normal_flow):
                line = f.readline().split()
                n_nodes = int(line[0])
                bound_type = int(line[1])
                
                bc = {
                    'id': bound_type,
                    'nodes': [],
                    'att': 0
                }
                
                # Different handling based on boundary type
                if bound_type in [0, 2, 10, 12, 20, 22]:
                    bc['des'] = 'External Boundary'
                    nodes = []
                    for _ in range(n_nodes):
                        nodes.append(int(f.readline().strip()))
                    bc['nodes'] = np.array(nodes)
                
                elif bound_type == 6:
                    bc['des'] = 'Radiation Boundary'
                    nodes = []
                    for _ in range(n_nodes):
                        nodes.append(int(f.readline().strip()))
                    bc['nodes'] = np.array(nodes)
                
                elif bound_type in [1, 11, 21]:
                    bc['des'] = 'Internal Boundary'
                    nodes = []
                    for _ in range(n_nodes):
                        nodes.append(int(f.readline().strip()))
                    bc['nodes'] = np.array(nodes)
                
                elif bound_type in [3, 13, 23]:
                    bc['des'] = 'External Barrier'
                    nodes = []
                    attrs = []
                    for _ in range(n_nodes):
                        line = f.readline().split()
                        nodes.append(int(line[0]))
                        attrs.append([float(line[1]), float(line[2])])
                    bc['nodes'] = np.array(nodes)
                    bc['att'] = np.array(attrs)
                
                # Add more boundary condition types as needed...
                
                boundary_condition.append(bc)
            
            # If there are no boundary conditions
            if n_open_ocean == 0 and n_normal_flow == 0:
                boundary_condition = {'id': [], 'des': [], 'nodes': [], 'att': []}
            
            self.points = points
            
            # Check if there are 1D elements
            idx_1d = np.any(conn_list[:, 3:] == 0, axis=1)
            if np.any(idx_1d):
                # Remove 1D elements from connectivity list
                self.connectivity_list = self.connectivity_list[~idx_1d]
                
                # Add 1D elements to boundary conditions
                if isinstance(boundary_condition, list):
                    bc_1d = {
                        'id': 19,
                        'des': '1D Elements',
                        'nodes': conn_list[idx_1d, 1:3],
                        'att': []
                    }
                    boundary_condition.append(bc_1d)
            
            return self, boundary_condition

    def read_gmsh(self, filename: str) -> Tuple['CHILmesh', Dict]:
        """
        Read in FE mesh from the Gmsh ASCII text file format.
        
        Parameters
        ----------
        filename : str
            Path to the Gmsh file
            
        Returns
        -------
        mesh : CHILmesh
            Updated mesh object
        bc : dict
            Boundary condition information
        """
        # This is a placeholder for the Gmsh file reader
        # A complete implementation would parse the .msh file format
        with open(filename, 'r') as f:
            # Read through header to find version
            version_line = None
            for i, line in enumerate(f):
                if '$MeshFormat' in line:
                    version_line = next(f)
                    break
            
            if version_line is None:
                raise ValueError("Invalid Gmsh format: $MeshFormat section not found")
            
            # Parse version
            version = float(version_line.split()[0])
            
            # Different handling based on Gmsh version
            # This is a simplified implementation
            nodes = []
            elements = []
            physical_names = {}
            
            # Reset file pointer and read content
            f.seek(0)
            reading_section = None
            
            for line in f:
                line = line.strip()
                
                # Check for section markers
                if line.startswith('$'):
                    if line == '$Nodes':
                        reading_section = 'nodes'
                        # Skip the count line
                        next(f)
                        continue
                    elif line == '$Elements':
                        reading_section = 'elements'
                        # Skip the count line
                        next(f)
                        continue
                    elif line == '$PhysicalNames':
                        reading_section = 'physical_names'
                        # Skip the count line
                        next(f)
                        continue
                    elif line.startswith('$End'):
                        reading_section = None
                        continue
                
                # Read data based on current section
                if reading_section == 'nodes':
                    parts = line.split()
                    if len(parts) >= 4:
                        # node_id = int(parts[0])
                        x = float(parts[1])
                        y = float(parts[2])
                        z = float(parts[3])
                        nodes.append([x, y, z])
                
                elif reading_section == 'elements':
                    parts = line.split()
                    if len(parts) >= 7:  # Simplified for triangle/quad elements
                        # elem_id = int(parts[0])
                        elem_type = int(parts[1])
                        # Handle different element types
                        if elem_type == 2:  # Triangle
                            # Skip tags
                            n_tags = int(parts[2])
                            node_indices = [int(parts[3+n_tags+i]) for i in range(3)]
                            elements.append([0] + node_indices)  # Add placeholder for type
                        elif elem_type == 3:  # Quadrilateral
                            # Skip tags
                            n_tags = int(parts[2])
                            node_indices = [int(parts[3+n_tags+i]) for i in range(4)]
                            elements.append([1] + node_indices)  # Add placeholder for type
                
                elif reading_section == 'physical_names':
                    parts = line.split()
                    if len(parts) >= 3:
                        dim = int(parts[0])
                        tag = int(parts[1])
                        name = parts[2].strip('"')
                        physical_names[(dim, tag)] = name
        
        # Convert to numpy arrays
        points = np.array(nodes)
        
        # Process connectivity list
        if elements:
            elem_array = np.array(elements, dtype=np.int32)
            # Separate by element type
            tri_mask = elem_array[:, 0] == 0
            quad_mask = elem_array[:, 0] == 1
            
            if np.any(tri_mask) and np.any(quad_mask):
                # Mixed elements - create a quad-compatible list with triangles having repeated vertices
                conn_list = np.zeros((len(elements), 4), dtype=np.int32)
                conn_list[tri_mask, :3] = elem_array[tri_mask, 1:4]
                conn_list[tri_mask, 3] = elem_array[tri_mask, 1]  # Repeat first vertex
                conn_list[quad_mask] = elem_array[quad_mask, 1:5]
                self.mesh_type = 'Mixed-Element'
            elif np.any(tri_mask):
                conn_list = elem_array[tri_mask, 1:4]
                self.mesh_type = 'Triangular'
            elif np.any(quad_mask):
                conn_list = elem_array[quad_mask, 1:5]
                self.mesh_type = 'Quadrangular'
            else:
                conn_list = np.array([])
        else:
            conn_list = np.array([])
        
        self.points = points
        self.connectivity_list = conn_list
        
        # Create boundary conditions from physical names
        boundary_condition = physical_names
        
        return self, boundary_condition

    def build_adjacencies(self):
        """
        Build the 6 adjacency lists for the mesh.
        
        This method computes:
        1. Elem2Vert - Vertices of each element
        2. Edge2Vert - Vertices of each edge
        3. Elem2Edge - Edges of each element
        4. Vert2Edge - Edges attached to each vertex
        5. Vert2Elem - Elements attached to each vertex
        6. Edge2Elem - Elements attached to each edge
        
        Returns
        -------
        self : CHILmesh
            Updated CHILmesh object
        """
        # Check if input is triangular or quadrilateral mesh
        n_cols = self.connectivity_list.shape[1]
        
        if n_cols == 3:
            self.mesh_type = 'Triangular'
        elif n_cols == 4:
            # Check if it's mixed (triangles represented as quads)
            has_triangles = False
            for i in range(self.n_elems):
                if (self.connectivity_list[i, 0] == self.connectivity_list[i, 1] or
                    self.connectivity_list[i, 1] == self.connectivity_list[i, 2] or
                    self.connectivity_list[i, 2] == self.connectivity_list[i, 3] or
                    self.connectivity_list[i, 3] == self.connectivity_list[i, 0] or
                    self.connectivity_list[i, 3] == 0):
                    has_triangles = True
                    break
            
            if has_triangles:
                self.mesh_type = 'Mixed-Element'
                # For triangles in a quad mesh, make the last vertex same as first
                tri_mask = (self.connectivity_list[:, 0] == self.connectivity_list[:, 1]) | \
                           (self.connectivity_list[:, 1] == self.connectivity_list[:, 2]) | \
                           (self.connectivity_list[:, 2] == self.connectivity_list[:, 3]) | \
                           (self.connectivity_list[:, 3] == self.connectivity_list[:, 0]) | \
                           (self.connectivity_list[:, 3] == 0)
                self.connectivity_list[tri_mask, 3] = self.connectivity_list[tri_mask, 0]
            else:
                self.mesh_type = 'Quadrangular'
        
        # Ensure connectivity list is properly typed
        self.connectivity_list = self.connectivity_list.astype(np.int32)
        
        # Step 1: Identify edges of the mesh
        edges = []
        if n_cols == 3:  # Triangular mesh
            for i in range(self.n_elems):
                e1 = sorted([self.connectivity_list[i, 0], self.connectivity_list[i, 1]])
                e2 = sorted([self.connectivity_list[i, 1], self.connectivity_list[i, 2]])
                e3 = sorted([self.connectivity_list[i, 2], self.connectivity_list[i, 0]])
                edges.extend([e1, e2, e3])
        else:  # Quadrilateral or mixed mesh
            for i in range(self.n_elems):
                e1 = sorted([self.connectivity_list[i, 0], self.connectivity_list[i, 1]])
                e2 = sorted([self.connectivity_list[i, 1], self.connectivity_list[i, 2]])
                e3 = sorted([self.connectivity_list[i, 2], self.connectivity_list[i, 3]])
                e4 = sorted([self.connectivity_list[i, 3], self.connectivity_list[i, 0]])
                edges.extend([e1, e2, e3, e4])
        
        # Get unique edges and create a mapping
        unique_edges = []
        edge_to_id = {}
        for edge in edges:
            edge_tuple = tuple(edge)
            if edge_tuple not in edge_to_id:
                edge_to_id[edge_tuple] = len(unique_edges)
                unique_edges.append(edge)
        
        self.n_edges = len(unique_edges)
        
        # Create adjacency lists
        adjacencies = {
            'elem2vert': self.connectivity_list,
            'edge2vert': np.array(unique_edges, dtype=np.int32),
        }
        
        # Initialize other adjacency lists
        elem2edge = np.zeros((self.n_elems, n_cols), dtype=np.int32)
        
        # Build elem2edge
        if n_cols == 3:  # Triangular mesh
            for i in range(self.n_elems):
                e1 = tuple(sorted([self.connectivity_list[i, 0], self.connectivity_list[i, 1]]))
                e2 = tuple(sorted([self.connectivity_list[i, 1], self.connectivity_list[i, 2]]))
                e3 = tuple(sorted([self.connectivity_list[i, 2], self.connectivity_list[i, 0]]))
                elem2edge[i, 0] = edge_to_id[e1]
                elem2edge[i, 1] = edge_to_id[e2]
                elem2edge[i, 2] = edge_to_id[e3]
        else:  # Quadrilateral or mixed mesh
            for i in range(self.n_elems):
                e1 = tuple(sorted([self.connectivity_list[i, 0], self.connectivity_list[i, 1]]))
                e2 = tuple(sorted([self.connectivity_list[i, 1], self.connectivity_list[i, 2]]))
                e3 = tuple(sorted([self.connectivity_list[i, 2], self.connectivity_list[i, 3]]))
                e4 = tuple(sorted([self.connectivity_list[i, 3], self.connectivity_list[i, 0]]))
                elem2edge[i, 0] = edge_to_id[e1]
                elem2edge[i, 1] = edge_to_id[e2]
                elem2edge[i, 2] = edge_to_id[e3]
                elem2edge[i, 3] = edge_to_id[e4]
        
        adjacencies['elem2edge'] = elem2edge
        
        # Build edge2elem
        edge2elem = np.zeros((self.n_edges, 2), dtype=np.int32)
        for i in range(self.n_elems):
            for j in range(n_cols):
                edge_id = elem2edge[i, j]
                if edge2elem[edge_id, 0] == 0:
                    edge2elem[edge_id, 0] = i + 1  # 1-indexed to avoid confusion with 0
                else:
                    edge2elem[edge_id, 1] = i + 1
        
        # Convert back to 0-indexed
        edge2elem = edge2elem - 1
        edge2elem[edge2elem < 0] = -1  # Mark invalid entries
        
        adjacencies['edge2elem'] = edge2elem
        
        # Build vert2edge and vert2elem
        vert2edge = [[] for _ in range(self.n_verts)]
        vert2elem = [[] for _ in range(self.n_verts)]
        
        for i, edge in enumerate(unique_edges):
            v1, v2 = edge
            vert2edge[v1].append(i)
            vert2edge[v2].append(i)
        
        for i in range(self.n_elems):
            for v in self.connectivity_list[i]:
                if v < self.n_verts:  # Guard against invalid indices
                    vert2elem[v].append(i)
        
        # Convert lists to arrays
        max_edges_per_vert = max(len(v) for v in vert2edge)
        max_elems_per_vert = max(len(v) for v in vert2elem)
        
        vert2edge_array = np.zeros((self.n_verts, max_edges_per_vert), dtype=np.int32) - 1
        vert2elem_array = np.zeros((self.n_verts, max_elems_per_vert), dtype=np.int32) - 1
        
        for i, edges in enumerate(vert2edge):
            vert2edge_array[i, :len(edges)] = edges
        
        for i, elems in enumerate(vert2elem):
            vert2elem_array[i, :len(elems)] = elems
        
        adjacencies['vert2edge'] = vert2edge_array
        adjacencies['vert2elem'] = vert2elem_array
        
        self.adjacencies = adjacencies
        return self

    def is_poly_ccw(self, elem_ids=None, output_type='index'):
        """
        Check if polygons have counter-clockwise vertex orientation and fix if needed.
        
        Parameters
        ----------
        elem_ids : array_like, optional
            Element IDs to check. If None, check all elements.
        output_type : str, optional
            Type of output indices: 'index' or 'logical'
            
        Returns
        -------
        iCW : array
            Indices of elements with clockwise orientation
        """
        if elem_ids is None:
            elem_ids = np.arange(self.n_elems)
        else:
            elem_ids = np.array(elem_ids)
        
        # Compute signed areas
        areas = self.signed_area(elem_ids=elem_ids)
        
        # Identify clockwise elements
        cw_mask = areas < 0
        
        if output_type.lower() == 'logical':
            iCW = cw_mask
        else:
            iCW = elem_ids[cw_mask]
        
        # Fix clockwise elements
        if np.any(cw_mask):
            if self.mesh_type == 'Triangular':
                # For triangles, swap the last two vertices
                self.connectivity_list[iCW, 1:3] = self.connectivity_list[iCW, 2:0:-1]
            else:
                # For quads, reverse the order keeping the first vertex fixed
                for i in iCW:
                    self.connectivity_list[i, 1:] = self.connectivity_list[i, -1:0:-1]
        
        return iCW

    def signed_area(self, sign='signed', elem_ids=None):
        """
        Compute the signed area of triangular or quadrilateral elements.
        
        Parameters
        ----------
        sign : str, optional
            'signed' or 'unsigned' to return signed or absolute areas
        elem_ids : array_like, optional
            Element IDs to compute areas for. If None, compute for all elements.
            
        Returns
        -------
        areas : ndarray
            Signed or unsigned areas of the elements
        """
        if elem_ids is None:
            elem_ids = np.arange(self.n_elems)
        else:
            elem_ids = np.array(elem_ids)
        
        # Identify triangles and quads
        tri_ids, quad_ids = self.elem_type(elem_ids)
        
        # Initialize areas
        areas = np.zeros(len(elem_ids))
        
        # Process triangles
        if len(tri_ids) > 0:
            tri_mask = np.isin(elem_ids, tri_ids)
            tri_indices = elem_ids[tri_mask]
            
            # Get vertex coordinates
            v1 = self.connectivity_list[tri_indices, 0]
            v2 = self.connectivity_list[tri_indices, 1]
            v3 = self.connectivity_list[tri_indices, 2]
            
            x1, y1 = self.points[v1, 0], self.points[v1, 1]
            x2, y2 = self.points[v2, 0], self.points[v2, 1]
            x3, y3 = self.points[v3, 0], self.points[v3, 1]
            
            # Compute signed area using the shoelace formula
            areas[tri_mask] = 0.5 * (
                (x1 * y2 - x2 * y1) +
                (x2 * y3 - x3 * y2) +
                (x3 * y1 - x1 * y3)
            )
        
        # Process quads
        if len(quad_ids) > 0:
            quad_mask = np.isin(elem_ids, quad_ids)
            quad_indices = elem_ids[quad_mask]
            
            # Get vertex coordinates
            v1 = self.connectivity_list[quad_indices, 0]
            v2 = self.connectivity_list[quad_indices, 1]
            v3 = self.connectivity_list[quad_indices, 2]
            v4 = self.connectivity_list[quad_indices, 3]
            
            x1, y1 = self.points[v1, 0], self.points[v1, 1]
            x2, y2 = self.points[v2, 0], self.points[v2, 1]
            x3, y3 = self.points[v3, 0], self.points[v3, 1]
            x4, y4 = self.points[v4, 0], self.points[v4, 1]
            
            #
