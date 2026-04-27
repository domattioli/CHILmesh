"""Bridge adapters for downstream projects.

This module provides convenient interfaces for integrating CHILmesh
with MADMESHR, ADMESH, and ADMESH-Domains. Each adapter adds
domain-specific convenience methods while delegating to the base
CHILmesh public API (CAI).
"""

from typing import Dict, Set, Tuple, List, Optional
import numpy as np


class MeshAdapterForMADMESHR:
    """
    Adapter for MADMESHR mesh adaptation research.

    Provides convenience methods for mesh refinement patterns
    specific to MADMESHR workflows, particularly element
    neighbor queries and local quality assessment.

    Example:
        >>> from chilmesh import CHILmesh
        >>> from chilmesh.bridge import MeshAdapterForMADMESHR
        >>> mesh = CHILmesh.read_from_fort14("domain.14")
        >>> adapter = MeshAdapterForMADMESHR(mesh)
        >>> neighbors = adapter.get_element_neighbors(elem_id)
    """

    def __init__(self, mesh):
        """
        Initialize adapter with base CHILmesh.

        Parameters:
            mesh: CHILmesh instance
        """
        self.mesh = mesh

    def get_element_neighbors(self, elem_id: int) -> Set[int]:
        """
        Get all elements adjacent to the given element.

        Returns all elements that share at least one edge with the
        given element, not including the element itself.

        Parameters:
            elem_id: Element index [0, n_elems)

        Returns:
            Set of adjacent element IDs (may be empty for isolated elements)

        Raises:
            ValueError: If elem_id out of range

        Complexity:
            Time: O(k) where k = element degree (typically 3-4)
            Space: O(k) for returned set

        Example:
            >>> neighbors = adapter.get_element_neighbors(0)
            >>> for elem in neighbors:
            ...     print(f"Adjacent element: {elem}")
        """
        if not 0 <= elem_id < self.mesh.n_elems:
            raise ValueError(f"Element {elem_id} out of range [0, {self.mesh.n_elems})")

        neighbors = set()
        edges = self.mesh.elem2edge(elem_id)

        for edge_id in edges:
            edge_elems = self.mesh.edge2elem(edge_id)[0]  # Get first (only) row
            e1, e2 = edge_elems
            if e1 == elem_id and e2 >= 0:
                neighbors.add(e2)
            elif e2 == elem_id and e1 >= 0:
                neighbors.add(e1)

        return neighbors

    def get_element_quality_neighborhood(self, elem_id: int) -> Dict:
        """
        Get quality metrics for element and all neighbors.

        Useful for assessing local mesh fitness before refinement decisions.

        Parameters:
            elem_id: Element index [0, n_elems)

        Returns:
            Dict with keys:
            - 'center': float - Quality of central element
            - 'neighbors': ndarray[k] - Quality of neighboring elements
            - 'mean': float - Mean quality in neighborhood
            - 'stats': Dict - Full statistics (min, max, std from elem_quality)

        Raises:
            ValueError: If elem_id out of range

        Complexity:
            Time: O(k) where k = neighborhood size (typically 4-8 elements)
            Space: O(k) for quality array

        Example:
            >>> metrics = adapter.get_element_quality_neighborhood(0)
            >>> print(f"Center quality: {metrics['center']:.3f}")
            >>> print(f"Neighborhood mean: {metrics['mean']:.3f}")
        """
        if not 0 <= elem_id < self.mesh.n_elems:
            raise ValueError(f"Element {elem_id} out of range [0, {self.mesh.n_elems})")

        neighbors = self.get_element_neighbors(elem_id)
        all_elems = [elem_id] + sorted(list(neighbors))

        quality, _, stats = self.mesh.elem_quality(all_elems)

        return {
            "center": quality[0],
            "neighbors": quality[1:] if len(quality) > 1 else np.array([]),
            "mean": float(np.mean(quality)),
            "stats": stats,
        }

    def get_refinement_region(
        self, elem_ids: List[int], include_neighbors: bool = True
    ) -> Set[int]:
        """
        Get refinement region from seed elements.

        Useful for mesh adaptation where you need to refine elements
        and their neighborhoods.

        Parameters:
            elem_ids: List of seed element IDs to refine
            include_neighbors: If True, include neighbors of seed elements

        Returns:
            Set of element IDs in refinement region

        Complexity:
            Time: O(m * k) where m = seed count, k = avg element degree
            Space: O(region_size)

        Example:
            >>> poor_elements = [0, 5, 10]  # Example poor quality elements
            >>> region = adapter.get_refinement_region(poor_elements)
            >>> print(f"Refine {len(region)} elements")
        """
        region = set(elem_ids)

        if include_neighbors:
            for elem_id in elem_ids:
                neighbors = self.get_element_neighbors(elem_id)
                region.update(neighbors)

        return region


class MeshAdapterForADMESH:
    """
    Adapter for ADMESH mesh adaptation framework.

    Provides convenience methods for mesh quality assessment and
    adaptation criteria specific to ADMESH workflows.

    Example:
        >>> from chilmesh import CHILmesh
        >>> from chilmesh.bridge import MeshAdapterForADMESH
        >>> mesh = CHILmesh.read_from_fort14("domain.14")
        >>> adapter = MeshAdapterForADMESH(mesh)
        >>> report = adapter.get_mesh_quality_report()
    """

    def __init__(self, mesh):
        """
        Initialize adapter with base CHILmesh.

        Parameters:
            mesh: CHILmesh instance
        """
        self.mesh = mesh

    def get_mesh_quality_report(self, quality_type: str = "skew") -> Dict:
        """
        Get comprehensive mesh quality report.

        Useful for overall mesh assessment before adaptation decisions.

        Parameters:
            quality_type: Quality metric type ('skew' recommended)

        Returns:
            Dict with keys:
            - 'mean': float - Mean quality [0, 1]
            - 'min': float - Minimum quality
            - 'max': float - Maximum quality
            - 'std': float - Standard deviation
            - 'poor_count': int - Number of elements below threshold
            - 'poor_fraction': float - Fraction below threshold
            - 'poor_elements': ndarray - Indices of poor elements

        Complexity:
            Time: O(n_elems)
            Space: O(n_elems) for quality array

        Example:
            >>> report = adapter.get_mesh_quality_report()
            >>> print(f"Mean quality: {report['mean']:.3f}")
            >>> print(f"Poor elements: {report['poor_count']}")
        """
        quality, _, stats = self.mesh.elem_quality(quality_type=quality_type)

        threshold = 0.3
        poor_mask = quality < threshold
        poor_elements = np.where(poor_mask)[0]

        return {
            "mean": float(stats["mean"]),
            "min": float(stats["min"]),
            "max": float(stats["max"]),
            "std": float(stats.get("std", 0.0)),
            "poor_count": int(np.sum(poor_mask)),
            "poor_fraction": float(np.mean(poor_mask)),
            "poor_elements": poor_elements,
        }

    def get_element_angles_summary(self, elem_ids: Optional[List[int]] = None) -> Dict:
        """
        Get angle-based quality summary.

        Useful for angle-based adaptation criteria.

        Parameters:
            elem_ids: List of element IDs, or None for all elements

        Returns:
            Dict with keys:
            - 'min_angle_degrees': float - Global minimum angle
            - 'max_angle_degrees': float - Global maximum angle
            - 'elements_with_acute': int - Count of elements with angles < 30°
            - 'elements_with_obtuse': int - Count of elements with angles > 120°

        Complexity:
            Time: O(n_elems * avg_vertices_per_element)
            Space: O(n_elems)

        Example:
            >>> summary = adapter.get_element_angles_summary()
            >>> print(f"Min angle: {summary['min_angle_degrees']:.1f}°")
        """
        angles = self.mesh.interior_angles(elem_ids)

        # Convert to degrees for reporting
        angles_deg = np.degrees(angles)

        min_angle = np.min(angles_deg)
        max_angle = np.max(angles_deg)

        # Count problematic angles
        acute_count = int(np.sum(angles_deg < 30))
        obtuse_count = int(np.sum(angles_deg > 120))

        return {
            "min_angle_degrees": float(min_angle),
            "max_angle_degrees": float(max_angle),
            "elements_with_acute": acute_count,
            "elements_with_obtuse": obtuse_count,
        }


class MeshAdapterForADMESHDomains:
    """
    Adapter for ADMESH-Domains multi-domain handling.

    Provides convenience methods for domain-level queries and
    boundary extraction specific to ADMESH-Domains.

    Example:
        >>> from chilmesh import CHILmesh
        >>> from chilmesh.bridge import MeshAdapterForADMESHDomains
        >>> mesh = CHILmesh.read_from_fort14("domain.14")
        >>> adapter = MeshAdapterForADMESHDomains(mesh)
        >>> boundaries = adapter.get_domain_boundaries()
    """

    def __init__(self, mesh):
        """
        Initialize adapter with base CHILmesh.

        Parameters:
            mesh: CHILmesh instance
        """
        self.mesh = mesh

    def get_domain_boundaries(self) -> Dict[int, Set[int]]:
        """
        Get mesh boundary as a domain boundary.

        Returns boundary vertices grouped by connected components.
        Useful for multi-domain decomposition.

        Returns:
            Dict mapping domain_id to Set of boundary vertex IDs
            For single-domain mesh, returns {0: boundary_vertices}

        Complexity:
            Time: O(n_boundary_edges)
            Space: O(n_boundary_vertices)

        Example:
            >>> boundaries = adapter.get_domain_boundaries()
            >>> for domain_id, verts in boundaries.items():
            ...     print(f"Domain {domain_id}: {len(verts)} boundary vertices")
        """
        boundary_edge_ids = self.mesh.boundary_edges()

        boundary_verts = set()
        for edge_id in boundary_edge_ids:
            edge_verts = self.mesh.edge2vert(edge_id)[0]  # Get first (only) row
            v1, v2 = edge_verts
            boundary_verts.add(v1)
            boundary_verts.add(v2)

        # For now, treat entire mesh as single domain
        return {0: boundary_verts}

    def get_mesh_connectivity_info(self) -> Dict:
        """
        Get high-level connectivity information.

        Useful for understanding mesh structure for domain splitting.

        Returns:
            Dict with keys:
            - 'n_verts': int
            - 'n_elems': int
            - 'n_edges': int
            - 'avg_vertex_degree': float - Average edges per vertex
            - 'avg_element_degree': float - Average neighbors per element
            - 'boundary_edge_count': int
            - 'interior_edge_count': int

        Complexity:
            Time: O(n_verts + n_edges)
            Space: O(1)

        Example:
            >>> info = adapter.get_mesh_connectivity_info()
            >>> print(f"Average vertex degree: {info['avg_vertex_degree']:.1f}")
        """
        boundary_edges = len(self.mesh.boundary_edges())
        interior_edges = self.mesh.n_edges - boundary_edges

        # Estimate average vertex degree
        total_vertex_degree = 0
        for v in range(self.mesh.n_verts):
            edges = self.mesh.get_vertex_edges(v)
            total_vertex_degree += len(edges)

        avg_vertex_degree = (
            total_vertex_degree / self.mesh.n_verts
            if self.mesh.n_verts > 0
            else 0.0
        )

        # Estimate average element degree (neighbors)
        total_neighbors = 0
        for e in range(min(self.mesh.n_elems, 1000)):  # Sample to avoid O(n²)
            edges = self.mesh.elem2edge(e)
            neighbors = 0
            for edge_id in edges:
                edge_elems = self.mesh.edge2elem(edge_id)[0]  # Get first (only) row
                e1, e2 = edge_elems
                if e1 == e and e2 >= 0:
                    neighbors += 1
                elif e2 == e and e1 >= 0:
                    neighbors += 1
            total_neighbors += neighbors

        sample_count = min(self.mesh.n_elems, 1000)
        avg_element_degree = (
            total_neighbors / sample_count if sample_count > 0 else 0.0
        )

        return {
            "n_verts": self.mesh.n_verts,
            "n_elems": self.mesh.n_elems,
            "n_edges": self.mesh.n_edges,
            "avg_vertex_degree": float(avg_vertex_degree),
            "avg_element_degree": float(avg_element_degree),
            "boundary_edge_count": int(boundary_edges),
            "interior_edge_count": int(interior_edges),
        }
