//! Quality analysis and mesh queries for CHILmesh.
//!
//! This module implements:
//! - Signed area computation (shoelace formula for arbitrary polygons)
//! - CCW orientation enforcement
//! - Vertex/element query methods (get_vertex_edges, get_vertex_elements, get_element_vertices)

use ndarray::{Array1, Array2};
use std::collections::HashSet;

/// Compute signed areas of elements using the shoelace formula.
///
/// For a polygon with vertices (x0, y0), (x1, y1), ..., (xn, yn), the signed area is:
///   0.5 * sum_{i=0}^{n-1} (x_i * (y_{i+1} - y_{i-1}))
///
/// Positive area indicates counter-clockwise (CCW) orientation.
/// Negative area indicates clockwise (CW) orientation.
///
/// For triangles (padded as [v0, v1, v2, v2]): uses only first 3 vertices.
/// For quads [v0, v1, v2, v3]: uses all 4 vertices.
///
/// # Arguments
/// * `points` - Array2<f64> [n_verts, 2] — vertex coordinates (x, y)
/// * `connectivity` - Array2<i32> [n_elems, 3|4] — element-to-vertex indices
///
/// # Returns
/// Array1<f64> [n_elems] — signed area for each element
pub fn compute_signed_areas(
    points: &Array2<f64>,
    connectivity: &Array2<i32>,
) -> Array1<f64> {
    let n_elems = connectivity.shape()[0];
    let elem_cols = connectivity.shape()[1];
    let mut areas = Array1::<f64>::zeros(n_elems);

    if elem_cols == 3 {
        // Pure triangular mesh — vectorized computation
        for elem_idx in 0..n_elems {
            let v0 = connectivity[[elem_idx, 0]] as usize;
            let v1 = connectivity[[elem_idx, 1]] as usize;
            let v2 = connectivity[[elem_idx, 2]] as usize;

            let x0 = points[[v0, 0]];
            let y0 = points[[v0, 1]];
            let x1 = points[[v1, 0]];
            let y1 = points[[v1, 1]];
            let x2 = points[[v2, 0]];
            let y2 = points[[v2, 1]];

            // Shoelace formula for triangle
            areas[elem_idx] = 0.5
                * (x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1));
        }
        return areas;
    }

    // 4-column connectivity: separate padded triangles from quads
    for elem_idx in 0..n_elems {
        let v0 = connectivity[[elem_idx, 0]];
        let v1 = connectivity[[elem_idx, 1]];
        let v2 = connectivity[[elem_idx, 2]];
        let v3 = connectivity[[elem_idx, 3]];

        // Check if this is a padded triangle (one vertex appears twice)
        let is_triangle = (v0 == v1)
            || (v1 == v2)
            || (v2 == v3)
            || (v3 == v0)
            || (v0 == v2)
            || (v1 == v3);

        if is_triangle {
            // Use only first 3 vertices
            let v0_idx = v0 as usize;
            let v1_idx = v1 as usize;
            let v2_idx = v2 as usize;

            let x0 = points[[v0_idx, 0]];
            let y0 = points[[v0_idx, 1]];
            let x1 = points[[v1_idx, 0]];
            let y1 = points[[v1_idx, 1]];
            let x2 = points[[v2_idx, 0]];
            let y2 = points[[v2_idx, 1]];

            areas[elem_idx] =
                0.5 * (x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1));
        } else {
            // Quadrilateral: use all 4 vertices
            let v0_idx = v0 as usize;
            let v1_idx = v1 as usize;
            let v2_idx = v2 as usize;
            let v3_idx = v3 as usize;

            let x0 = points[[v0_idx, 0]];
            let y0 = points[[v0_idx, 1]];
            let x1 = points[[v1_idx, 0]];
            let y1 = points[[v1_idx, 1]];
            let x2 = points[[v2_idx, 0]];
            let y2 = points[[v2_idx, 1]];
            let x3 = points[[v3_idx, 0]];
            let y3 = points[[v3_idx, 1]];

            // Shoelace formula for quadrilateral
            areas[elem_idx] = 0.5
                * (x0 * (y1 - y3) + x1 * (y2 - y0) + x2 * (y3 - y1) + x3 * (y0 - y2));
        }
    }

    areas
}

/// Ensure counter-clockwise (CCW) orientation for all elements.
///
/// Reorders vertices of any element with negative signed area (CW orientation)
/// to flip the area sign. For mixed-element meshes:
/// - Triangles [v0, v1, v2, v2]: reverse to [v0, v2, v1, v2]
/// - Quads [v0, v1, v2, v3]: reverse to [v0, v3, v2, v1]
///
/// # Arguments
/// * `points` - Array2<f64> [n_verts, 2] — vertex coordinates (read-only)
/// * `connectivity` - Array2<i32> [n_elems, 3|4] — element-to-vertex (MUTABLE)
///
/// After execution, all elements will have positive signed area.
pub fn ensure_ccw_orientation(
    points: &Array2<f64>,
    connectivity: &mut Array2<i32>,
) {
    let areas = compute_signed_areas(points, connectivity);
    let n_elems = connectivity.shape()[0];
    let elem_cols = connectivity.shape()[1];

    // Find all elements with negative (CW) area
    for elem_idx in 0..n_elems {
        if areas[elem_idx] < 0.0 {
            // Reverse vertex order to flip orientation
            if elem_cols == 3 {
                // Triangle: [v0, v1, v2] → [v0, v2, v1]
                let temp = connectivity[[elem_idx, 1]];
                connectivity[[elem_idx, 1]] = connectivity[[elem_idx, 2]];
                connectivity[[elem_idx, 2]] = temp;
            } else {
                // Quad: [v0, v1, v2, v3] → [v0, v3, v2, v1]
                let temp1 = connectivity[[elem_idx, 1]];
                let temp3 = connectivity[[elem_idx, 3]];
                connectivity[[elem_idx, 1]] = temp3;
                connectivity[[elem_idx, 3]] = temp1;
            }
        }
    }
}

/// Get all edges incident to a vertex.
///
/// For each edge in the edges array, check if it contains the given vertex.
/// Return deduplicated list of incident edge IDs in sorted order.
///
/// # Arguments
/// * `v` - Vertex ID
/// * `edges` - Array2<i32> [n_edges, 2] — edge connectivity
///
/// # Returns
/// Vec<usize> — sorted list of incident edge IDs
pub fn get_vertex_edges(v: usize, edges: &Array2<i32>) -> Vec<usize> {
    let n_edges = edges.shape()[0];
    let mut incident_edges = HashSet::new();

    for edge_idx in 0..n_edges {
        let v0 = edges[[edge_idx, 0]] as usize;
        let v1 = edges[[edge_idx, 1]] as usize;

        if v0 == v || v1 == v {
            incident_edges.insert(edge_idx);
        }
    }

    let mut result: Vec<usize> = incident_edges.into_iter().collect();
    result.sort_unstable();
    result
}

/// Get all elements incident to a vertex.
///
/// For each element in the connectivity array, check if it contains the given vertex.
/// Return deduplicated list of incident element IDs in sorted order.
///
/// # Arguments
/// * `v` - Vertex ID
/// * `connectivity` - Array2<i32> [n_elems, 3|4] — element-to-vertex
///
/// # Returns
/// Vec<usize> — sorted list of incident element IDs
pub fn get_vertex_elements(v: usize, connectivity: &Array2<i32>) -> Vec<usize> {
    let n_elems = connectivity.shape()[0];
    let elem_cols = connectivity.shape()[1];
    let mut incident_elems = HashSet::new();

    for elem_idx in 0..n_elems {
        for col_idx in 0..elem_cols {
            if connectivity[[elem_idx, col_idx]] as usize == v {
                incident_elems.insert(elem_idx);
                break;
            }
        }
    }

    let mut result: Vec<usize> = incident_elems.into_iter().collect();
    result.sort_unstable();
    result
}

/// Get all vertices of an element as a 4-element array.
///
/// Extracts the element row from connectivity. For triangles (padded as [v0, v1, v2, v2]),
/// the padding is preserved in the returned array.
///
/// # Arguments
/// * `e` - Element ID
/// * `connectivity` - Array2<i32> [n_elems, 3|4] — element-to-vertex
///
/// # Returns
/// Array1<i32> — 4-element array [v0, v1, v2, v3] (triangles have v3 = v2)
pub fn get_element_vertices(e: usize, connectivity: &Array2<i32>) -> Array1<i32> {
    let elem_cols = connectivity.shape()[1];

    if elem_cols == 3 {
        // Triangle: expand to 4 elements by duplicating the last vertex
        let mut result = Array1::<i32>::zeros(4);
        result[0] = connectivity[[e, 0]];
        result[1] = connectivity[[e, 1]];
        result[2] = connectivity[[e, 2]];
        result[3] = connectivity[[e, 2]]; // Padding
        result
    } else {
        // Quad: already 4 elements
        let mut result = Array1::<i32>::zeros(4);
        result[0] = connectivity[[e, 0]];
        result[1] = connectivity[[e, 1]];
        result[2] = connectivity[[e, 2]];
        result[3] = connectivity[[e, 3]];
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_signed_area_triangle_ccw() {
        let points = Array2::from_shape_vec(
            (3, 2),
            vec![0.0, 0.0, 1.0, 0.0, 0.5, 1.0],
        )
        .unwrap();
        let connectivity = Array2::from_shape_vec((1, 3), vec![0, 1, 2]).unwrap();
        let areas = compute_signed_areas(&points, &connectivity);
        // CCW triangle: area should be positive
        assert!(areas[0] > 0.0);
        // Expected area: 0.5 * 1 * 1 = 0.5
        assert!((areas[0] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_signed_area_triangle_cw() {
        let points = Array2::from_shape_vec(
            (3, 2),
            vec![0.0, 0.0, 0.5, 1.0, 1.0, 0.0],
        )
        .unwrap();
        let connectivity = Array2::from_shape_vec((1, 3), vec![0, 1, 2]).unwrap();
        let areas = compute_signed_areas(&points, &connectivity);
        // CW triangle: area should be negative
        assert!(areas[0] < 0.0);
        // Expected area: -0.5
        assert!((areas[0] + 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_get_vertex_edges() {
        // Two edges: (0,1) and (1,2)
        let edges = Array2::from_shape_vec((2, 2), vec![0, 1, 1, 2]).unwrap();
        let incident = get_vertex_edges(1, &edges);
        assert_eq!(incident, vec![0, 1]);
    }

    #[test]
    fn test_get_vertex_elements() {
        // Two triangles: [0,1,2] and [1,2,3]
        let connectivity = Array2::from_shape_vec((2, 3), vec![0, 1, 2, 1, 2, 3]).unwrap();
        let incident = get_vertex_elements(1, &connectivity);
        assert_eq!(incident, vec![0, 1]);
    }

    #[test]
    fn test_get_element_vertices() {
        let connectivity = Array2::from_shape_vec((1, 3), vec![0, 1, 2]).unwrap();
        let verts = get_element_vertices(0, &connectivity);
        assert_eq!(verts[0], 0);
        assert_eq!(verts[1], 1);
        assert_eq!(verts[2], 2);
        assert_eq!(verts[3], 2); // Padding
    }
}
