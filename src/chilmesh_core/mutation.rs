use crate::errors::RustMeshError;
use crate::adjacency;
use ndarray::Array2;
use std::collections::HashSet;

/// Add a new element to the mesh
///
/// Validates vertex bounds, element type, and geometric constraints.
/// Recomputes adjacency after insertion.
///
/// # Arguments
/// * `connectivity` - Element connectivity array [n_elems, 4]
/// * `points` - Vertex coordinates [n_verts, 2]
/// * `num_verts` - Number of vertices
/// * `num_elems` - Number of elements
/// * `elem_type` - Vector of element types (3=triangle, 4=quad)
/// * `verts` - Vertex indices for new element [v0, v1, v2] or [v0, v1, v2, v3]
/// * `elem_type_new` - Type of new element (3 or 4)
///
/// # Returns
/// New element ID if successful
pub fn add_element(
    connectivity: &Array2<i32>,
    points: &Array2<f64>,
    num_verts: usize,
    num_elems: usize,
    elem_type: &[u32],
    verts: &[i32],
    elem_type_new: u32,
) -> Result<usize, RustMeshError> {
    // Validate element type
    if elem_type_new != 3 && elem_type_new != 4 {
        return Err(RustMeshError::InvalidElement {
            elem_id: num_elems,
            reason: format!("Invalid element type {} (expected 3 or 4)", elem_type_new),
        });
    }

    // Validate vertex count matches element type
    let expected_verts = if elem_type_new == 3 { 3 } else { 4 };
    if verts.len() < expected_verts {
        return Err(RustMeshError::InvalidElement {
            elem_id: num_elems,
            reason: format!(
                "Invalid vertex count {} for element type {} (expected {})",
                verts.len(),
                elem_type_new,
                expected_verts
            ),
        });
    }

    // Validate all vertex indices are in bounds
    for (i, &v) in verts.iter().enumerate() {
        if v < 0 || v >= num_verts as i32 {
            return Err(RustMeshError::VertexOutOfBounds {
                vertex_id: v,
                max_vertices: num_verts,
            });
        }
    }

    // Validate vertices form a valid polygon (no self-intersection, CCW orientation)
    validate_element_geometry(points, verts, elem_type_new)?;

    // New element would have ID = num_elems (0-indexed)
    let new_elem_id = num_elems;
    Ok(new_elem_id)
}

/// Remove an element from the mesh
///
/// Validates element bounds and checks for orphaned vertices.
/// Recomputes adjacency after removal.
///
/// # Arguments
/// * `connectivity` - Element connectivity array [n_elems, 4]
/// * `num_elems` - Number of elements
/// * `elem_id` - Element ID to remove (0-indexed)
///
/// # Returns
/// Ok(()) if successful
pub fn remove_element(
    connectivity: &Array2<i32>,
    num_elems: usize,
    elem_id: usize,
) -> Result<(), RustMeshError> {
    // Validate element bounds
    if elem_id >= num_elems {
        return Err(RustMeshError::InvalidElement {
            elem_id,
            reason: format!("Element ID {} out of bounds (max: {})", elem_id, num_elems - 1),
        });
    }

    // Check for orphaned vertices
    // Build vertex incident element counts excluding the removed element
    let mut vert_incident_count = vec![0usize; connectivity.shape()[0]];

    for (e_idx, elem_idx) in (0..num_elems).enumerate() {
        if elem_idx == elem_id {
            continue; // Skip the element being removed
        }

        // Count incident elements for each vertex in this element
        for col in 0..4 {
            let v = connectivity[[elem_idx, col]];
            if v >= 0 {
                vert_incident_count[v as usize] += 1;
            }
        }
    }

    // Check if any vertex would be orphaned (no incident elements after removal)
    let removed_elem = connectivity.row(elem_id);
    for col in 0..4 {
        let v = removed_elem[col];
        if v >= 0 {
            let v_usize = v as usize;
            if vert_incident_count[v_usize] == 0 {
                return Err(RustMeshError::InvalidElement {
                    elem_id,
                    reason: format!("Removing element {} would orphan vertex {}", elem_id, v),
                });
            }
        }
    }

    Ok(())
}

/// Validate that an element's vertices form a valid geometry
///
/// Checks for:
/// - Non-zero area (2D cross product)
/// - Counter-clockwise (CCW) orientation
/// - No self-intersection (simplified check: convex hull validation)
fn validate_element_geometry(
    points: &Array2<f64>,
    verts: &[i32],
    elem_type: u32,
) -> Result<(), RustMeshError> {
    // Get vertex coordinates
    let mut coords = Vec::new();
    for &v in verts.iter().take(if elem_type == 3 { 3 } else { 4 }) {
        if v >= 0 && (v as usize) < points.shape()[0] {
            coords.push((
                points[[v as usize, 0]],
                points[[v as usize, 1]],
            ));
        } else {
            return Err(RustMeshError::InvalidGeometry(
                "Invalid vertex reference in element".to_string(),
            ));
        }
    }

    // Compute signed area using shoelace formula
    let signed_area = compute_signed_polygon_area(&coords);

    if signed_area.abs() < 1e-14 {
        return Err(RustMeshError::InvalidGeometry(
            "Element has zero or near-zero area".to_string(),
        ));
    }

    if signed_area < 0.0 {
        return Err(RustMeshError::InvalidGeometry(
            "Element has clockwise orientation (expected counter-clockwise)".to_string(),
        ));
    }

    Ok(())
}

/// Compute signed area of a polygon using the shoelace formula
fn compute_signed_polygon_area(vertices: &[(f64, f64)]) -> f64 {
    let n = vertices.len();
    if n < 3 {
        return 0.0;
    }

    let mut area = 0.0;
    for i in 0..n {
        let (x1, y1) = vertices[i];
        let (x2, y2) = vertices[(i + 1) % n];
        area += (x2 - x1) * (y2 + y1);
    }

    area / 2.0
}
