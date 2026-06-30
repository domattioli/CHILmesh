//! Quad-edge topology backend for CHILmesh (Rust port).
//!
//! This module implements quad-edge construction and adjacency converters,
//! ported from Python reference implementation in mesh_topology_quadegg.py.
//!
//! The quad-edge structure stores directed edges with 4 components:
//!     edges: ndarray[n_edges, 4] with columns:
//!     - [0] origin_vertex: vertex ID at start of directed edge
//!     - [1] next_cw: next edge clockwise around incident face
//!     - [2] next_ccw: next edge counter-clockwise around opposite face
//!     - [3] opposite_idx: index of opposite directed edge (-1 for boundary)
//!
//! Sentinel conventions:
//!     - Boundary edges: opposite_idx == -1, next_ccw == -1
//!     - No per-object Python overhead; all lookups are array indexing

use ndarray::Array2;
use std::collections::HashMap;

/// Construct quad-edge topology from element-to-vertex connectivity.
///
/// O(n) algorithm:
/// 1. Phase 1: Create all directed edges per element in iteration order
/// 2. Phase 2: Hash-table lookup to pair opposite edges and assign next pointers
///
/// Returns edges: Array2<i32> [n_edges, 4] with (origin, next_cw, next_ccw, opposite_idx)
pub fn build_quadegg_from_connectivity(
    elem2vert: &Array2<i32>,
    _n_verts: usize,
) -> Array2<i32> {
    let n_elems = elem2vert.shape()[0];
    let elem_cols = elem2vert.shape()[1];

    // Determine element types (3 for triangle, 4 for quad)
    let mut elem_type = vec![0i32; n_elems];
    for i in 0..n_elems {
        if elem_cols == 3 {
            elem_type[i] = 3;
        } else {
            // If elem_cols >= 4, check if last vertex equals first (padding)
            elem_type[i] = if elem2vert[[i, 3]] != elem2vert[[i, 0]] {
                4
            } else {
                3
            };
        }
    }

    // Phase 1: Create all directed edges in element order
    let mut edges_list: Vec<Vec<i32>> = Vec::new();
    let mut edge_map: HashMap<(i32, i32), usize> = HashMap::new();
    let mut elem_to_edge_indices: Vec<Vec<usize>> = vec![Vec::new(); n_elems];

    let mut edge_idx_counter = 0usize;

    for elem_idx in 0..n_elems {
        let n_verts_in_elem = elem_type[elem_idx] as usize;

        for i in 0..n_verts_in_elem {
            let v_origin = elem2vert[[elem_idx, i]];
            let v_dest = elem2vert[[elem_idx, (i + 1) % n_verts_in_elem]];

            // Create edge: [origin, next_cw, next_ccw, opposite_idx]
            // next_cw will be assigned in same element loop below
            // next_ccw and opposite_idx assigned during Phase 2
            let edge = vec![v_origin, -1, -1, -1];
            edges_list.push(edge);
            elem_to_edge_indices[elem_idx].push(edge_idx_counter);

            // Record directed edge for opposite-pairing lookup
            let directed_edge = (v_origin, v_dest);
            edge_map.insert(directed_edge, edge_idx_counter);

            edge_idx_counter += 1;
        }
    }

    let mut edges_array = Array2::zeros((edges_list.len(), 4));
    for (i, edge) in edges_list.iter().enumerate() {
        for (j, &val) in edge.iter().enumerate() {
            edges_array[[i, j]] = val;
        }
    }

    // Assign next_cw within each element's edges
    for elem_idx in 0..n_elems {
        let edge_indices = &elem_to_edge_indices[elem_idx];
        let n_edges_in_elem = edge_indices.len();
        for i in 0..n_edges_in_elem {
            let curr_edge_idx = edge_indices[i];
            let next_edge_idx = edge_indices[(i + 1) % n_edges_in_elem];
            edges_array[[curr_edge_idx, 1]] = next_edge_idx as i32;
        }
    }

    // Phase 2: Pair opposite edges and assign next_ccw
    for ((v_origin, v_dest), edge_idx) in edge_map.iter() {
        let reverse_edge = (*v_dest, *v_origin);
        if let Some(&opposite_edge_idx) = edge_map.get(&reverse_edge) {
            // Found opposite edge
            edges_array[[*edge_idx, 3]] = opposite_edge_idx as i32;
            edges_array[[opposite_edge_idx, 3]] = *edge_idx as i32;

            // Assign next_ccw pointers (going around opposite face)
            // next_ccw of current edge = next_cw of opposite edge
            let next_cw_of_opposite = edges_array[[opposite_edge_idx, 1]];
            edges_array[[*edge_idx, 2]] = next_cw_of_opposite;

            // Symmetrically: next_ccw of opposite edge = next_cw of current edge
            let next_cw_of_current = edges_array[[*edge_idx, 1]];
            edges_array[[opposite_edge_idx, 2]] = next_cw_of_current;
        } else {
            // Boundary edge: opposite_idx = -1, next_ccw = -1
            edges_array[[*edge_idx, 3]] = -1;
            edges_array[[*edge_idx, 2]] = -1;
        }
    }

    edges_array
}

/// Extract Edge2Vert: unique undirected edges in first-encounter order.
///
/// Returns unique undirected edges in canonical form (min_vert, max_vert),
/// in first-encounter element-traversal order, matching Python/C++ edge IDs.
pub fn to_edge2vert(elem2vert: &Array2<i32>) -> Array2<i32> {
    let edges = first_encounter_canonical_edges(elem2vert);
    let n_edges = edges.len();
    let mut result_data = Vec::with_capacity(n_edges * 2);
    for (v0, v1) in &edges {
        result_data.push(*v0);
        result_data.push(*v1);
    }
    Array2::from_shape_vec((n_edges, 2), result_data).expect("Edge list shape mismatch")
}

/// Extract Elem2Edge: element to incident edge indices.
///
/// Returns Array2<i32> [n_elems, 3|4] with edge IDs for each element.
pub fn to_elem2edge(
    elem2vert: &Array2<i32>,
) -> Array2<i32> {
    let n_elems = elem2vert.shape()[0];
    let elem_cols = elem2vert.shape()[1];

    // Build canonical edge list (same order as to_edge2vert)
    let edges = to_canonical_edge_list(elem2vert);

    // Create edge_to_id map
    let mut edge_to_id: HashMap<(i32, i32), usize> = HashMap::new();
    for (i, edge) in edges.iter().enumerate() {
        edge_to_id.insert(*edge, i);
    }

    // For each element, collect its edge IDs
    let mut elem2edge: Vec<Vec<i32>> = Vec::new();
    for elem_idx in 0..n_elems {
        let elem_type = if elem_cols == 3 {
            3
        } else {
            if elem2vert[[elem_idx, 3]] == elem2vert[[elem_idx, 0]] {
                3
            } else {
                4
            }
        };

        let mut elem_edges = Vec::new();
        for j in 0..elem_type {
            let v1 = elem2vert[[elem_idx, j]];
            let v2 = elem2vert[[elem_idx, (j + 1) % elem_type]];
            let edge = (v1.min(v2), v1.max(v2));
            let edge_id = edge_to_id.get(&edge).copied().unwrap_or(usize::MAX);
            elem_edges.push(edge_id as i32);
        }
        elem2edge.push(elem_edges);
    }

    // Determine max edges per element and pad with -1
    let max_edges = elem2edge.iter().map(|e| e.len()).max().unwrap_or(0);
    let mut result = Array2::from_elem((n_elems, max_edges), -1i32);

    for (i, edges) in elem2edge.iter().enumerate() {
        for (j, &edge_id) in edges.iter().enumerate() {
            result[[i, j]] = edge_id;
        }
    }

    result
}

/// Extract Edge2Elem: edge to adjacent elements.
///
/// Returns Array2<i32> [n_edges, 2] with element indices (-1 for boundary).
pub fn to_edge2elem(elem2vert: &Array2<i32>) -> Array2<i32> {
    let n_elems = elem2vert.shape()[0];
    let elem_cols = elem2vert.shape()[1];

    // Build canonical edge list
    let edges = to_canonical_edge_list(elem2vert);

    // Build undirected-edge-to-elements map
    let mut edge_to_elems: HashMap<(i32, i32), Vec<usize>> = HashMap::new();
    for edge in edges.iter() {
        edge_to_elems.insert(*edge, Vec::new());
    }

    // Iterate elements to build mapping
    for elem_idx in 0..n_elems {
        let elem_type = if elem_cols == 3 {
            3
        } else {
            if elem2vert[[elem_idx, 3]] == elem2vert[[elem_idx, 0]] {
                3
            } else {
                4
            }
        };

        for i in 0..elem_type {
            let v1 = elem2vert[[elem_idx, i]];
            let v2 = elem2vert[[elem_idx, (i + 1) % elem_type]];
            let undirected_edge = (v1.min(v2), v1.max(v2));

            if let Some(elems) = edge_to_elems.get_mut(&undirected_edge) {
                if elems.len() < 2 {
                    elems.push(elem_idx);
                }
            }
        }
    }

    // Convert to result array with padding for boundary edges
    let mut result_data = Vec::new();
    for edge in edges.iter() {
        let elems = &edge_to_elems[edge];
        result_data.push(elems.get(0).copied().unwrap_or(usize::MAX) as i32);
        result_data.push(elems.get(1).copied().unwrap_or(usize::MAX) as i32);
    }

    let n_edges = edges.len();
    Array2::from_shape_vec((n_edges, 2), result_data).expect("Edge2Elem shape mismatch")
}

/// Extract Vert2Edge: vertex to incident edge indices.
///
/// Returns Vec<Vec<usize>> (List[List[int]] in Python for PyO3 conversion).
pub fn to_vert2edge(elem2vert: &Array2<i32>, n_verts: usize) -> Vec<Vec<usize>> {
    let edges = to_canonical_edge_list(elem2vert);

    // Create edge_to_id map
    let mut edge_to_id: HashMap<(i32, i32), usize> = HashMap::new();
    for (i, edge) in edges.iter().enumerate() {
        edge_to_id.insert(*edge, i);
    }

    // Initialize result with empty sets (for deduplication)
    let mut vert2edge_sets: Vec<std::collections::HashSet<usize>> =
        vec![std::collections::HashSet::new(); n_verts];

    // For each edge, add it to both endpoint vertices (using set for deduplication)
    for edge in edges.iter() {
        let edge_id = edge_to_id[edge];
        vert2edge_sets[edge.0 as usize].insert(edge_id);
        vert2edge_sets[edge.1 as usize].insert(edge_id);
    }

    // Convert sets to sorted vectors
    let vert2edge: Vec<Vec<usize>> = vert2edge_sets
        .iter()
        .map(|set| {
            let mut vec: Vec<usize> = set.iter().copied().collect();
            vec.sort();
            vec
        })
        .collect();

    vert2edge
}

/// Extract Vert2Elem: vertex to incident elements.
///
/// Returns Vec<Vec<usize>> (List[List[int]] in Python for PyO3 conversion).
pub fn to_vert2elem(elem2vert: &Array2<i32>, n_verts: usize) -> Vec<Vec<usize>> {
    let n_elems = elem2vert.shape()[0];
    let elem_cols = elem2vert.shape()[1];

    // Initialize result with empty sets (for deduplication)
    let mut vert2elem_sets: Vec<std::collections::HashSet<usize>> =
        vec![std::collections::HashSet::new(); n_verts];

    // For each element, add it to all incident vertices (using set for deduplication)
    for elem_idx in 0..n_elems {
        let elem_type = if elem_cols == 3 {
            3
        } else {
            if elem2vert[[elem_idx, 3]] == elem2vert[[elem_idx, 0]] {
                3
            } else {
                4
            }
        };

        for i in 0..elem_type {
            let v = elem2vert[[elem_idx, i]] as usize;
            if v < n_verts {
                vert2elem_sets[v].insert(elem_idx);
            }
        }
    }

    // Convert sets to sorted vectors
    let vert2elem: Vec<Vec<usize>> = vert2elem_sets
        .iter()
        .map(|set| {
            let mut vec: Vec<usize> = set.iter().copied().collect();
            vec.sort();
            vec
        })
        .collect();

    vert2elem
}

/// Extract Elem2Vert: return connectivity array as-is.
///
/// This is a passthrough converter (no computation needed).
pub fn to_elem2vert(elem2vert: &Array2<i32>) -> Array2<i32> {
    elem2vert.clone()
}

/// Unique undirected edges in first-encounter order.
///
/// Iterates elements in order and slots in order, taking each edge's canonical
/// (min,max) form the first time it is seen. Mirrors the Python reference
/// `_build_adjacencies` (element-major, slot-minor first-occurrence numbering);
/// the C++ backend uses the same ordering, so edge IDs — and everything keyed on
/// them (bEdgeIDs, Vert2Edge) — match across all three backends. (#163)
fn first_encounter_canonical_edges(elem2vert: &Array2<i32>) -> Vec<(i32, i32)> {
    let n_elems = elem2vert.shape()[0];
    let elem_cols = elem2vert.shape()[1];

    let mut edges: Vec<(i32, i32)> = Vec::new();
    let mut seen: std::collections::HashSet<(i32, i32)> = std::collections::HashSet::new();

    for elem_idx in 0..n_elems {
        let elem_type = if elem_cols == 3 {
            3
        } else if elem2vert[[elem_idx, 3]] == elem2vert[[elem_idx, 0]] {
            3
        } else {
            4
        };

        for i in 0..elem_type {
            let v0 = elem2vert[[elem_idx, i]];
            let v1 = elem2vert[[elem_idx, (i + 1) % elem_type]];
            if v0 < 0 || v1 < 0 || v0 == v1 {
                continue;
            }
            let edge = (v0.min(v1), v0.max(v1));
            if seen.insert(edge) {
                edges.push(edge);
            }
        }
    }
    edges
}

/// Extract canonical edge list in first-encounter order.
///
/// Returns list of unique undirected edges as (min_vert, max_vert) tuples,
/// in first-encounter element-traversal order, matching Python/C++ edge IDs.
fn to_canonical_edge_list(elem2vert: &Array2<i32>) -> Vec<(i32, i32)> {
    first_encounter_canonical_edges(elem2vert)
}
