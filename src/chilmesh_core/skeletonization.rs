//! Mesh skeletonization — layer-by-layer medial axis extraction.
//!
//! This module implements layer-by-layer boundary removal, extracting concentric
//! layers from the inside-out. Each layer is classified into:
//! - OE (Oriented Elements): boundary elements
//! - IE (Inner Elements): interior elements
//! - OV (Outer Vertices): boundary vertices
//! - IV (Inner Vertices): interior vertices
//!
//! The algorithm maintains invariants:
//! - sum(OE + IE) across all layers == total elements
//! - OE and IE are disjoint
//! - OV and IV are disjoint
//! - All elements classified exactly once

use ndarray::Array2;
use std::collections::HashSet;

/// Represents a single layer in mesh skeletonization
#[derive(Clone, Debug)]
pub struct Layer {
    pub oe: Vec<usize>,        // Oriented elements (boundary)
    pub ie: Vec<usize>,        // Inner elements
    pub ov: Vec<usize>,        // Outer vertices
    pub iv: Vec<usize>,        // Inner vertices
    pub b_edge_ids: Vec<usize>, // Boundary edges defining this layer's frontier
}

/// Translate per-iteration subset boundary-edge indices into full-mesh edge IDs,
/// sorted ascending to match Python's `np.where`-ordered `bEdgeIDs` (#163).
fn to_full_mesh_edge_ids(
    subset_boundary_ids: &[usize],
    subset_edge2vert: &[(i32, i32)],
    full_edge_id: &std::collections::HashMap<(i32, i32), usize>,
) -> Vec<usize> {
    let mut out: Vec<usize> = subset_boundary_ids
        .iter()
        .filter_map(|&bid| subset_edge2vert.get(bid).and_then(|key| full_edge_id.get(key).copied()))
        .collect();
    out.sort_unstable();
    out
}

/// Skeletonize a mesh by layer-by-layer boundary removal.
///
/// Algorithm:
/// 1. Initialize tracking: remaining_elems, remaining_verts, layers
/// 2. While remaining_elems not empty:
///    a. Identify boundary edges (edges with ≤1 adjacent element)
///    b. Mark boundary elements as OE
///    c. Identify outer vertices (incident to boundary edges) as OV
///    d. Identify inner elements (no boundary edges) as IE
///    e. Identify inner vertices as IV = remaining_verts - OV
///    f. Create Layer { oe, ie, ov, iv }
///    g. Remove boundary elements from remaining_elems
///    h. Remove orphaned vertices
///    i. Recompute adjacency for next iteration
///
/// Returns: Vec<Layer> with proper invariants maintained
pub fn skeletonize_medial_axis(
    connectivity: &Array2<i32>,
    points: &Array2<f64>,
    _quality_threshold: Option<f64>,
) -> Result<Vec<Layer>, String> {
    let n_elems = connectivity.shape()[0];
    let _n_verts = points.shape()[0];

    // Full-mesh canonical edge IDs (first-encounter order; matches Python Edge2Vert
    // and thus Python Layers["bEdgeIDs"]. Used to translate the per-iteration
    // subset boundary-edge indices into stable full-mesh edge IDs (#163).
    let full_edges = crate::adjacency::first_encounter_canonical_edges(connectivity);
    let mut full_edge_id: std::collections::HashMap<(i32, i32), usize> =
        std::collections::HashMap::with_capacity(full_edges.len());
    for (i, &e) in full_edges.iter().enumerate() {
        full_edge_id.insert(e, i);
    }

    // Step 0: Build initial aligned (edge2vert, edge2elem) — shared edge-id order
    let mut remaining_elems: HashSet<usize> = (0..n_elems).collect();
    let (mut edge2vert, mut edge2elem) = build_edge_adjacency_subset(connectivity, &remaining_elems);
    let mut layers: Vec<Layer> = Vec::new();

    while !remaining_elems.is_empty() {
        // Step 1a: Identify boundary edges (≤1 adjacent element in remaining set)
        let boundary_edge_ids = identify_boundary_edges(&edge2elem, &remaining_elems);

        if boundary_edge_ids.is_empty() {
            // All remaining elements form a single core (no boundary)
            // Classify them all as OE
            let mut core_elems: Vec<usize> = remaining_elems.iter().copied().collect();
            core_elems.sort();
            if !core_elems.is_empty() {
                let (core_ov, core_iv) = classify_vertices_for_layer(
                    &[],  // No boundary edges for the core
                    &edge2vert,
                    &core_elems,
                    connectivity,
                );
                layers.push(Layer {
                    oe: core_elems,
                    ie: Vec::new(),
                    ov: core_ov,
                    iv: core_iv,
                    b_edge_ids: Vec::new(),
                });
                remaining_elems.clear();
            }
            break;
        }

        // Step 1b: Mark boundary elements as OE
        let oe_elems = classify_boundary_elements(&boundary_edge_ids, &edge2elem, &remaining_elems);

        // Step 1c: Identify outer vertices (incident to boundary edges)
        let ov_verts = classify_outer_vertices(&boundary_edge_ids, &edge2vert);

        // Step 1d: Find all edges touching any OV vertex
        let ov_edge_indices = find_edges_with_vertices(&edge2vert, &ov_verts);

        // Step 1e: Identify inner elements (adjacent to OV edges but not OE)
        let ie_elems = classify_inner_elements(
            &ov_edge_indices,
            &edge2elem,
            &oe_elems,
            &remaining_elems,
        );

        // Step 1f: Identify outer and inner vertices for this layer
        // Outer vertices = vertices on boundary edges (from Step 1c)
        // Inner vertices = vertices of (OE ∪ IE) that are not OV
        let all_layer_elems: Vec<usize> = oe_elems
            .iter()
            .chain(ie_elems.iter())
            .copied()
            .collect();
        let (layer_ov, layer_iv) = classify_vertices_for_layer(
            &boundary_edge_ids,
            &edge2vert,
            &all_layer_elems,
            connectivity,
        );

        // Step 1g: Remove OE and IE from remaining_elems
        for &elem in &oe_elems {
            remaining_elems.remove(&elem);
        }
        for &elem in &ie_elems {
            remaining_elems.remove(&elem);
        }

        // Create and store layer
        let layer_b_edge_ids = to_full_mesh_edge_ids(&boundary_edge_ids, &edge2vert, &full_edge_id);
        layers.push(Layer {
            oe: oe_elems,
            ie: ie_elems,
            ov: layer_ov,
            iv: layer_iv,
            b_edge_ids: layer_b_edge_ids,
        });

        // Step 1i: Recompute aligned adjacency for the remaining elements
        if !remaining_elems.is_empty() {
            let (ev, ee) = build_edge_adjacency_subset(connectivity, &remaining_elems);
            edge2vert = ev;
            edge2elem = ee;
        }
    }

    // Validate coverage invariant
    validate_layers(&layers, n_elems)?;

    Ok(layers)
}

/// Build aligned (edge2vert, edge2elem) restricted to `remaining` elements.
///
/// Both returned vecs share ONE edge ordering: index `i` denotes the same
/// physical edge in `edge2vert[i]` and `edge2elem[i]`. Order is the sorted
/// canonical (min,max) vertex key, so it is deterministic and HashMap-order
/// independent. `edge2elem[i]` is `[elemA, elemB]` for an interior edge or
/// `[-1, elem]` for a boundary edge (≤1 incident remaining element).
fn build_edge_adjacency_subset(
    connectivity: &Array2<i32>,
    remaining: &HashSet<usize>,
) -> (Vec<(i32, i32)>, Vec<Vec<i32>>) {
    let mut edge_map: std::collections::HashMap<(i32, i32), Vec<i32>> =
        std::collections::HashMap::new();
    for &elem in remaining {
        for edge in get_element_edges(elem, connectivity) {
            let key = if edge[0] <= edge[1] {
                (edge[0], edge[1])
            } else {
                (edge[1], edge[0])
            };
            edge_map.entry(key).or_insert_with(Vec::new).push(elem as i32);
        }
    }
    let mut keys: Vec<(i32, i32)> = edge_map.keys().copied().collect();
    keys.sort();
    let mut edge2vert = Vec::with_capacity(keys.len());
    let mut edge2elem = Vec::with_capacity(keys.len());
    for key in keys {
        let elems = &edge_map[&key];
        edge2vert.push(key);
        if elems.len() == 1 {
            edge2elem.push(vec![-1, elems[0]]);
        } else {
            edge2elem.push(vec![elems[0], elems[1]]);
        }
    }
    (edge2vert, edge2elem)
}

/// Get edges of an element (handles triangles and quads)
fn get_element_edges(elem_idx: usize, connectivity: &Array2<i32>) -> Vec<[i32; 2]> {
    let row = connectivity.row(elem_idx);
    let v0 = row[0];
    let v1 = row[1];
    let v2 = row[2];
    let v3 = if row.len() > 3 { row[3] } else { v2 };

    let mut edges = vec![[v0, v1], [v1, v2]];
    if v3 != v2 {
        // Quad
        edges.push([v2, v3]);
        edges.push([v3, v0]);
    } else {
        // Triangle (padded)
        edges.push([v2, v0]);
    }
    edges
}

/// Identify boundary edges: edges with ≤1 adjacent element in remaining set
fn identify_boundary_edges(
    edge2elem: &[Vec<i32>],
    remaining: &HashSet<usize>,
) -> Vec<usize> {
    let mut boundary = Vec::new();
    for (edge_id, elems) in edge2elem.iter().enumerate() {
        let active_count = elems
            .iter()
            .filter(|&&e| e >= 0 && remaining.contains(&(e as usize)))
            .count();
        if active_count <= 1 {
            boundary.push(edge_id);
        }
    }
    boundary
}

/// Classify vertices: OV from boundary edges, IV from all layer elements minus OV
fn classify_vertices_for_layer(
    boundary_edge_ids: &[usize],
    edge2vert: &[(i32, i32)],
    layer_elems: &[usize],
    connectivity: &Array2<i32>,
) -> (Vec<usize>, Vec<usize>) {
    // Outer vertices: vertices on boundary edges
    let mut ov_set = HashSet::new();
    for &edge_id in boundary_edge_ids {
        if edge_id < edge2vert.len() {
            let (v1, v2) = edge2vert[edge_id];
            ov_set.insert(v1 as usize);
            ov_set.insert(v2 as usize);
        }
    }

    // Collect all vertices from layer elements
    let mut layer_verts = HashSet::new();
    for &elem in layer_elems {
        if elem < connectivity.shape()[0] {
            let row = connectivity.row(elem);
            for &v in row.iter() {
                if v >= 0 {
                    layer_verts.insert(v as usize);
                }
            }
        }
    }

    // IV = all layer vertices - OV
    let mut iv_set: HashSet<_> = layer_verts.difference(&ov_set).copied().collect();

    let mut ov: Vec<_> = ov_set.into_iter().collect();
    let mut iv: Vec<_> = iv_set.into_iter().collect();
    ov.sort();
    iv.sort();

    (ov, iv)
}

/// Classify boundary elements: elements adjacent to boundary edges
fn classify_boundary_elements(
    boundary_edge_ids: &[usize],
    edge2elem: &[Vec<i32>],
    remaining: &HashSet<usize>,
) -> Vec<usize> {
    let mut oe_set = HashSet::new();
    for &edge_id in boundary_edge_ids {
        if edge_id < edge2elem.len() {
            for &elem_id in &edge2elem[edge_id] {
                if elem_id >= 0 && remaining.contains(&(elem_id as usize)) {
                    oe_set.insert(elem_id as usize);
                }
            }
        }
    }
    let mut oe: Vec<_> = oe_set.into_iter().collect();
    oe.sort();
    oe
}

/// Classify outer vertices: vertices incident to boundary edges
fn classify_outer_vertices(
    boundary_edge_ids: &[usize],
    edge2vert: &[(i32, i32)],
) -> HashSet<usize> {
    let mut ov = HashSet::new();
    for &edge_id in boundary_edge_ids {
        if edge_id < edge2vert.len() {
            let (v1, v2) = edge2vert[edge_id];
            ov.insert(v1 as usize);
            ov.insert(v2 as usize);
        }
    }
    ov
}

/// Find all edges that contain any of the given vertices
fn find_edges_with_vertices(
    edge2vert: &[(i32, i32)],
    vertices: &HashSet<usize>,
) -> Vec<usize> {
    let mut edge_ids = Vec::new();
    for (edge_id, &(v1, v2)) in edge2vert.iter().enumerate() {
        if vertices.contains(&(v1 as usize)) || vertices.contains(&(v2 as usize)) {
            edge_ids.push(edge_id);
        }
    }
    edge_ids
}

/// Classify inner elements: elements adjacent to OV-incident edges but not in OE
fn classify_inner_elements(
    ov_edge_indices: &[usize],
    edge2elem: &[Vec<i32>],
    oe_elems: &[usize],
    remaining: &HashSet<usize>,
) -> Vec<usize> {
    let mut ie_set = HashSet::new();
    let oe_set: HashSet<_> = oe_elems.iter().copied().collect();

    for &edge_id in ov_edge_indices {
        if edge_id < edge2elem.len() {
            for &elem_id in &edge2elem[edge_id] {
                if elem_id >= 0 {
                    let elem_usize = elem_id as usize;
                    if remaining.contains(&elem_usize) && !oe_set.contains(&elem_usize) {
                        ie_set.insert(elem_usize);
                    }
                }
            }
        }
    }

    let mut ie: Vec<_> = ie_set.into_iter().collect();
    ie.sort();
    ie
}


/// Validate layer invariants: coverage, disjoint sets
fn validate_layers(layers: &[Layer], total_elems: usize) -> Result<(), String> {
    let mut elem_count = 0;
    let mut all_elems = HashSet::new();

    for (layer_idx, layer) in layers.iter().enumerate() {
        // Check disjoint OE/IE
        let oe_set: HashSet<_> = layer.oe.iter().copied().collect();
        let ie_set: HashSet<_> = layer.ie.iter().copied().collect();

        let overlap: Vec<_> = oe_set.intersection(&ie_set).copied().collect();
        if !overlap.is_empty() {
            return Err(format!(
                "Layer {}: OE/IE overlap detected: {:?}",
                layer_idx, overlap
            ));
        }

        // Check coverage
        elem_count += layer.oe.len() + layer.ie.len();

        // Check no duplicate elements across layers
        for &elem in layer.oe.iter().chain(layer.ie.iter()) {
            if !all_elems.insert(elem) {
                return Err(format!(
                    "Element {} appears in multiple layers",
                    elem
                ));
            }
        }
    }

    // Verify total coverage
    if elem_count != total_elems {
        return Err(format!(
            "Coverage invariant violated: {} elements classified, {} total",
            elem_count, total_elems
        ));
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    /// Structured NxN quad grid: nodes (N+1)x(N+1), CCW quads.
    fn grid(n: usize) -> (Array2<f64>, Array2<i32>) {
        let np = n + 1;
        let mut pts = Vec::with_capacity(np * np * 2);
        for i in 0..np {
            for j in 0..np {
                pts.push(i as f64);
                pts.push(j as f64);
            }
        }
        let points = Array2::from_shape_vec((np * np, 2), pts).unwrap();
        let mut conn = Vec::with_capacity(n * n * 4);
        let id = |i: usize, j: usize| (i * np + j) as i32;
        for i in 0..n {
            for j in 0..n {
                conn.push(id(i, j));
                conn.push(id(i + 1, j));
                conn.push(id(i + 1, j + 1));
                conn.push(id(i, j + 1));
            }
        }
        let connectivity = Array2::from_shape_vec((n * n, 4), conn).unwrap();
        (points, connectivity)
    }

    #[test]
    fn skeletonize_peels_multiple_layers() {
        // Regression for #163: array-fed skeletonize previously returned
        // n_layers==2 for every mesh due to misaligned edge-id spaces.
        let (points, connectivity) = grid(8); // 64 quads
        let layers = skeletonize_medial_axis(&connectivity, &points, None).unwrap();
        assert!(
            layers.len() > 2,
            "expected concentric peel to yield >2 layers, got {}",
            layers.len()
        );
        // Coverage: every element classified exactly once.
        let total: usize = layers.iter().map(|l| l.oe.len() + l.ie.len()).sum();
        assert_eq!(total, connectivity.shape()[0]);
    }

    #[test]
    fn skeletonize_layer_count_grows_with_size() {
        let small = skeletonize_medial_axis(&grid(4).1, &grid(4).0, None).unwrap().len();
        let large = skeletonize_medial_axis(&grid(12).1, &grid(12).0, None).unwrap().len();
        assert!(
            large > small,
            "bigger grid must peel more layers: grid(12)={} !> grid(4)={}",
            large, small
        );
    }
}
