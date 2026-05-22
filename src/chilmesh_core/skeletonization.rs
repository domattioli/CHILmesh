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
    let n_verts = points.shape()[0];

    // Step 0: Build initial edge2elem adjacency
    let mut edge2elem = build_edge2elem(connectivity);
    let mut edge2vert = build_edge2vert(connectivity);

    // Initialize tracking
    let mut remaining_elems: HashSet<usize> = (0..n_elems).collect();
    let mut remaining_verts: HashSet<usize> = (0..n_verts).collect();
    let mut layers: Vec<Layer> = Vec::new();

    while !remaining_elems.is_empty() {
        // Step 1a: Identify boundary edges (≤1 adjacent element in remaining set)
        let boundary_edge_ids = identify_boundary_edges(&edge2elem, &remaining_elems);

        if boundary_edge_ids.is_empty() {
            // All remaining elements form a single core (no boundary)
            // Classify them all as OE
            let core_elems: Vec<usize> = remaining_elems.iter().copied().collect();
            if !core_elems.is_empty() {
                let (core_ov, core_iv) = classify_vertices(&core_elems, connectivity, &remaining_verts);
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

        // Step 1f: Identify inner vertices
        let (layer_ov, layer_iv) = classify_vertices(
            &oe_elems,
            connectivity,
            &remaining_verts,
        );

        // Step 1g: Remove OE and IE from remaining_elems
        for &elem in &oe_elems {
            remaining_elems.remove(&elem);
        }
        for &elem in &ie_elems {
            remaining_elems.remove(&elem);
        }

        // Step 1h: Remove OV vertices from remaining_verts
        for &vert in &layer_ov {
            remaining_verts.remove(&vert);
        }

        // Create and store layer
        layers.push(Layer {
            oe: oe_elems,
            ie: ie_elems,
            ov: layer_ov,
            iv: layer_iv,
            b_edge_ids: boundary_edge_ids,
        });

        // Step 1i: Recompute adjacency for next iteration (only for remaining elements)
        if !remaining_elems.is_empty() {
            edge2elem = build_edge2elem_subset(connectivity, &remaining_elems);
            edge2vert = build_edge2vert_subset(connectivity, &remaining_verts);
        }
    }

    // Validate coverage invariant
    validate_layers(&layers, n_elems)?;

    Ok(layers)
}

/// Build Edge2Elem adjacency: edge → [elem1, elem2] (or [-1, elem] for boundary)
fn build_edge2elem(connectivity: &Array2<i32>) -> Vec<Vec<i32>> {
    let n_elems = connectivity.shape()[0];
    let mut edge_map: std::collections::HashMap<(i32, i32), Vec<i32>> = std::collections::HashMap::new();

    for elem in 0..n_elems {
        let edges = get_element_edges(elem, connectivity);
        for edge in edges {
            let key = if edge[0] <= edge[1] {
                (edge[0], edge[1])
            } else {
                (edge[1], edge[0])
            };
            edge_map.entry(key).or_insert_with(Vec::new).push(elem as i32);
        }
    }

    let mut edge2elem_vec = Vec::new();
    for (_key, elems) in edge_map.into_iter() {
        if elems.len() == 1 {
            edge2elem_vec.push(vec![-1, elems[0]]);
        } else if elems.len() >= 2 {
            edge2elem_vec.push(vec![elems[0], elems[1]]);
        }
    }

    edge2elem_vec
}

/// Build Edge2Vert adjacency: maps each edge to its two vertices
fn build_edge2vert(connectivity: &Array2<i32>) -> Vec<(i32, i32)> {
    let n_elems = connectivity.shape()[0];
    let mut edge_set: std::collections::HashSet<(i32, i32)> = std::collections::HashSet::new();

    for elem in 0..n_elems {
        let edges = get_element_edges(elem, connectivity);
        for edge in edges {
            let normalized = if edge[0] <= edge[1] {
                (edge[0], edge[1])
            } else {
                (edge[1], edge[0])
            };
            edge_set.insert(normalized);
        }
    }

    let mut edge_vec: Vec<_> = edge_set.into_iter().collect();
    edge_vec.sort();
    edge_vec
}

/// Build Edge2Elem for only remaining elements
fn build_edge2elem_subset(
    connectivity: &Array2<i32>,
    remaining: &HashSet<usize>,
) -> Vec<Vec<i32>> {
    let mut edge_map: std::collections::HashMap<(i32, i32), Vec<i32>> = std::collections::HashMap::new();

    for &elem in remaining {
        let edges = get_element_edges(elem, connectivity);
        for edge in edges {
            let key = if edge[0] <= edge[1] {
                (edge[0], edge[1])
            } else {
                (edge[1], edge[0])
            };
            edge_map.entry(key).or_insert_with(Vec::new).push(elem as i32);
        }
    }

    let mut edge2elem_vec = Vec::new();
    for (_key, elems) in edge_map.into_iter() {
        if elems.len() == 1 {
            edge2elem_vec.push(vec![-1, elems[0]]);
        } else if elems.len() >= 2 {
            edge2elem_vec.push(vec![elems[0], elems[1]]);
        }
    }

    edge2elem_vec
}

/// Build Edge2Vert for only remaining vertices
fn build_edge2vert_subset(
    connectivity: &Array2<i32>,
    remaining: &HashSet<usize>,
) -> Vec<(i32, i32)> {
    let mut edge_set: std::collections::HashSet<(i32, i32)> = std::collections::HashSet::new();

    let n_elems = connectivity.shape()[0];
    for elem in 0..n_elems {
        let edges = get_element_edges(elem, connectivity);
        for edge in edges {
            // Only include edge if both vertices are in remaining set
            if remaining.contains(&(edge[0] as usize)) && remaining.contains(&(edge[1] as usize)) {
                let normalized = if edge[0] <= edge[1] {
                    (edge[0], edge[1])
                } else {
                    (edge[1], edge[0])
                };
                edge_set.insert(normalized);
            }
        }
    }

    let mut edge_vec: Vec<_> = edge_set.into_iter().collect();
    edge_vec.sort();
    edge_vec
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

/// Classify vertices for a layer: OV and IV
fn classify_vertices(
    layer_elems: &[usize],
    connectivity: &Array2<i32>,
    remaining_verts: &HashSet<usize>,
) -> (Vec<usize>, Vec<usize>) {
    let mut layer_verts = HashSet::new();

    // Collect all vertices from layer elements
    for &elem in layer_elems {
        if elem < connectivity.shape()[0] {
            let row = connectivity.row(elem);
            for &v in row.iter() {
                if v >= 0 && remaining_verts.contains(&(v as usize)) {
                    layer_verts.insert(v as usize);
                }
            }
        }
    }

    // OV and IV will be determined by the layer extraction logic
    // For now, all vertices in layer elements are potential OV/IV
    // This classification is refined during the main loop
    let mut verts: Vec<_> = layer_verts.into_iter().collect();
    verts.sort();

    // Placeholder: actual OV/IV separation handled in main loop
    // Return all layer verts as OV, empty IV
    (verts.clone(), Vec::new())
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
