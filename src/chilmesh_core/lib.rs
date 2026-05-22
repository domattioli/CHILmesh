pub mod errors;
pub mod io;
pub mod adjacency;
pub mod queries;
pub mod skeletonization;
pub mod mutation;

use pyo3::prelude::*;
use ndarray::{Array1, Array2};
use numpy::{PyArray1, PyArray2};

/// RustMesh represents a 2D mesh with mixed element types (triangles and quads)
#[pyclass]
pub struct RustMesh {
    pub points: Array2<f64>,         // [n_verts, 2] — vertex coordinates
    pub connectivity: Array2<i32>,   // [n_elems, 3|4] — element vertex IDs
    pub elem_type: Vec<u32>,         // Element types (3=triangle, 4=quad)
    pub num_verts: usize,
    pub num_elems: usize,
    pub edges: Option<Array2<i32>>,  // Quad-edge topology: [n_edges, 4]
    pub areas: Option<Array1<f64>>,  // Signed areas: [n_elems] (computed on demand)
    pub layers: Option<Vec<skeletonization::Layer>>, // Skeletonization layers (computed on demand)
}

#[pymethods]
impl RustMesh {
    #[new]
    fn new() -> Self {
        RustMesh {
            points: Array2::zeros((0, 2)),
            connectivity: Array2::zeros((0, 4)),
            elem_type: Vec::new(),
            num_verts: 0,
            num_elems: 0,
            edges: None,
            areas: None,
            layers: None,
        }
    }

    fn read_from_fort14(&mut self, path: &str) -> PyResult<()> {
        let mesh = io::parse_fort14(path)?;
        self.points = mesh.points;
        self.connectivity = mesh.connectivity;
        self.elem_type = mesh.elem_type;
        self.num_verts = mesh.num_verts;
        self.num_elems = mesh.num_elems;
        Ok(())
    }

    fn write_fort14(&self, path: &str) -> PyResult<()> {
        io::write_fort14(self, path)?;
        Ok(())
    }

    fn read_from_2dm(&mut self, path: &str) -> PyResult<()> {
        let mesh = io::parse_2dm(path)?;
        self.points = mesh.points;
        self.connectivity = mesh.connectivity;
        self.elem_type = mesh.elem_type;
        self.num_verts = mesh.num_verts;
        self.num_elems = mesh.num_elems;
        Ok(())
    }

    #[getter]
    fn n_verts(&self) -> usize {
        self.num_verts
    }

    #[getter]
    fn n_elems(&self) -> usize {
        self.num_elems
    }

    /// Set points from a numpy array (for testing)
    fn set_points(&mut self, py: Python, points: Py<PyArray2<f64>>) -> PyResult<()> {
        let array = points.as_ref(py).to_owned_array();
        self.points = array;
        self.num_verts = self.points.shape()[0];
        Ok(())
    }

    /// Set connectivity from a numpy array (for testing)
    fn set_connectivity(&mut self, py: Python, connectivity: Py<PyArray2<i32>>) -> PyResult<()> {
        let array = connectivity.as_ref(py).to_owned_array();
        self.connectivity = array;
        self.num_elems = self.connectivity.shape()[0];
        Ok(())
    }

    /// Get points as numpy array
    fn get_points<'py>(&self, py: Python<'py>) -> PyResult<Py<PyArray2<f64>>> {
        let result = self.points.clone();
        Ok(PyArray2::from_owned_array(py, result).to_owned())
    }

    /// Get connectivity as numpy array
    fn get_connectivity<'py>(&self, py: Python<'py>) -> PyResult<Py<PyArray2<i32>>> {
        let result = self.connectivity.clone();
        Ok(PyArray2::from_owned_array(py, result).to_owned())
    }

    /// Build adjacencies (quad-edge construction + converters)
    fn build_adjacencies(&mut self) -> PyResult<()> {
        let edges = adjacency::build_quadegg_from_connectivity(&self.connectivity, self.num_verts);
        self.edges = Some(edges);
        Ok(())
    }

    /// Get Edge2Vert: returns ndarray [n_edges, 2]
    fn get_edge2vert<'py>(&self, py: Python<'py>) -> PyResult<Py<PyArray2<i32>>> {
        let result = adjacency::to_edge2vert(&self.connectivity);
        Ok(PyArray2::from_owned_array(py, result).to_owned())
    }

    /// Get Elem2Edge: returns ndarray [n_elems, 3|4]
    fn get_elem2edge<'py>(&self, py: Python<'py>) -> PyResult<Py<PyArray2<i32>>> {
        let result = adjacency::to_elem2edge(&self.connectivity);
        Ok(PyArray2::from_owned_array(py, result).to_owned())
    }

    /// Get Vert2Edge: returns PyList of lists
    fn get_vert2edge(&self) -> PyResult<Vec<Vec<usize>>> {
        let result = adjacency::to_vert2edge(&self.connectivity, self.num_verts);
        Ok(result)
    }

    /// Get Vert2Elem: returns PyList of lists
    fn get_vert2elem(&self) -> PyResult<Vec<Vec<usize>>> {
        let result = adjacency::to_vert2elem(&self.connectivity, self.num_verts);
        Ok(result)
    }

    /// Get Edge2Elem: returns ndarray [n_edges, 2]
    fn get_edge2elem<'py>(&self, py: Python<'py>) -> PyResult<Py<PyArray2<i32>>> {
        let result = adjacency::to_edge2elem(&self.connectivity);
        Ok(PyArray2::from_owned_array(py, result).to_owned())
    }

    /// Get Elem2Vert: returns connectivity ndarray [n_elems, 3|4]
    fn get_elem2vert<'py>(&self, py: Python<'py>) -> PyResult<Py<PyArray2<i32>>> {
        let result = adjacency::to_elem2vert(&self.connectivity);
        Ok(PyArray2::from_owned_array(py, result).to_owned())
    }

    /// Compute signed areas for all elements using shoelace formula
    fn compute_quality(&mut self) -> PyResult<()> {
        let areas = queries::compute_signed_areas(&self.points, &self.connectivity);
        self.areas = Some(areas);
        Ok(())
    }

    /// Ensure counter-clockwise orientation for all elements
    fn ensure_ccw(&mut self) -> PyResult<()> {
        queries::ensure_ccw_orientation(&self.points, &mut self.connectivity);
        // Recompute areas after reorientation
        let areas = queries::compute_signed_areas(&self.points, &self.connectivity);
        self.areas = Some(areas);
        Ok(())
    }

    /// Get signed areas array (must call compute_quality first)
    fn get_signed_areas<'py>(&self, py: Python<'py>) -> PyResult<Py<PyArray1<f64>>> {
        match &self.areas {
            Some(areas) => {
                let result = areas.clone();
                Ok(PyArray1::from_owned_array(py, result).to_owned())
            }
            None => Err(pyo3::exceptions::PyRuntimeError::new_err(
                "Signed areas not computed. Call compute_quality() first.",
            )),
        }
    }

    /// Get edges incident to a vertex (requires edges to be computed)
    fn get_vertex_edges(&self, v: usize) -> PyResult<Vec<usize>> {
        match &self.edges {
            Some(_edges_array) => {
                // Convert quad-edge format [n_edges, 4] to edge list [n_edges, 2]
                let edge2vert = adjacency::to_edge2vert(&self.connectivity);
                Ok(queries::get_vertex_edges(v, &edge2vert))
            }
            None => Err(pyo3::exceptions::PyRuntimeError::new_err(
                "Edges not computed. Call build_adjacencies() first.",
            )),
        }
    }

    /// Get elements incident to a vertex
    fn get_vertex_elements(&self, v: usize) -> PyResult<Vec<usize>> {
        Ok(queries::get_vertex_elements(v, &self.connectivity))
    }

    /// Get vertices of an element as 4-element padded array
    fn get_element_vertices<'py>(
        &self,
        py: Python<'py>,
        e: usize,
    ) -> PyResult<Py<PyArray1<i32>>> {
        let result = queries::get_element_vertices(e, &self.connectivity);
        Ok(PyArray1::from_owned_array(py, result).to_owned())
    }

    /// Skeletonize the mesh by layer-by-layer boundary removal
    fn skeletonize(&mut self, quality_threshold: Option<f64>) -> PyResult<()> {
        let layers = skeletonization::skeletonize_medial_axis(
            &self.connectivity,
            &self.points,
            quality_threshold,
        ).map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;

        validate_layers(&layers, self.num_elems)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;

        self.layers = Some(layers);
        Ok(())
    }

    /// Get number of layers (must call skeletonize first)
    fn get_num_layers(&self) -> PyResult<usize> {
        match &self.layers {
            Some(layers) => Ok(layers.len()),
            None => Err(pyo3::exceptions::PyRuntimeError::new_err(
                "Layers not computed. Call skeletonize() first.",
            )),
        }
    }

    /// Get a specific layer as a Python dict
    fn get_layer<'py>(&self, py: Python<'py>, layer_idx: usize) -> PyResult<Py<pyo3::types::PyDict>> {
        match &self.layers {
            Some(layers) => {
                if layer_idx >= layers.len() {
                    return Err(pyo3::exceptions::PyValueError::new_err(
                        format!("Layer index {} out of range [0, {}]", layer_idx, layers.len() - 1)
                    ));
                }

                let layer = &layers[layer_idx];

                // Create a Python dict with the layer data
                let dict = pyo3::types::PyDict::new_bound(py);

                // Convert Vec<usize> to PyList
                let oe: Vec<usize> = layer.oe.clone();
                let ie: Vec<usize> = layer.ie.clone();
                let ov: Vec<usize> = layer.ov.clone();
                let iv: Vec<usize> = layer.iv.clone();
                let b_edge_ids: Vec<usize> = layer.b_edge_ids.clone();

                dict.set_item("OE", oe)?;
                dict.set_item("IE", ie)?;
                dict.set_item("OV", ov)?;
                dict.set_item("IV", iv)?;
                dict.set_item("bEdgeIDs", b_edge_ids)?;

                Ok(dict.unbind())
            }
            None => Err(pyo3::exceptions::PyRuntimeError::new_err(
                "Layers not computed. Call skeletonize() first.",
            )),
        }
    }

    /// Get all layers as a list of dicts
    fn get_all_layers<'py>(&self, py: Python<'py>) -> PyResult<Vec<Py<pyo3::types::PyDict>>> {
        match &self.layers {
            Some(layers) => {
                let mut result = Vec::new();
                for (layer_idx, _layer) in layers.iter().enumerate() {
                    // Use get_layer to construct each layer dict
                    result.push(self.get_layer(py, layer_idx)?);
                }
                Ok(result)
            }
            None => Err(pyo3::exceptions::PyRuntimeError::new_err(
                "Layers not computed. Call skeletonize() first.",
            )),
        }
    }

    /// Add a new element to the mesh
    ///
    /// Validates vertex bounds, element type, and geometric constraints.
    /// Returns the new element ID.
    ///
    /// # Arguments
    /// * `verts` - Vertex indices [v0, v1, v2] or [v0, v1, v2, v3]
    /// * `elem_type_new` - Element type (3 for triangle, 4 for quad)
    fn add_element(&self, verts: Vec<i32>, elem_type_new: u32) -> PyResult<usize> {
        let result = mutation::add_element(
            &self.connectivity,
            &self.points,
            self.num_verts,
            self.num_elems,
            &self.elem_type,
            &verts,
            elem_type_new,
        )?;
        Ok(result)
    }

    /// Remove an element from the mesh
    ///
    /// Validates element bounds and checks for orphaned vertices.
    ///
    /// # Arguments
    /// * `elem_id` - Element ID to remove (0-indexed)
    fn remove_element(&self, elem_id: usize) -> PyResult<()> {
        mutation::remove_element(&self.connectivity, self.num_elems, elem_id)?;
        Ok(())
    }
}

/// Validate layer invariants
fn validate_layers(layers: &[skeletonization::Layer], total_elems: usize) -> Result<(), String> {
    let mut elem_count = 0;
    let mut all_elems = std::collections::HashSet::new();

    for (layer_idx, layer) in layers.iter().enumerate() {
        // Check disjoint OE/IE
        let oe_set: std::collections::HashSet<_> = layer.oe.iter().copied().collect();
        let ie_set: std::collections::HashSet<_> = layer.ie.iter().copied().collect();

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

#[pymodule]
fn chilmesh_core(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<RustMesh>()?;
    Ok(())
}
