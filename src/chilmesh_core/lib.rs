pub mod errors;
pub mod io;
pub mod adjacency;

use pyo3::prelude::*;
use ndarray::Array2;
use numpy::PyArray2;

/// RustMesh represents a 2D mesh with mixed element types (triangles and quads)
#[pyclass]
pub struct RustMesh {
    pub points: Array2<f64>,         // [n_verts, 2] — vertex coordinates
    pub connectivity: Array2<i32>,   // [n_elems, 3|4] — element vertex IDs
    pub elem_type: Vec<u32>,         // Element types (3=triangle, 4=quad)
    pub num_verts: usize,
    pub num_elems: usize,
    pub edges: Option<Array2<i32>>,  // Quad-edge topology: [n_edges, 4]
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

    #[getter]
    fn n_verts(&self) -> usize {
        self.num_verts
    }

    #[getter]
    fn n_elems(&self) -> usize {
        self.num_elems
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
}

#[pymodule]
fn chilmesh_core(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<RustMesh>()?;
    Ok(())
}
