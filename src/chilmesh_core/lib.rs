pub mod errors;
pub mod io;

use pyo3::prelude::*;
use ndarray::Array2;

/// RustMesh represents a 2D mesh with mixed element types (triangles and quads)
#[pyclass]
pub struct RustMesh {
    pub points: Array2<f64>,         // [n_verts, 2] — vertex coordinates
    pub connectivity: Array2<i32>,   // [n_elems, 3|4] — element vertex IDs
    pub elem_type: Vec<u32>,         // Element types (3=triangle, 4=quad)
    pub num_verts: usize,
    pub num_elems: usize,
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
}

#[pymodule]
fn chilmesh_core(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<RustMesh>()?;
    Ok(())
}
