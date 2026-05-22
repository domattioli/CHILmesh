use thiserror::Error;
use pyo3::prelude::*;
use pyo3::exceptions;

/// Error types for RustMesh operations
#[derive(Error, Debug)]
pub enum RustMeshError {
    #[error("Fort.14 parse error at line {line}: {reason}")]
    ParseError { line: usize, reason: String },

    #[error("Invalid geometry: {0}")]
    InvalidGeometry(String),

    #[error("Vertex {vertex_id} out of bounds [0, {max_vertices})")]
    VertexOutOfBounds { vertex_id: i32, max_vertices: usize },

    #[error("Invalid element {elem_id}: {reason}")]
    InvalidElement { elem_id: usize, reason: String },

    #[error("I/O error: {0}")]
    IOError(String),
}

impl From<RustMeshError> for PyErr {
    fn from(err: RustMeshError) -> PyErr {
        match err {
            RustMeshError::ParseError { line, reason } => {
                let msg = format!("Fort.14 parse error at line {}: {}", line, reason);
                PyErr::new::<exceptions::PyOSError, _>(msg)
            }
            RustMeshError::IOError(msg) => PyErr::new::<exceptions::PyOSError, _>(msg),
            RustMeshError::InvalidGeometry(msg) => {
                PyErr::new::<exceptions::PyValueError, _>(msg)
            }
            RustMeshError::VertexOutOfBounds {
                vertex_id,
                max_vertices,
            } => {
                let msg = format!(
                    "Vertex {} out of bounds [0, {})",
                    vertex_id, max_vertices
                );
                PyErr::new::<exceptions::PyValueError, _>(msg)
            }
            RustMeshError::InvalidElement { elem_id, reason } => {
                let msg = format!("Invalid element {}: {}", elem_id, reason);
                PyErr::new::<exceptions::PyValueError, _>(msg)
            }
        }
    }
}

impl From<std::io::Error> for RustMeshError {
    fn from(err: std::io::Error) -> Self {
        RustMeshError::IOError(err.to_string())
    }
}
