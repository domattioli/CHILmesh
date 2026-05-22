use crate::errors::RustMeshError;
use crate::RustMesh;
use ndarray::Array2;
use std::fs;

/// Parse Fort.14 (ADCIRC ASCII format) mesh file
pub fn parse_fort14(path: &str) -> Result<RustMesh, RustMeshError> {
    let content = fs::read_to_string(path).map_err(|e| RustMeshError::IOError(e.to_string()))?;
    let lines: Vec<&str> = content.lines().map(|l| l.trim()).collect();

    if lines.len() < 3 {
        return Err(RustMeshError::ParseError {
            line: 0,
            reason: "File too short (need header, counts, and data)".to_string(),
        });
    }

    // Line 0: Title (skip)
    // Line 1: Counts (NE NP)
    let counts_parts: Vec<&str> = lines[1].split_whitespace().collect();
    if counts_parts.len() < 2 {
        return Err(RustMeshError::ParseError {
            line: 1,
            reason: "Invalid counts format (expected NE NP)".to_string(),
        });
    }

    let n_elems: u32 = counts_parts[0].parse().map_err(|_| RustMeshError::ParseError {
        line: 1,
        reason: format!("Invalid NE (expected u32): {}", counts_parts[0]),
    })?;

    let n_verts: u32 = counts_parts[1].parse().map_err(|_| RustMeshError::ParseError {
        line: 1,
        reason: format!("Invalid NP (expected u32): {}", counts_parts[1]),
    })?;

    if n_elems == 0 || n_verts == 0 {
        return Err(RustMeshError::ParseError {
            line: 1,
            reason: format!("Invalid header: NE={}, NP={} (both must be > 0)", n_elems, n_verts),
        });
    }

    let n_elems = n_elems as usize;
    let n_verts = n_verts as usize;

    // Initialize arrays
    let mut connectivity = Array2::zeros((n_elems, 4));
    let mut elem_type = Vec::with_capacity(n_elems);
    let mut points = Array2::zeros((n_verts, 2));

    // Auto-detect format: check line 3 (index 2) to see if it's a vertex or element
    // Vertex line: ID X Y [Z]
    // Element line: ID ITYPE V1 V2 [V3 [V4]]
    // The key difference: element line has ITYPE which should be 3 or 4
    let test_parts: Vec<&str> = lines[2].split_whitespace().collect();
    let vertices_first = if test_parts.len() >= 3 {
        // Try to parse the second field as a number to determine type
        if let Ok(second_val) = test_parts[1].parse::<u32>() {
            // If it's 3 or 4, it's likely ITYPE (element), so elements first
            // Otherwise, it's likely X coordinate (vertex), so vertices first
            second_val != 3 && second_val != 4
        } else {
            // Can't parse as integer, assume it's X coordinate (vertex)
            true
        }
    } else {
        // Default to vertices first if we can't determine
        true
    };

    if vertices_first {
        // Parse nodes/coordinates first (starting at line 2, 0-indexed)
        for vert_idx in 0..n_verts {
            let line_num = 2 + vert_idx;
            if line_num >= lines.len() {
                return Err(RustMeshError::ParseError {
                    line: line_num,
                    reason: "Unexpected end of file while parsing coordinates".to_string(),
                });
            }

            let parts: Vec<&str> = lines[line_num].split_whitespace().collect();
            if parts.len() < 3 {
                return Err(RustMeshError::ParseError {
                    line: line_num,
                    reason: format!("Invalid coordinate line (expected ID X Y [Z]): {}", lines[line_num]),
                });
            }

            let x: f64 = parts[1].parse().map_err(|_| RustMeshError::ParseError {
                line: line_num,
                reason: format!("Invalid X coordinate: {}", parts[1]),
            })?;

            let y: f64 = parts[2].parse().map_err(|_| RustMeshError::ParseError {
                line: line_num,
                reason: format!("Invalid Y coordinate: {}", parts[2]),
            })?;

            // Validate coordinates are finite
            if !x.is_finite() || !y.is_finite() {
                return Err(RustMeshError::InvalidGeometry(format!(
                    "Non-finite coordinate at vertex {}: ({}, {})",
                    vert_idx, x, y
                )));
            }

            points[[vert_idx, 0]] = x;
            points[[vert_idx, 1]] = y;
        }

        // Parse element connectivity (after coordinates)
        for elem_idx in 0..n_elems {
            let line_num = 2 + n_verts + elem_idx;
            if line_num >= lines.len() {
                return Err(RustMeshError::ParseError {
                    line: line_num,
                    reason: "Unexpected end of file while parsing elements".to_string(),
                });
            }

            parse_element_line(elem_idx, line_num, &lines[line_num], n_verts, &mut connectivity, &mut elem_type)?;
        }
    } else {
        // Parse element connectivity first (starting at line 2, 0-indexed)
        for elem_idx in 0..n_elems {
            let line_num = 2 + elem_idx;
            if line_num >= lines.len() {
                return Err(RustMeshError::ParseError {
                    line: line_num,
                    reason: "Unexpected end of file while parsing elements".to_string(),
                });
            }

            parse_element_line(elem_idx, line_num, &lines[line_num], n_verts, &mut connectivity, &mut elem_type)?;
        }

        // Parse nodes/coordinates (after elements)
        for vert_idx in 0..n_verts {
            let line_num = 2 + n_elems + vert_idx;
            if line_num >= lines.len() {
                return Err(RustMeshError::ParseError {
                    line: line_num,
                    reason: "Unexpected end of file while parsing coordinates".to_string(),
                });
            }

            let parts: Vec<&str> = lines[line_num].split_whitespace().collect();
            if parts.len() < 3 {
                return Err(RustMeshError::ParseError {
                    line: line_num,
                    reason: format!("Invalid coordinate line (expected ID X Y [Z]): {}", lines[line_num]),
                });
            }

            let x: f64 = parts[1].parse().map_err(|_| RustMeshError::ParseError {
                line: line_num,
                reason: format!("Invalid X coordinate: {}", parts[1]),
            })?;

            let y: f64 = parts[2].parse().map_err(|_| RustMeshError::ParseError {
                line: line_num,
                reason: format!("Invalid Y coordinate: {}", parts[2]),
            })?;

            // Validate coordinates are finite
            if !x.is_finite() || !y.is_finite() {
                return Err(RustMeshError::InvalidGeometry(format!(
                    "Non-finite coordinate at vertex {}: ({}, {})",
                    vert_idx, x, y
                )));
            }

            points[[vert_idx, 0]] = x;
            points[[vert_idx, 1]] = y;
        }
    }

    Ok(RustMesh {
        points,
        connectivity,
        elem_type,
        num_verts: n_verts,
        num_elems: n_elems,
        edges: None,
        areas: None,
        layers: None,
    })
}

/// Helper function to parse an element line
fn parse_element_line(
    elem_idx: usize,
    line_num: usize,
    line: &str,
    n_verts: usize,
    connectivity: &mut Array2<i32>,
    elem_type: &mut Vec<u32>,
) -> Result<(), RustMeshError> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 4 {
        return Err(RustMeshError::ParseError {
            line: line_num,
            reason: format!(
                "Invalid element line (expected IE ITYPE V1 V2 [V3 [V4]]): {}",
                line
            ),
        });
    }

    // Parse element type
    let itype: u32 = parts[1].parse().map_err(|_| RustMeshError::ParseError {
        line: line_num,
        reason: format!("Invalid ITYPE (expected 3 or 4): {}", parts[1]),
    })?;

    if itype != 3 && itype != 4 {
        return Err(RustMeshError::InvalidElement {
            elem_id: elem_idx,
            reason: format!("Invalid element type {} (expected 3 or 4)", itype),
        });
    }

    elem_type.push(itype);

    // Parse vertex indices (1-indexed in file, convert to 0-indexed)
    // Note: Some Fort.14 files have vertex IDs as floats (e.g., "1.000000" instead of "1"),
    // so we parse as f64 first, then convert to i32
    let v0: i32 = (parts[2]
        .parse::<f64>()
        .map_err(|_| RustMeshError::ParseError {
            line: line_num,
            reason: format!("Invalid vertex ID: {}", parts[2]),
        })? as i32)
        - 1;

    let v1: i32 = (parts[3]
        .parse::<f64>()
        .map_err(|_| RustMeshError::ParseError {
            line: line_num,
            reason: format!("Invalid vertex ID: {}", parts[3]),
        })? as i32)
        - 1;

    // Validate vertex bounds for v0, v1
    if v0 < 0 || v0 >= n_verts as i32 {
        return Err(RustMeshError::VertexOutOfBounds {
            vertex_id: v0 + 1,
            max_vertices: n_verts,
        });
    }
    if v1 < 0 || v1 >= n_verts as i32 {
        return Err(RustMeshError::VertexOutOfBounds {
            vertex_id: v1 + 1,
            max_vertices: n_verts,
        });
    }

    connectivity[[elem_idx, 0]] = v0;
    connectivity[[elem_idx, 1]] = v1;

    if itype == 3 {
        // Triangle: parse V3, pad with V3 (repeated vertex)
        if parts.len() < 5 {
            return Err(RustMeshError::ParseError {
                line: line_num,
                reason: "Triangle element missing V3".to_string(),
            });
        }

        let v2: i32 = (parts[4]
            .parse::<f64>()
            .map_err(|_| RustMeshError::ParseError {
                line: line_num,
                reason: format!("Invalid vertex ID: {}", parts[4]),
            })? as i32)
            - 1;

        if v2 < 0 || v2 >= n_verts as i32 {
            return Err(RustMeshError::VertexOutOfBounds {
                vertex_id: v2 + 1,
                max_vertices: n_verts,
            });
        }

        connectivity[[elem_idx, 2]] = v2;
        connectivity[[elem_idx, 3]] = v2; // Padded with V3 per spec
    } else {
        // Quad: parse V3 and V4
        if parts.len() < 6 {
            return Err(RustMeshError::ParseError {
                line: line_num,
                reason: "Quad element missing V3 or V4".to_string(),
            });
        }

        let v2: i32 = (parts[4]
            .parse::<f64>()
            .map_err(|_| RustMeshError::ParseError {
                line: line_num,
                reason: format!("Invalid vertex ID: {}", parts[4]),
            })? as i32)
            - 1;

        let v3: i32 = (parts[5]
            .parse::<f64>()
            .map_err(|_| RustMeshError::ParseError {
                line: line_num,
                reason: format!("Invalid vertex ID: {}", parts[5]),
            })? as i32)
            - 1;

        if v2 < 0 || v2 >= n_verts as i32 {
            return Err(RustMeshError::VertexOutOfBounds {
                vertex_id: v2 + 1,
                max_vertices: n_verts,
            });
        }
        if v3 < 0 || v3 >= n_verts as i32 {
            return Err(RustMeshError::VertexOutOfBounds {
                vertex_id: v3 + 1,
                max_vertices: n_verts,
            });
        }

        connectivity[[elem_idx, 2]] = v2;
        connectivity[[elem_idx, 3]] = v3;
    }

    Ok(())
}

/// Write RustMesh to Fort.14 format
pub fn write_fort14(mesh: &RustMesh, path: &str) -> Result<(), RustMeshError> {
    let mut output = String::new();

    // Write header
    output.push_str(&format!("{} {}\n", mesh.num_elems, mesh.num_verts));

    // Write elements
    for elem_idx in 0..mesh.num_elems {
        let itype = mesh.elem_type[elem_idx];
        let v0 = mesh.connectivity[[elem_idx, 0]] + 1; // Convert to 1-indexed
        let v1 = mesh.connectivity[[elem_idx, 1]] + 1;
        let v2 = mesh.connectivity[[elem_idx, 2]] + 1;

        if itype == 3 {
            // Triangle: write v0, v1, v2 (v3 is redundant)
            output.push_str(&format!("{} {} {} {} {}\n", elem_idx + 1, itype, v0, v1, v2));
        } else {
            // Quad: write all four vertices
            let v3 = mesh.connectivity[[elem_idx, 3]] + 1;
            output.push_str(&format!("{} {} {} {} {} {}\n", elem_idx + 1, itype, v0, v1, v2, v3));
        }
    }

    // Write coordinates
    for vert_idx in 0..mesh.num_verts {
        let x = mesh.points[[vert_idx, 0]];
        let y = mesh.points[[vert_idx, 1]];
        output.push_str(&format!("{:.8e} {:.8e}\n", x, y));
    }

    // Write to file
    fs::write(path, output).map_err(|e| RustMeshError::IOError(e.to_string()))?;
    Ok(())
}
