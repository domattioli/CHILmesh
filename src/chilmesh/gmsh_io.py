"""Gmsh ASCII format I/O (.msh) for CHILmesh.

Supports Gmsh format versions 2.2 and 4.1 with triangular and quadrilateral elements.
"""
from __future__ import annotations

import numpy as np
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .CHILmesh import CHILmesh


class GmshParseError(Exception):
    """Raised when parsing a Gmsh file encounters an error."""
    pass


def _require_int(token: str, what: str) -> int:
    """Parse and validate an integer token.

    Parameters:
        token: String to parse as integer
        what: Human-readable description of the field (for error message)

    Returns:
        Parsed integer value

    Raises:
        GmshParseError: If token cannot be parsed as integer
    """
    try:
        return int(token)
    except ValueError as exc:
        raise GmshParseError(f"malformed {what}: non-numeric {token!r}") from exc


def read_msh(full_file_name: str) -> "CHILmesh":
    """
    Load a mesh from a Gmsh ASCII .msh file.

    Supports both format 2.2 and 4.1. Parses nodes and elements, supporting
    triangular (type 2) and quadrilateral (type 3) elements. Other element types
    are silently skipped. Triangles in a mixed-element mesh use the padded
    convention [v0, v1, v2, v0].

    Parameters:
        full_file_name: Path to the .msh file

    Returns:
        A CHILmesh object with connectivity and points, fast-initialized
        (compute_layers=False, compute_adjacencies=False).

    Raises:
        GmshParseError: If format is unsupported, required sections missing,
            or file is malformed.
    """

    with open(full_file_name, 'r', encoding='utf-8') as f:
        lines = [line.strip() for line in f]

    # Find $MeshFormat section
    mf_idx = None
    for i, line in enumerate(lines):
        if line == "$MeshFormat":
            mf_idx = i
            break

    if mf_idx is None:
        raise GmshParseError("Missing $MeshFormat section in .msh file")

    # Parse version from the first line after $MeshFormat
    if mf_idx + 1 >= len(lines):
        raise GmshParseError("$MeshFormat section incomplete")

    format_line = lines[mf_idx + 1].split()
    if not format_line:
        raise GmshParseError("$MeshFormat section has no version line")

    version_str = format_line[0]
    if version_str not in ("2.2", "2", "4.1", "4"):
        raise GmshParseError(f"Unsupported Gmsh format version: {version_str}")

    version = "2.2" if version_str in ("2.2", "2") else "4.1"

    if version == "2.2":
        return _read_msh_v2_2(lines, full_file_name)
    else:
        return _read_msh_v4_1(lines, full_file_name)


def _find_node_and_element_sections(lines: list) -> tuple:
    """Return the line indices of the last ``$Nodes`` and ``$Elements`` headers.

    The last occurrence wins, matching the original single-pass scan.

    Raises
    ------
    GmshParseError
        If either section header is absent.
    """
    nodes_idx = None
    elems_idx = None
    for i, line in enumerate(lines):
        if line == "$Nodes":
            nodes_idx = i
        elif line == "$Elements":
            elems_idx = i

    if nodes_idx is None:
        raise GmshParseError("Missing $Nodes section in .msh file")
    if elems_idx is None:
        raise GmshParseError("Missing $Elements section in .msh file")
    return nodes_idx, elems_idx


def _parse_v2_count(lines: list, idx: int, section: str) -> int:
    """Read the entry count that follows a v2.2 section header at ``idx``."""
    if idx + 1 >= len(lines):
        raise GmshParseError(f"{section} section incomplete")
    try:
        return int(lines[idx + 1])
    except ValueError:
        raise GmshParseError(f"{section} section count is not an integer")


def _parse_v2_node_line(line: str) -> tuple:
    """Parse one v2.2 node line into ``(node_id, [x, y, z])``."""
    parts = line.split()
    if len(parts) < 4:
        raise GmshParseError(f"Node line malformed: {line}")
    try:
        node_id = int(parts[0])
        return node_id, [float(parts[1]), float(parts[2]), float(parts[3])]
    except ValueError:
        raise GmshParseError(f"Node line has non-numeric values: {line}")


def _parse_v2_nodes(lines: list, nodes_idx: int) -> dict:
    """Parse the v2.2 ``$Nodes`` section into ``{node_id: [x, y, z]}``."""
    n_nodes = _parse_v2_count(lines, nodes_idx, "$Nodes")
    nodes = {}
    for i in range(n_nodes):
        line_idx = nodes_idx + 2 + i
        if line_idx >= len(lines):
            raise GmshParseError(f"$Nodes section incomplete: expected {n_nodes} nodes")
        node_id, coords = _parse_v2_node_line(lines[line_idx])
        nodes[node_id] = coords
    return nodes


def _parse_v2_element_line(line: str):
    """Parse one v2.2 element line.

    Returns
    -------
    tuple or None
        ``("tri", [v1, v2, v3])`` or ``("quad", [v1, v2, v3, v4])``; ``None``
        for unsupported element types, which are skipped.
    """
    parts = line.split()
    if len(parts) < 4:
        raise GmshParseError(f"Element line malformed: {line}")
    try:
        _require_int(parts[0], "element ID")
        elem_type = int(parts[1])
        n_tags = int(parts[2])
        # Skip tag section, extract node list
        node_start = 3 + n_tags
        if len(parts) < node_start + (3 if elem_type == 2 else 4 if elem_type == 3 else 0):
            raise GmshParseError(f"Element line too short: {line}")
        # Index (not slice) so a negative node_start from a negative tag count
        # behaves exactly as before.
        if elem_type == 2:  # Triangle
            return ("tri", [int(parts[node_start + k]) for k in range(3)])
        if elem_type == 3:  # Quad
            return ("quad", [int(parts[node_start + k]) for k in range(4)])
        return None  # Other types are silently skipped
    except (ValueError, IndexError):
        raise GmshParseError(f"Element line malformed or non-numeric: {line}")


def _parse_v2_elements(lines: list, elems_idx: int) -> list:
    """Parse the v2.2 ``$Elements`` section into ``(kind, vertices)`` tuples."""
    n_elems = _parse_v2_count(lines, elems_idx, "$Elements")
    elements = []
    elem_idx = elems_idx + 2
    for i in range(n_elems):
        if elem_idx >= len(lines):
            raise GmshParseError(f"$Elements section incomplete: expected {n_elems} elements")
        parsed = _parse_v2_element_line(lines[elem_idx])
        if parsed is not None:
            elements.append(parsed)
        elem_idx += 1
    return elements


def _read_msh_v2_2(lines: list, filename: str) -> "CHILmesh":
    """Parse Gmsh format 2.2."""
    nodes_idx, elems_idx = _find_node_and_element_sections(lines)

    nodes = _parse_v2_nodes(lines, nodes_idx)
    if not nodes:
        raise GmshParseError("No nodes found in .msh file")

    elements = _parse_v2_elements(lines, elems_idx)
    if not elements:
        raise GmshParseError("No supported elements (triangles or quads) found in .msh file")

    return _build_mesh(nodes, elements)


def _parse_v41_section_header(lines: list, idx: int, section: str, total_label: str) -> int:
    """Read a v4.1 ``$Nodes``/``$Elements`` header and return the block count.

    The total entity count (second field) is validated but not used.
    """
    if idx + 1 >= len(lines):
        raise GmshParseError(f"{section} section incomplete")
    try:
        header_parts = lines[idx + 1].split()
        if len(header_parts) < 4:
            raise GmshParseError(f"{section} header malformed")
        num_blocks = int(header_parts[0])
        _require_int(header_parts[1], total_label)
    except ValueError:
        raise GmshParseError(f"{section} header has non-numeric values")
    return num_blocks


def _parse_v41_node_block_header(lines: list, line_idx: int, block_idx: int) -> int:
    """Return the node count of the v4.1 node block header at ``line_idx``."""
    if line_idx >= len(lines):
        raise GmshParseError("$Nodes section incomplete")
    block_header = lines[line_idx].split()
    if len(block_header) < 4:
        raise GmshParseError(f"Node block {block_idx} header malformed")
    try:
        _require_int(block_header[0], f"node block {block_idx} dimension")
        _require_int(block_header[1], f"node block {block_idx} tag")
        _require_int(block_header[2], f"node block {block_idx} parametric")
        return int(block_header[3])
    except ValueError:
        raise GmshParseError(f"Node block {block_idx} header has non-numeric values")


def _parse_v41_node_tags(lines: list, line_idx: int, block_idx: int, count: int) -> list:
    """Read ``count`` node-tag lines starting at ``line_idx``."""
    node_tags = []
    for i in range(count):
        if line_idx + i >= len(lines):
            raise GmshParseError(f"Node block {block_idx} node-tags incomplete")
        try:
            node_tags.append(int(lines[line_idx + i]))
        except ValueError:
            raise GmshParseError(f"Node tag line malformed: {lines[line_idx + i]}")
    return node_tags


def _parse_v41_node_coords(lines: list, line_idx: int, block_idx: int, node_tags: list) -> list:
    """Read one coordinate line per tag; returns ``[(node_id, [x, y, z]), ...]``."""
    out = []
    for i, node_id in enumerate(node_tags):
        if line_idx + i >= len(lines):
            raise GmshParseError(f"Node block {block_idx} coordinates incomplete")
        line = lines[line_idx + i]
        parts = line.split()
        if len(parts) < 3:
            raise GmshParseError(f"Node coordinate line malformed: {line}")
        try:
            out.append((node_id, [float(parts[0]), float(parts[1]), float(parts[2])]))
        except ValueError:
            raise GmshParseError(f"Node coordinate line has non-numeric values: {line}")
    return out


def _parse_v41_nodes(lines: list, nodes_idx: int) -> dict:
    """Parse the v4.1 ``$Nodes`` section into ``{node_id: [x, y, z]}``."""
    num_blocks = _parse_v41_section_header(lines, nodes_idx, "$Nodes", "number of nodes")
    nodes = {}
    line_idx = nodes_idx + 2
    for block_idx in range(num_blocks):
        num_in_block = _parse_v41_node_block_header(lines, line_idx, block_idx)
        line_idx += 1
        node_tags = _parse_v41_node_tags(lines, line_idx, block_idx, num_in_block)
        # len(node_tags), not num_in_block: a negative count must not move the cursor back.
        line_idx += len(node_tags)
        for node_id, coords in _parse_v41_node_coords(lines, line_idx, block_idx, node_tags):
            nodes[node_id] = coords
        line_idx += len(node_tags)
    return nodes


def _parse_v41_element_block_header(lines: list, line_idx: int, block_idx: int) -> tuple:
    """Return ``(elem_type, num_in_block)`` for the v4.1 element block header."""
    if line_idx >= len(lines):
        raise GmshParseError("$Elements section incomplete")
    block_header = lines[line_idx].split()
    if len(block_header) < 4:
        raise GmshParseError(f"Element block {block_idx} header malformed")
    try:
        _require_int(block_header[0], f"element block {block_idx} dimension")
        _require_int(block_header[1], f"element block {block_idx} tag")
        return int(block_header[2]), int(block_header[3])
    except ValueError:
        raise GmshParseError(f"Element block {block_idx} header has non-numeric values")


def _parse_v41_element_line(line: str, elem_type: int):
    """Parse one v4.1 element line; ``None`` for unsupported types (skipped)."""
    parts = line.split()
    try:
        int(parts[0])  # validation only: ValueError -> enclosing except keeps the "Element line malformed" contract
        if elem_type == 2:  # Triangle
            if len(parts) < 4:
                raise GmshParseError(f"Triangle element line too short: {line}")
            return ("tri", [int(p) for p in parts[1:4]])
        if elem_type == 3:  # Quad
            if len(parts) < 5:
                raise GmshParseError(f"Quad element line too short: {line}")
            return ("quad", [int(p) for p in parts[1:5]])
        return None  # Other types silently skipped
    except (ValueError, IndexError):
        raise GmshParseError(f"Element line malformed: {line}")


def _parse_v41_elements(lines: list, elems_idx: int) -> list:
    """Parse the v4.1 ``$Elements`` section into ``(kind, vertices)`` tuples."""
    elem_blocks = _parse_v41_section_header(
        lines, elems_idx, "$Elements", "total number of elements"
    )
    elements = []
    line_idx = elems_idx + 2
    for block_idx in range(elem_blocks):
        elem_type, num_in_block = _parse_v41_element_block_header(lines, line_idx, block_idx)
        line_idx += 1
        for i in range(num_in_block):
            if line_idx >= len(lines):
                raise GmshParseError(f"Element block {block_idx} incomplete")
            parsed = _parse_v41_element_line(lines[line_idx], elem_type)
            if parsed is not None:
                elements.append(parsed)
            line_idx += 1
    return elements


def _read_msh_v4_1(lines: list, filename: str) -> "CHILmesh":
    """Parse Gmsh format 4.1."""
    nodes_idx, elems_idx = _find_node_and_element_sections(lines)

    nodes = _parse_v41_nodes(lines, nodes_idx)
    if not nodes:
        raise GmshParseError("No nodes found in .msh file")

    elements = _parse_v41_elements(lines, elems_idx)
    if not elements:
        raise GmshParseError("No supported elements (triangles or quads) found in .msh file")

    return _build_mesh(nodes, elements)


def _build_mesh(nodes: dict, elements: list) -> "CHILmesh":
    """Build CHILmesh from parsed nodes and elements.

    Parameters:
        nodes: Dict mapping node_id to [x, y, z]
        elements: List of (elem_type, vertex_list) tuples

    Returns:
        CHILmesh object
    """
    from .CHILmesh import CHILmesh

    # Sort node IDs and build index map
    node_ids = sorted(nodes.keys())
    points = np.array([nodes[nid] for nid in node_ids], dtype=float)
    id_to_idx = {nid: idx for idx, nid in enumerate(node_ids)}

    # Determine element array width
    has_quads = any(elem_type == "quad" for elem_type, _ in elements)
    elem_width = 4 if has_quads else 3

    # Build connectivity array
    connectivity = np.zeros((len(elements), elem_width), dtype=int)
    for i, (elem_type, vertices) in enumerate(elements):
        if elem_type == "tri":
            v1, v2, v3 = vertices
            idx1, idx2, idx3 = id_to_idx[v1], id_to_idx[v2], id_to_idx[v3]
            if elem_width == 4:
                # Mixed mesh: use padding convention [v0, v1, v2, v0]
                connectivity[i] = [idx1, idx2, idx3, idx1]
            else:
                connectivity[i] = [idx1, idx2, idx3]
        else:  # quad
            v1, v2, v3, v4 = vertices
            connectivity[i] = [id_to_idx[v1], id_to_idx[v2], id_to_idx[v3], id_to_idx[v4]]

    return CHILmesh(
        connectivity=connectivity,
        points=points,
        grid_name="Gmsh",
        compute_layers=False,
        compute_adjacencies=False,
    )


def write_msh(
    filename: str,
    points: np.ndarray,
    connectivity: np.ndarray,
    grid_name: str = "CHILmesh Grid",
    version: str = "4.1",
) -> bool:
    """
    Write mesh data to a Gmsh ASCII .msh file.

    Supports triangular, quadrilateral, and mixed-element meshes. Triangles in a
    4-column connectivity array using the padded convention [v0, v1, v2, v0] are
    written as 3-node triangles. All node IDs and element indices are written
    1-based (Gmsh convention).

    Parameters:
        filename: Output path
        points: (n_nodes, 2 or 3) numpy array of node coordinates
        connectivity: (n_elems, 3 or 4) array of vertex indices (0-based)
        grid_name: Optional title/name for the mesh
        version: Gmsh format version, "2.2" or "4.1" (default: "4.1")

    Returns:
        True on success.

    Raises:
        ValueError: If version is not "2.2" or "4.1".
    """
    if version not in ("2.2", "4.1"):
        raise ValueError(f"Unsupported Gmsh version: {version}. Use '2.2' or '4.1'.")

    if version == "2.2":
        return _write_msh_v2_2(filename, points, connectivity, grid_name)
    else:
        return _write_msh_v4_1(filename, points, connectivity, grid_name)


def _write_msh_v2_2(
    filename: str,
    points: np.ndarray,
    connectivity: np.ndarray,
    grid_name: str,
) -> bool:
    """Write Gmsh format 2.2."""
    try:
        with open(filename, 'w', encoding='utf-8', newline='\n') as f:
            f.write("$MeshFormat\n")
            f.write("2.2 0 8\n")
            f.write("$EndMeshFormat\n")

            # Write nodes
            f.write("$Nodes\n")
            f.write(f"{len(points)}\n")
            for i, pt in enumerate(points, start=1):
                x, y = pt[:2]
                z = pt[2] if len(pt) >= 3 else 0.0
                f.write(f"{i} {x:.8e} {y:.8e} {z:.8e}\n")
            f.write("$EndNodes\n")

            # Write elements
            f.write("$Elements\n")
            f.write(f"{len(connectivity)}\n")
            for i, elem in enumerate(connectivity, start=1):
                if len(elem) == 3:
                    # 3-column array: all triangles
                    n1, n2, n3 = elem + 1
                    f.write(f"{i} 2 2 0 0 {n1} {n2} {n3}\n")
                else:
                    # 4-column array: check if padded triangle or quad
                    if elem[3] == elem[0]:
                        # Padded triangle: [v0, v1, v2, v0]
                        n1, n2, n3 = (elem[:3] + 1)
                        f.write(f"{i} 2 2 0 0 {n1} {n2} {n3}\n")
                    else:
                        # Quad: [v0, v1, v2, v3]
                        n1, n2, n3, n4 = elem + 1
                        f.write(f"{i} 3 2 0 0 {n1} {n2} {n3} {n4}\n")
            f.write("$EndElements\n")

        return True
    except Exception as e:
        print(f"Error writing Gmsh file {filename}: {e}")
        return False


def _write_msh_v4_1(
    filename: str,
    points: np.ndarray,
    connectivity: np.ndarray,
    grid_name: str,
) -> bool:
    """Write Gmsh format 4.1."""
    try:
        with open(filename, 'w', encoding='utf-8', newline='\n') as f:
            f.write("$MeshFormat\n")
            f.write("4.1 0 8\n")
            f.write("$EndMeshFormat\n")

            # Write nodes in a single block
            f.write("$Nodes\n")
            n_nodes = len(points)
            f.write(f"1 {n_nodes} 1 {n_nodes}\n")  # 1 block, n_nodes nodes, min_tag=1, max_tag=n_nodes
            f.write("0 1 0 {}\n".format(n_nodes))  # Node block header: dim=0, tag=1, parametric=0, count

            # Node tags (1-based)
            for i in range(1, n_nodes + 1):
                f.write(f"{i}\n")

            # Node coordinates
            for pt in points:
                x, y = pt[:2]
                z = pt[2] if len(pt) >= 3 else 0.0
                f.write(f"{x:.8e} {y:.8e} {z:.8e}\n")

            # Write elements
            f.write("$Elements\n")

            # Categorize elements by type
            tri_elems = []
            quad_elems = []
            for i, elem in enumerate(connectivity, start=1):
                if len(elem) == 3:
                    # Pure triangle
                    tri_elems.append((i, elem))
                else:
                    # Check if padded triangle or quad
                    if elem[3] == elem[0]:
                        # Padded triangle
                        tri_elems.append((i, elem[:3]))
                    else:
                        # Quad
                        quad_elems.append((i, elem))

            # Count blocks and total elements
            n_blocks = (1 if tri_elems else 0) + (1 if quad_elems else 0)
            total_elems = len(tri_elems) + len(quad_elems)
            f.write(f"{n_blocks} {total_elems} 1 {total_elems}\n")

            elem_tag = 1
            # Write triangle block if present
            if tri_elems:
                f.write(f"2 1 2 {len(tri_elems)}\n")  # dim=2, tag=1, elemType=2 (tri), count
                for i, elem in tri_elems:
                    n1, n2, n3 = elem + 1
                    f.write(f"{elem_tag} {n1} {n2} {n3}\n")
                    elem_tag += 1

            # Write quad block if present
            if quad_elems:
                f.write(f"2 2 3 {len(quad_elems)}\n")  # dim=2, tag=2, elemType=3 (quad), count
                for i, elem in quad_elems:
                    n1, n2, n3, n4 = elem + 1
                    f.write(f"{elem_tag} {n1} {n2} {n3} {n4}\n")
                    elem_tag += 1

            f.write("$EndElements\n")

        return True
    except Exception as e:
        print(f"Error writing Gmsh file {filename}: {e}")
        return False
