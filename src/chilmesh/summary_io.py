"""Lazy header-only mesh metadata reading for CHILmesh.

Supports fast metadata extraction from mesh files (fort.14, 2dm) without
loading full mesh data. For CHILmesh objects, returns mesh properties directly.
"""
from __future__ import annotations

from pathlib import Path


class SummaryError(Exception):
    """Raised when summary extraction encounters an error."""
    pass


def summary(path_or_mesh, *, deep: bool = False) -> dict:
    """Extract metadata from a mesh file or CHILmesh object.

    Performs lazy header-only reading for file paths (fort.14, 2dm).
    For CHILmesh objects, returns object properties; optionally computes
    bounding box and element type via full load when deep=True.

    Parameters
    ----------
    path_or_mesh : str, Path, or CHILmesh
        Either a filesystem path (fort.14, 2dm, .grd) or a CHILmesh instance.
    deep : bool, optional
        If True and path_or_mesh is a file path, load the full mesh to extract
        element_type and bbox. If False (default), skip these for faster
        performance. Ignored for CHILmesh objects.

    Returns
    -------
    dict
        Metadata dictionary. Always includes:
        - n_nodes: number of vertices
        - n_elems: number of elements
        - format: "fort14", "2dm", or "mesh-object"
        - path: str file path, or None for mesh objects

        When deep=True or for mesh objects:
        - element_type: "Triangular", "Quadrilateral", or "Mixed-Element"
        - bbox: dict with keys min_x, min_y, max_x, max_y

        For fort.14 files only:
        - grid_name: mesh title line

    Raises
    ------
    SummaryError
        If file format is unrecognized, file cannot be read, or mesh object
        is invalid.
    """
    # Check if it's a CHILmesh object (duck-typed to avoid circular import)
    if type(path_or_mesh).__name__ == 'CHILmesh':
        return _summary_from_mesh(path_or_mesh)

    # Otherwise treat as a path
    path = Path(path_or_mesh)
    return _summary_from_file(path, deep=deep)


def _summary_from_mesh(mesh) -> dict:
    """Extract metadata from a CHILmesh object."""
    try:
        # Compute bounding box
        points = mesh.points
        bbox = {
            'min_x': float(points[:, 0].min()),
            'max_x': float(points[:, 0].max()),
            'min_y': float(points[:, 1].min()),
            'max_y': float(points[:, 1].max()),
        }

        return {
            'n_nodes': mesh.n_verts,
            'n_elems': mesh.n_elems,
            'element_type': mesh.type,
            'bbox': bbox,
            'path': None,
            'format': 'mesh-object',
        }
    except (AttributeError, IndexError, ValueError) as e:
        raise SummaryError(f"Failed to extract metadata from CHILmesh object: {e}")


def _summary_from_file(path: Path, *, deep: bool = False) -> dict:
    """Extract metadata from a mesh file."""
    # Infer format from file suffix
    suffix = path.suffix.lower()

    if suffix in ('.14', '.grd'):
        fmt = 'fort14'
    elif suffix in ('.fort14',):
        fmt = 'fort14'
    elif suffix == '.2dm':
        fmt = '2dm'
    else:
        raise SummaryError(f"Unknown mesh format: {suffix}")

    # Get file size
    try:
        file_bytes = path.stat().st_size
    except (OSError, ValueError) as e:
        raise SummaryError(f"Cannot stat file {path}: {e}")

    result = {
        'path': str(path),
        'format': fmt,
        'file_bytes': file_bytes,
    }

    # Format-specific lazy reading
    if fmt == 'fort14':
        _read_fort14_header(path, result)
    elif fmt == '2dm':
        _read_2dm_header(path, result)

    # If deep=True, load the full mesh for element_type and bbox
    if deep:
        try:
            # Import here to avoid circular dependency
            from chilmesh import CHILmesh
            mesh = CHILmesh.load(path, compute_layers=False, compute_adjacencies=False)
            result['element_type'] = mesh.type
            points = mesh.points
            result['bbox'] = {
                'min_x': float(points[:, 0].min()),
                'max_x': float(points[:, 0].max()),
                'min_y': float(points[:, 1].min()),
                'max_y': float(points[:, 1].max()),
            }
        except Exception as e:
            raise SummaryError(f"Failed to load mesh for deep summary: {e}")

    return result


def _read_fort14_header(path: Path, result: dict) -> None:
    """Read fort.14 file header (2-line lazy read)."""
    try:
        with open(path, 'r', encoding='utf-8') as f:
            grid_name = f.readline().strip()
            header_line = f.readline().strip()

        # Parse second line: "n_elems n_nodes" format
        parts = header_line.split()
        if len(parts) < 2:
            raise SummaryError(f"fort.14 header malformed: {header_line}")

        try:
            n_elems = int(parts[0])
            n_nodes = int(parts[1])
        except ValueError as e:
            raise SummaryError(f"fort.14 header parse error: {e}")

        result['grid_name'] = grid_name
        result['n_elems'] = n_elems
        result['n_nodes'] = n_nodes
    except IOError as e:
        raise SummaryError(f"Cannot read fort.14 file {path}: {e}")


def _read_2dm_header(path: Path, result: dict) -> None:
    """Read 2dm file header (streaming line-by-line scan)."""
    try:
        n_nodes = 0
        n_elems = 0
        grid_name = None

        with open(path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                # Check for MESHNAME
                if line.startswith('MESHNAME'):
                    parts = line.split(None, 1)
                    if len(parts) > 1:
                        grid_name = parts[1]

                # Count ND (node) lines
                elif line.startswith('ND '):
                    n_nodes += 1

                # Count E3T (triangles) and E4Q (quads)
                elif line.startswith('E3T ') or line.startswith('E4Q '):
                    n_elems += 1

        result['n_nodes'] = n_nodes
        result['n_elems'] = n_elems
        if grid_name is not None:
            result['grid_name'] = grid_name
    except IOError as e:
        raise SummaryError(f"Cannot read 2dm file {path}: {e}")




__all__ = ["summary", "SummaryError"]
