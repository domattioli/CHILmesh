"""Lazy header-only mesh metadata reading for CHILmesh.

Supports fast metadata extraction from mesh files (fort.14, 2dm, fort.13, fort.15, .msh, .npy, .npz) without
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
    FileNotFoundError
        If the given path does not exist.
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
    # Missing file → FileNotFoundError so CLI exit code matches `info`
    # (cli.main catches FileNotFoundError → exit 2). See issue #235.
    if not path.exists():
        raise FileNotFoundError(f"No such file: {path}")

    # Infer format from file suffix
    suffix = path.suffix.lower()

    if suffix in ('.14', '.grd'):
        fmt = 'fort14'
    elif suffix in ('.fort14',):
        fmt = 'fort14'
    elif suffix == '.2dm':
        fmt = '2dm'
    elif suffix == '.13':
        fmt = 'fort13'
    elif suffix == '.15':
        fmt = 'fort15'
    elif suffix == '.npy':
        fmt = 'npy'
    elif suffix == '.npz':
        fmt = 'npz'
    elif suffix == '.msh':
        fmt = 'gmsh'
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
    elif fmt == 'fort13':
        _read_fort13_header(path, result)
    elif fmt == 'fort15':
        _read_fort15_header(path, result)
    elif fmt == 'npy':
        _read_npy_header(path, result)
    elif fmt == 'npz':
        _read_npz_header(path, result)
    elif fmt == 'gmsh':
        _read_msh_header(path, result)

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


def _read_fort13_header(path: Path, result: dict) -> None:
    """Read fort.13 header (grid name + node count + attribute count).

    fort.13 layout: line 1 = AGRID (grid name), line 2 = NumOfNodes,
    line 3 = NAttr (number of nodal attributes). Read only these 3 lines —
    the per-node non-default overlay blocks (potentially millions of rows)
    are never loaded. Blank lines skipped.
    """
    try:
        vals = []
        with open(path, 'r', encoding='utf-8') as f:
            for line in f:
                s = line.strip()
                if not s:
                    continue
                vals.append(s)
                if len(vals) >= 3:
                    break
        if len(vals) < 3:
            raise SummaryError(f"fort.13 header malformed: expected 3 lines, got {len(vals)}")
        try:
            n_nodes = int(vals[1].split()[0])
            n_attributes = int(vals[2].split()[0])
        except (ValueError, IndexError) as e:
            raise SummaryError(f"fort.13 header parse error: {e}")
        result['grid_name'] = vals[0]
        result['n_nodes'] = n_nodes
        result['n_attributes'] = n_attributes
    except IOError as e:
        raise SummaryError(f"Cannot read fort.13 file {path}: {e}")


def _read_fort15_header(path: Path, result: dict) -> None:
    """Read fort.15 header (run description + run identification).

    fort.15 is ADCIRC model run-config, not a mesh — it carries no node/element
    counts. Line 1 = RUNDES, line 2 = RUNID; trailing inline ``! comment`` is
    stripped (mirrors read_fort15). Only 2 lines are read.
    """
    try:
        with open(path, 'r', encoding='utf-8') as f:
            line1 = f.readline()
            line2 = f.readline()
        if not line1 or not line2:
            raise SummaryError("fort.15 header malformed: fewer than 2 lines")
        result['rundes'] = line1.split('!', 1)[0].rstrip()
        result['runid'] = line2.split('!', 1)[0].rstrip()
    except IOError as e:
        raise SummaryError(f"Cannot read fort.15 file {path}: {e}")


def _read_npy_header(path: Path, result: dict) -> None:
    """Read a .npy header (shape + dtype) without loading the array body.

    Uses numpy.load with mmap_mode='r' to memory-map the file lazily,
    extracting metadata without allocating the full array in RAM.
    """
    import numpy as np
    try:
        # Memory-map the array without loading it into RAM
        mmap = np.load(path, allow_pickle=False, mmap_mode='r')
        result['shape'] = list(mmap.shape)
        result['dtype'] = str(mmap.dtype)
        result['fortran_order'] = mmap.flags['F_CONTIGUOUS']
    except Exception as e:
        raise SummaryError(f"Cannot read .npy header {path}: {e}")


def _read_npz_header(path: Path, result: dict) -> None:
    """Read a .npz archive's member headers (name -> shape/dtype) lazily.

    A .npz is a zip of .npy members; each member's header is extracted via
    numpy.load with mmap_mode, never loading array bodies into RAM.
    """
    import zipfile
    import numpy as np
    try:
        arrays = {}
        with zipfile.ZipFile(path) as z:
            for entry in z.namelist():
                # Remove .npy suffix to get the array name
                name = entry[:-4] if entry.endswith('.npy') else entry
                # Stream only the .npy header off the zip member — the array
                # body is never read into memory (open(), not read()).
                with z.open(entry) as member:
                    major, _minor = np.lib.format.read_magic(member)
                    if major == 2:
                        shape, _fortran_order, dtype = np.lib.format.read_array_header_2_0(member)
                    else:
                        shape, _fortran_order, dtype = np.lib.format.read_array_header_1_0(member)
                arrays[name] = {'shape': list(shape), 'dtype': str(dtype)}
        result['members'] = sorted(arrays.keys())
        result['n_members'] = len(arrays)
        result['arrays'] = arrays
    except Exception as e:
        raise SummaryError(f"Cannot read .npz header {path}: {e}")


def _read_msh_header(path: Path, result: dict) -> None:
    """Read a Gmsh .msh header (version + node/element counts) by streaming.

    Scans the file line-by-line without allocating node/element arrays (same
    streaming, no-array approach as the .2dm path). Supports Gmsh ASCII
    versions 2.2 and 4.1. In v2.2 the count line after ``$Nodes`` /
    ``$Elements`` is a single integer; in v4.1 it is
    ``numBlocks numEntities minTag maxTag`` and the entity count is the second
    token. The scan stops once the element count is read, so element bodies are
    never traversed.
    """
    version = None
    try:
        with open(path, 'r', encoding='utf-8') as f:
            for line in f:
                s = line.strip()
                if s == '$MeshFormat':
                    fmt_line = f.readline().split()
                    if fmt_line:
                        version = fmt_line[0]
                        result['gmsh_version'] = version
                elif s == '$Nodes':
                    toks = f.readline().split()
                    if toks:
                        idx = 1 if (version and version.startswith('4')) else 0
                        result['n_nodes'] = int(toks[idx])
                elif s == '$Elements':
                    toks = f.readline().split()
                    if toks:
                        idx = 1 if (version and version.startswith('4')) else 0
                        result['n_elems'] = int(toks[idx])
                    break
    except (IOError, ValueError, IndexError) as e:
        raise SummaryError(f"Cannot read .msh header {path}: {e}")
    if version is None:
        raise SummaryError(f"gmsh .msh header malformed: no $MeshFormat in {path}")


__all__ = ["summary", "SummaryError"]
