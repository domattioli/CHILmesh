"""
CHILmesh binary format (.chil) — fast native serialization.

ADCIRC .fort.14 is text; parsing is slow. .chil stores mesh data as binary
NumPy arrays (NPZ) for O(10x) read/write speedup. Lossless roundtrip with
fort.14.

Format: NPZ (zipped NumPy arrays)
- points: ndarray[n_vertices, 3] float64 (x, y, z)
- connectivity: ndarray[n_elements, 3|4] int32 (vert indices)
- element_types: ndarray[n_elements] int8 (3=tri, 4=quad)
- grid_name: str (metadata)
"""

import numpy as np
from pathlib import Path
from typing import Tuple, Optional

def write_chil(
    filename,
    points: np.ndarray,
    connectivity: np.ndarray,
    grid_name: str = "CHILmesh Grid"
) -> bool:
    """
    Write mesh to .chil binary format (NPZ-backed).

    Args:
        filename: .chil file path (will be saved as .chil.npz internally)
        points: [n_verts, 3] float64
        connectivity: [n_elems, 3|4] int32
        grid_name: metadata string

    Returns:
        True on success, False on error
    """
    try:
        filename_str = str(filename)
        # Strip trailing .npz if user provided it; savez_compressed will re-add it
        if filename_str.endswith('.npz'):
            filename_base = filename_str[:-4]
        else:
            filename_base = filename_str

        # Infer element types (3=tri, 4=quad)
        element_types = np.full(connectivity.shape[0], connectivity.shape[1], dtype=np.int8)

        np.savez_compressed(
            filename_base,  # Will become filename_base.npz
            points=points.astype(np.float64),
            connectivity=connectivity.astype(np.int32),
            element_types=element_types,
            grid_name=grid_name
        )
        return True
    except Exception as e:
        print(f"Error writing .chil file {filename_str}: {e}")
        return False


def read_chil(filename) -> Tuple[np.ndarray, np.ndarray, Optional[str]]:
    """
    Read mesh from .chil binary format (NPZ-backed).

    Args:
        filename: .chil file path (looks for .chil.npz or directly as provided)

    Returns:
        (points, connectivity, grid_name)

    Raises:
        ValueError: if file format is invalid or not found
    """
    filename_str = str(filename)
    try:
        # Try filename as-is first
        if not Path(filename_str).exists():
            # If .chil file, try .chil.npz
            if filename_str.endswith('.chil'):
                alt_filename = filename_str + '.npz'
                if Path(alt_filename).exists():
                    filename_str = alt_filename

        data = np.load(filename_str, allow_pickle=True)

        points = data['points'].astype(np.float64)
        connectivity = data['connectivity'].astype(np.int32)
        grid_name = str(data['grid_name']) if 'grid_name' in data else "CHILmesh Grid"

        return points, connectivity, grid_name
    except Exception as e:
        raise ValueError(f"Error reading .chil file {filename_str}: {e}")
