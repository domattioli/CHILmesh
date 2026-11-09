"""Synthetic tri, quad and mixed grids shared by the benchmark-harness tests."""
from __future__ import annotations

import numpy as np


def _grid_points(nx: int, ny: int) -> np.ndarray:
    xs, ys = np.meshgrid(np.arange(nx + 1, dtype=float), np.arange(ny + 1, dtype=float))
    return np.column_stack([xs.ravel(), ys.ravel()])


def _node(i: int, j: int, nx: int) -> int:
    return j * (nx + 1) + i


def tri_grid(nx: int = 4, ny: int = 4):
    pts = _grid_points(nx, ny)
    conn = []
    for j in range(ny):
        for i in range(nx):
            a, b = _node(i, j, nx), _node(i + 1, j, nx)
            c, d = _node(i + 1, j + 1, nx), _node(i, j + 1, nx)
            conn.append([a, b, c])
            conn.append([a, c, d])
    return pts, np.array(conn, dtype=int)


def quad_grid(nx: int = 4, ny: int = 4):
    pts = _grid_points(nx, ny)
    conn = []
    for j in range(ny):
        for i in range(nx):
            conn.append([_node(i, j, nx), _node(i + 1, j, nx),
                         _node(i + 1, j + 1, nx), _node(i, j + 1, nx)])
    return pts, np.array(conn, dtype=int)


def mixed_grid(nx: int = 4, ny: int = 4):
    """Left half quads, right half triangles. Triangles are padded with the first vertex."""
    pts = _grid_points(nx, ny)
    conn = []
    for j in range(ny):
        for i in range(nx):
            a, b = _node(i, j, nx), _node(i + 1, j, nx)
            c, d = _node(i + 1, j + 1, nx), _node(i, j + 1, nx)
            if i < nx // 2:
                conn.append([a, b, c, d])
            else:
                conn.append([a, b, c, a])
                conn.append([a, c, d, a])
    return pts, np.array(conn, dtype=int)


GRIDS = {"tri": tri_grid, "quad": quad_grid, "mixed": mixed_grid}
