"""Element coloring over node-sharing adjacency (spec 005 US2)."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["ElementColoring", "color_elements", "check_coloring"]


@dataclass
class ElementColoring:
    """Result of :func:`color_elements`: ``colors[i]`` is element ``i``'s color."""

    colors: np.ndarray
    n_colors: int


def _element_vertices(mesh, e: int) -> set:
    row = mesh.connectivity_list[e]
    if len(row) == 4 and row[3] == row[0]:
        row = row[:3]
    return {int(v) for v in row}


def _neighbors(mesh) -> list:
    v2e = mesh.adjacencies["Vert2Elem"]
    nbrs = []
    for e in range(mesh.n_elems):
        s: set = set()
        for v in _element_vertices(mesh, e):
            s |= v2e[v]
        s.discard(e)
        nbrs.append(s)
    return nbrs


def color_elements(mesh) -> ElementColoring:
    """Greedy-color elements so no two sharing a node get the same color.

    Elements are visited largest node-sharing degree first, ties by element
    id, each taking the smallest free color. Deterministic; uses at most
    ``max_degree + 1`` colors.

    Example:
        >>> from chilmesh import examples
        >>> mesh = examples.donut()
        >>> c = color_elements(mesh)
        >>> check_coloring(mesh, c)
        True

    Tested by: tests/test_coloring.py::test_mixed_mesh, ::test_deterministic, ::test_all_bundled_meshes.
    """
    n = mesh.n_elems
    if "Vert2Elem" not in mesh.adjacencies:
        mesh._build_adjacencies()
    nbrs = _neighbors(mesh)
    order = sorted(range(n), key=lambda e: (-len(nbrs[e]), e))
    colors = np.full(n, -1, dtype=int)
    for e in order:
        used = {int(colors[m]) for m in nbrs[e] if colors[m] >= 0}
        c = 0
        while c in used:
            c += 1
        colors[e] = c
    return ElementColoring(colors=colors, n_colors=int(colors.max()) + 1 if n else 0)


def check_coloring(mesh, colors) -> bool:
    """Return True when no two node-sharing elements share a color.

    Example:
        >>> from chilmesh import examples
        >>> mesh = examples.donut()
        >>> check_coloring(mesh, color_elements(mesh))
        True
        >>> check_coloring(mesh, [0] * mesh.n_elems)
        False

    Tested by: tests/test_coloring.py::test_check_coloring_detects_bad.
    """
    c = colors.colors if isinstance(colors, ElementColoring) else np.asarray(colors)
    if len(c) != mesh.n_elems:
        return False
    if "Vert2Elem" not in mesh.adjacencies:
        mesh._build_adjacencies()
    for e, nb in enumerate(_neighbors(mesh)):
        if any(c[m] == c[e] for m in nb):
            return False
    return True
