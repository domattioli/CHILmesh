"""ADCIRC non-periodic boundary forcing: fort.19, fort.20, and fort.24 (spec 005 US8).

fort.19 (elevation) and fort.20 (normal flux) hold a time increment followed by
one value per matching boundary node per time step, one value per line; text
after the first number on a line (for example ``! comment``) is ignored.
fort.24 (self-attraction and earth load tide) holds, per constituent, four
header lines and one ``node amplitude phase`` line per mesh node.

Matching boundary nodes, taken from ``mesh.boundary_segments`` in file order
(node order matters because the files carry no node ids):

* fort.19: all ``open`` (elevation-specified) segments.
* fort.20: ``flow`` segments whose IBTYPE is in :data:`FORT20_IBTYPES`.
* fort.24: all mesh nodes.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

__all__ = [
    "BoundaryForcingSeries", "TidalConstituent", "FORT20_IBTYPES",
    "read_fort19", "write_fort19", "read_fort20", "write_fort20",
    "read_fort24", "write_fort24",
]

#: Non-periodic normal-flow IBTYPE values that read fort.20. Source: ADCIRC
#: docs (fort20.rst: "IBTYPE=2, 12, or 22"); see research.md R8.
FORT20_IBTYPES = (2, 12, 22)


@dataclass
class BoundaryForcingSeries:
    """Time increment (seconds) and values of a fort.19 or fort.20 file.

    ``values`` is ``(nt, nbn)`` when the per-step node count is known (from a
    mesh or ``n_nodes``). With neither, the file cannot be split into steps and
    ``values`` is the flat 1-D series in file order.
    """

    time_increment: float
    values: np.ndarray


@dataclass
class TidalConstituent:
    """One fort.24 constituent: header text, frequency, and per-node amplitude/phase.

    ``nodal_factor`` is the numeric third header line (always 1 in practice; the
    ADCIRC reader skips it). ``phase`` is in degrees.
    """

    name: str
    alpha: str
    frequency: float
    nodal_factor: float
    amp: np.ndarray = field(repr=False)
    phase: np.ndarray = field(repr=False)


def _float(tok: str) -> float:
    return float(tok.replace("D", "E").replace("d", "e"))


def _first_numbers(path) -> List[float]:
    out = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for lineno, line in enumerate(f, 1):
            tok = line.split()
            if not tok:
                continue
            try:
                out.append(_float(tok[0].rstrip(",")))
            except ValueError:
                raise ValueError(f"{path}:{lineno}: expected a number, got {tok[0]!r}") from None
    return out


def _open_count(mesh) -> int:
    return sum(len(s["nodes"]) for s in (mesh.boundary_segments or []) if s.get("kind") == "open")


def _flux_count(mesh) -> int:
    return sum(len(s["nodes"]) for s in (mesh.boundary_segments or [])
               if s.get("kind") == "flow" and s.get("ibtype") in FORT20_IBTYPES)


def _read_forcing(path, mesh, n_nodes, count_fn, what: str) -> BoundaryForcingSeries:
    nums = _first_numbers(path)
    if not nums:
        raise ValueError(f"{path}: empty file")
    dt, vals = nums[0], np.array(nums[1:], dtype=float)
    n = n_nodes
    if mesh is not None:
        expected = count_fn(mesh)
        if n is not None and n != expected:
            raise ValueError(f"n_nodes={n} but mesh has {expected} {what} boundary nodes")
        n = expected
        if n == 0:
            raise ValueError(f"mesh has 0 {what} boundary nodes but the file holds {vals.size} values")
    if n is None:
        return BoundaryForcingSeries(dt, vals)
    if n <= 0 or vals.size % n:
        raise ValueError(
            f"file holds {vals.size} values, which is not a whole number of steps of the "
            f"{n} {what} boundary nodes expected")
    return BoundaryForcingSeries(dt, vals.reshape(-1, n))


def _write_forcing(series: BoundaryForcingSeries, path) -> None:
    lines = [repr(float(series.time_increment))]
    lines += [repr(float(v)) for v in np.asarray(series.values, dtype=float).ravel()]
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def read_fort19(path, mesh=None, n_nodes: Optional[int] = None) -> BoundaryForcingSeries:
    """Read a fort.19 non-periodic elevation file.

    With ``mesh``, the per-step count is the number of open-boundary nodes and
    a file whose value count is not a whole number of steps raises ``ValueError``
    naming both counts. Without a mesh, pass ``n_nodes`` to split steps; with
    neither the flat series is returned (see :class:`BoundaryForcingSeries`).

    Example
    -------
    >>> s = read_fort19("fort.19", n_nodes=5)  # doctest: +SKIP
    >>> s.values.shape  # doctest: +SKIP
    (3, 5)
    """
    return _read_forcing(path, mesh, n_nodes, _open_count, "open")


def write_fort19(series: BoundaryForcingSeries, path) -> None:
    """Write a fort.19 file: time increment, then one value per line."""
    _write_forcing(series, path)


def read_fort20(path, mesh=None, n_nodes: Optional[int] = None) -> BoundaryForcingSeries:
    """Read a fort.20 non-periodic normal-flow file.

    Same rules as :func:`read_fort19`; matching nodes are the ``flow`` boundary
    nodes with IBTYPE in :data:`FORT20_IBTYPES`. Trailing ``!`` comments on value
    lines are ignored.

    Example
    -------
    >>> read_fort20("fort.20").time_increment  # doctest: +SKIP
    31536000.0
    """
    return _read_forcing(path, mesh, n_nodes, _flux_count, "normal-flow")


def write_fort20(series: BoundaryForcingSeries, path) -> None:
    """Write a fort.20 file: time increment, then one value per line."""
    _write_forcing(series, path)


def read_fort24(path, mesh=None) -> List[TidalConstituent]:
    """Read a fort.24 self-attraction/load-tide file.

    Node lines must be numbered 1..N in order. With ``mesh``, N must equal
    ``mesh.n_verts`` or ``ValueError`` names both counts.

    Example
    -------
    >>> c = read_fort24("fort.24")  # doctest: +SKIP
    >>> c[0].name, c[0].amp.shape  # doctest: +SKIP
    ('M2', (27330,))
    """
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.read().split("\n")
    out, pos = [], 0
    while pos < len(lines):
        if not lines[pos].strip():
            pos += 1
            continue
        if pos + 4 > len(lines):
            raise ValueError(f"{path}: truncated constituent header at line {pos + 1}")
        alpha = lines[pos].strip()
        try:
            freq = _float(lines[pos + 1].split()[0])
            nodal = _float(lines[pos + 2].split()[0])
        except (ValueError, IndexError):
            raise ValueError(f"{path}:{pos + 1}: malformed constituent header") from None
        name = lines[pos + 3].strip()
        pos += 4
        amp, pha = [], []
        while pos < len(lines):
            tok = lines[pos].split()
            if len(tok) < 3:
                break
            try:
                node = int(tok[0])
                a, p = _float(tok[1]), _float(tok[2])
            except ValueError:
                break
            if node != len(amp) + 1:
                break
            amp.append(a)
            pha.append(p)
            pos += 1
        if mesh is not None and len(amp) != mesh.n_verts:
            raise ValueError(
                f"constituent {name} has {len(amp)} nodes but mesh has {mesh.n_verts} nodes")
        out.append(TidalConstituent(name, alpha, freq, nodal, np.array(amp), np.array(pha)))
    if not out:
        raise ValueError(f"{path}: no constituents found")
    return out


def write_fort24(constituents: List[TidalConstituent], path) -> None:
    """Write fort.24: per constituent the alpha line, frequency, nodal factor, name, node lines."""
    lines = []
    for c in constituents:
        if len(c.amp) != len(c.phase):
            raise ValueError(f"constituent {c.name}: {len(c.amp)} amplitudes but {len(c.phase)} phases")
        lines += [c.alpha, repr(float(c.frequency)), repr(float(c.nodal_factor)), c.name]
        lines += [f"{i + 1} {float(a)!r} {float(p)!r}" for i, (a, p) in enumerate(zip(c.amp, c.phase))]
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
