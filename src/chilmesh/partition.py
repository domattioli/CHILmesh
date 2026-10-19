"""Partition quality report for element labelings (spec 005 US3)."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["PartitionReport", "partition_report", "partition_mesh"]


@dataclass
class PartitionReport:
    """Per-part sizes, imbalance (max/mean), cut edges and connected pieces."""

    sizes: np.ndarray
    imbalance: float
    cut_edges: int
    pieces: np.ndarray


def _edge2elem(mesh) -> np.ndarray:
    # Meshes loaded with compute_layers=False carry no adjacencies (#293).
    if "Edge2Elem" not in mesh.adjacencies:
        mesh._build_adjacencies()
    return mesh.adjacencies["Edge2Elem"]


def partition_report(mesh, labels) -> PartitionReport:
    """Summarize a partition of mesh elements.

    ``labels`` has one integer per element, 0 <= label < k, with
    1 < k <= n_elems where k = max label + 1. Pieces count connected
    components of each part over edge (Edge2Elem) adjacency.

    Example:
        >>> from chilmesh import examples
        >>> mesh = examples.donut()
        >>> rep = partition_report(mesh, partition_mesh(mesh, 2))
        >>> int(rep.sizes.sum()) == mesh.n_elems, rep.pieces.tolist()
        (True, [1, 1])

    Tested by: tests/test_partition_report.py::test_sizes_imbalance_cut_pieces, ::test_disconnected_part, ::test_k_bounds.
    """
    labels = np.asarray(labels)
    n = mesh.n_elems
    if labels.shape != (n,):
        raise ValueError(f"labels length {labels.size} does not match n_elems {n}")
    if labels.size and labels.min() < 0:
        raise ValueError("labels must be non-negative")
    k = int(labels.max()) + 1 if n else 0
    if not (1 < k <= n):
        raise ValueError(f"k={k} out of range: need 1 < k <= n_elems ({n})")
    sizes = np.bincount(labels, minlength=k)
    imbalance = float(sizes.max() / sizes.mean())
    e2e = _edge2elem(mesh)
    inter = e2e[(e2e[:, 0] >= 0) & (e2e[:, 1] >= 0)]
    cut = int(np.count_nonzero(labels[inter[:, 0]] != labels[inter[:, 1]]))
    nbrs = [[] for _ in range(n)]
    for a, b in inter:
        if labels[a] == labels[b]:
            nbrs[a].append(int(b))
            nbrs[b].append(int(a))
    pieces = np.zeros(k, dtype=int)
    seen = np.zeros(n, dtype=bool)
    for s in range(n):
        if seen[s]:
            continue
        pieces[labels[s]] += 1
        stack = [s]
        seen[s] = True
        while stack:
            u = stack.pop()
            for v in nbrs[u]:
                if not seen[v]:
                    seen[v] = True
                    stack.append(v)
    return PartitionReport(sizes=sizes, imbalance=imbalance, cut_edges=cut, pieces=pieces)


def _dual_graph(mesh):
    import scipy.sparse as sp

    n = mesh.n_elems
    e2e = _edge2elem(mesh)
    inter = e2e[(e2e[:, 0] >= 0) & (e2e[:, 1] >= 0)]
    a, b = inter[:, 0], inter[:, 1]
    g = sp.coo_matrix((np.ones(2 * len(a)), (np.r_[a, b], np.r_[b, a])), shape=(n, n)).tocsr()
    g.data[:] = 1.0
    return g


def _fiedler(sub):
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla

    n = sub.shape[0]
    deg = np.asarray(sub.sum(axis=1)).ravel()
    lap = sp.diags(deg) - sub
    if n <= 600:
        w, v = np.linalg.eigh(lap.toarray())
        f = v[:, 1]
    else:
        v0 = np.linspace(1.0, 2.0, n)
        w, v = spla.eigsh(lap.tocsc(), k=2, sigma=-1e-3, which="LM", v0=v0)
        f = v[:, np.argsort(w)[1]]
    f = np.where(np.abs(f) < 1e-12, 0.0, f)
    if f[np.argmax(np.abs(f))] < 0:
        f = -f
    return f


def _repair(sub, side):
    """Move non-largest components of each side to the other side."""
    from scipy.sparse.csgraph import connected_components

    for _ in range(20):
        changed = False
        for val in (True, False):
            idx = np.where(side == val)[0]
            if len(idx) == 0:
                continue
            nc, comp = connected_components(sub[idx][:, idx], directed=False)
            if nc > 1:
                big = np.argmax(np.bincount(comp))
                side[idx[comp != big]] = not val
                changed = True
        if not changed:
            break
    return side


def _bisect(g, ids, k, labels, base):
    if k == 1:
        labels[ids] = base
        return
    k1 = k // 2
    sub = g[ids][:, ids]
    f = _fiedler(sub)
    order = np.lexsort((np.arange(len(ids)), f))
    n1 = int(round(len(ids) * k1 / k))
    n1 = min(max(n1, k1), len(ids) - (k - k1))
    side = np.zeros(len(ids), dtype=bool)
    side[order[:n1]] = True
    side = _repair(sub, side)
    _bisect(g, ids[side], k1, labels, base)
    _bisect(g, ids[~side], k - k1, labels, base + k1)


def partition_mesh(mesh, k: int, method: str = "spectral") -> np.ndarray:
    """Partition mesh elements into ``k`` parts; returns int labels in [0, k).

    ``method="spectral"``: recursive spectral bisection of the element dual
    graph (edge-sharing adjacency, scipy only), splitting the Fiedler vector at
    the size-proportional position, then reassigning stray components so each
    part is connected. Deterministic. ``method="metis"`` needs the optional
    ``partition`` extra (``pip install chilmesh[partition]``).

    Example:
        >>> from chilmesh import examples
        >>> mesh = examples.donut()
        >>> labels = partition_mesh(mesh, 4)
        >>> partition_report(mesh, labels).imbalance <= 1.25
        True

    Tested by: tests/test_partition_mesh.py::test_spectral_threshold, ::test_k_range, ::test_metis_missing.
    """
    n = mesh.n_elems
    if method not in ("spectral", "metis"):
        raise ValueError(f"unknown method {method!r}; use 'spectral' or 'metis'")
    if not (1 <= k <= n):
        raise ValueError(f"k={k} out of range: need 1 <= k <= n_elems ({n})")
    g = _dual_graph(mesh)
    if method == "metis":
        try:
            import pymetis
        except ImportError as exc:
            raise ImportError(
                "method='metis' needs pymetis; install the 'partition' extra: "
                "pip install chilmesh[partition]"
            ) from exc
        adj = [g.indices[g.indptr[i]:g.indptr[i + 1]].tolist() for i in range(n)]
        _, parts = pymetis.part_graph(k, adjacency=adj)
        return np.asarray(parts, dtype=int)
    labels = np.zeros(n, dtype=int)
    _bisect(g, np.arange(n), k, labels, 0)
    return labels
