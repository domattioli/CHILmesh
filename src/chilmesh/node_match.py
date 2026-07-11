"""Cross-mesh node matching + nodal-field deltas (#239).

A numpy-native, renumbering-tolerant node matcher between two meshes, plus a
generic per-node field-delta helper. Ports the grid-hash nearest-neighbour
matcher from domattioli/Valence (``packages/valence_domains/diff.py``
``match_nodes``) into CHILmesh, using scipy's cKDTree for the spatial query
(CHILmesh already depends on scipy) instead of a hand-rolled grid hash.

The matcher survives node renumbering because it matches on coordinates, not
ids; ``nodal_field_delta`` then compares any per-node scalar (bathymetry,
representative edge length / resolution, ...) between the two meshes over the
matched pairs.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

DEFAULT_TOLERANCE_FRACTION = 0.5


@dataclass
class NodeMatch:
    """Result of matching mesh A nodes onto mesh B nodes.

    Indices are 0-based positional indices into the respective point arrays.
    """
    matched: dict[int, int]   # a_idx -> b_idx
    added: list[int]          # b indices with no a counterpart
    removed: list[int]        # a indices with no b counterpart
    method: str               # "identity" | "spatial"


def _as_xy(points) -> np.ndarray:
    arr = np.asarray(points, dtype=float)
    if arr.ndim != 2 or arr.shape[1] < 2:
        raise ValueError("points must be an (N, >=2) array of coordinates")
    return arr[:, :2]


def _median_nn_distance(xy: np.ndarray) -> float:
    """Median nearest-neighbour distance (mesh spacing proxy)."""
    if len(xy) < 2:
        return 0.0
    tree = cKDTree(xy)
    # k=2: column 0 is the point itself (dist 0), column 1 is the neighbour.
    dists, _ = tree.query(xy, k=2)
    return float(np.median(dists[:, 1]))


def derive_tolerance(
    points_a, points_b, *, fraction: float = DEFAULT_TOLERANCE_FRACTION
) -> float:
    """Data-derived match radius = ``fraction`` * finer mesh's median NN spacing.

    A node that moved less than this is treated as the same node displaced, not
    an add + remove. Mirrors Valence's ``derive_tolerance`` but derives spacing
    from nearest-neighbour distances (no connectivity required).
    """
    xa, xb = _as_xy(points_a), _as_xy(points_b)
    spacings = [s for s in (_median_nn_distance(xa), _median_nn_distance(xb)) if s > 0]
    spacing = min(spacings) if spacings else 1.0
    tol = fraction * spacing
    return tol if tol > 0 else 1e-9


def _identity_holds(xa: np.ndarray, xb: np.ndarray, tolerance: float) -> bool:
    if xa.shape != xb.shape or len(xa) == 0:
        return False
    return bool(np.all(np.hypot(xb[:, 0] - xa[:, 0], xb[:, 1] - xa[:, 1]) <= tolerance))


def match_nodes(points_a, points_b, *, tolerance: float) -> NodeMatch:
    """Match A nodes onto B nodes within ``tolerance`` (renumbering-tolerant).

    Identity fast-path: when the arrays have equal length and every aligned pair
    agrees within ``tolerance``, they are the same nodes in the same order.
    Otherwise a greedy spatial nearest-neighbour match (ascending A index, each B
    node used at most once) via cKDTree over the (up to) 8 nearest candidates.
    """
    xa, xb = _as_xy(points_a), _as_xy(points_b)

    if _identity_holds(xa, xb, tolerance):
        return NodeMatch({i: i for i in range(len(xa))}, [], [], "identity")

    if len(xb) == 0:
        return NodeMatch({}, [], list(range(len(xa))), "spatial")

    tree = cKDTree(xb)
    k = min(8, len(xb))
    dists, idxs = tree.query(xa, k=k)
    if k == 1:  # cKDTree drops the trailing dim when k == 1
        dists = dists[:, None]
        idxs = idxs[:, None]

    matched: dict[int, int] = {}
    removed: list[int] = []
    used: set[int] = set()
    for aid in range(len(xa)):
        chosen = None
        for d, bid in zip(dists[aid], idxs[aid]):
            if d <= tolerance and int(bid) not in used:
                chosen = int(bid)
                break
        if chosen is None:
            removed.append(aid)
        else:
            matched[aid] = chosen
            used.add(chosen)
    added = sorted(set(range(len(xb))) - used)
    return NodeMatch(matched, added, removed, "spatial")


def _stats(values: np.ndarray) -> dict:
    if len(values) == 0:
        return {"count": 0, "max": 0.0, "median": 0.0, "mean": 0.0}
    return {
        "count": int(len(values)),
        "max": float(np.max(values)),
        "median": float(np.median(values)),
        "mean": float(np.mean(values)),
    }


def nodal_field_delta(match: NodeMatch, field_a, field_b) -> dict:
    """Per-matched-node delta ``field_b[b] - field_a[a]`` with signed rollups.

    ``field_a`` / ``field_b`` are per-node scalar arrays aligned with the point
    arrays passed to :func:`match_nodes`. Use for bathymetry (positive-down: +Δ
    deepened), representative edge length (resolution: +Δ coarser), or any other
    per-node scalar. Returns matched-count, ``delta`` (a_idx -> signed Δ), and
    ``increased`` / ``decreased`` magnitude stats.
    """
    fa = np.asarray(field_a, dtype=float)
    fb = np.asarray(field_b, dtype=float)
    delta: dict[int, float] = {}
    inc: list[float] = []
    dec: list[float] = []
    for aid, bid in match.matched.items():
        d = float(fb[bid] - fa[aid])
        if d > 0.0:
            delta[aid] = d
            inc.append(d)
        elif d < 0.0:
            delta[aid] = d
            dec.append(-d)
    return {
        "matched_nodes": len(match.matched),
        "changed_nodes": len(delta),
        "delta": delta,
        "increased": _stats(np.asarray(inc)),
        "decreased": _stats(np.asarray(dec)),
    }


__all__ = [
    "NodeMatch",
    "match_nodes",
    "derive_tolerance",
    "nodal_field_delta",
    "DEFAULT_TOLERANCE_FRACTION",
]
