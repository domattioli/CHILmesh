"""Boundary-node displacement check (FR-006).

Boundary nodes are fixed by the *input* mesh. Their coordinates must be
bit-identical after the run (``np.array_equal``); any move, even one ulp, is
``invalid``. With ``frozen_topology=True`` (Track 1) node and element counts must
also be unchanged.
"""
from __future__ import annotations

import numpy as np

from .quality import boundary_nodes


def check_boundary(points_before, conn_before, points_after, conn_after, *, frozen_topology: bool = True) -> dict:
    before = np.asarray(points_before, dtype=float)[:, :2]
    after = np.asarray(points_after, dtype=float)[:, :2]
    out = {"status": "ok", "reason": "", "n_boundary": 0, "n_boundary_moved": 0,
           "max_boundary_displacement": 0.0}
    if frozen_topology:
        if len(before) != len(after):
            return {**out, "status": "invalid", "reason": f"node count changed ({len(before)} -> {len(after)}) in Track 1"}
        if len(np.asarray(conn_before)) != len(np.asarray(conn_after)):
            return {**out, "status": "invalid",
                    "reason": f"element count changed ({len(np.asarray(conn_before))} -> {len(np.asarray(conn_after))}) in Track 1"}
    b = boundary_nodes(conn_before)
    b = b[b < min(len(before), len(after))]
    out["n_boundary"] = int(len(b))
    if np.array_equal(before[b], after[b]):
        return out
    moved = np.any(before[b] != after[b], axis=1)
    disp = np.linalg.norm(before[b] - after[b], axis=1)
    out.update(status="invalid", n_boundary_moved=int(moved.sum()),
               max_boundary_displacement=float(disp.max()),
               reason=f"{int(moved.sum())} boundary node(s) moved (max {float(disp.max()):.3g})")
    return out
