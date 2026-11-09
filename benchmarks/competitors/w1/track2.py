"""Track 2 (secondary): full clean. Vertices may be deleted by competitors; CHILmesh deletes none."""
from __future__ import annotations

from .. import records
from . import track1

TRACK = 2
CHILMESH_NOTE = "CHILmesh deletes no vertices (coordinates move, topology is constant)"


def run_chilmesh(smoother, points, conn, depths, *, mesh, crs, dt, courant_max=1.0, cache="unknown",
                 repeat_index=0):
    rec = track1.run_chilmesh(
        smoother, points, conn, depths, mesh=mesh, crs=crs, dt=dt, courant_max=courant_max, track=TRACK,
        frozen=False, cache=cache, repeat_index=repeat_index,
        extra_metrics_fn=lambda before, after: {"vertices_deleted": records.metric(int(len(before) - len(after)))})
    if rec["status"] != "failed":
        rec["reason"] = "; ".join(x for x in (rec["reason"], CHILMESH_NOTE) if x)
    return rec


def run_external(tool, points, conn, depths, **kw):
    return track1.run_external(tool, points, conn, depths, track=TRACK, frozen=False, **kw)
