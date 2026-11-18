"""W1 runner: every CHILmesh smoother and the external tools, Track 1 then Track 2."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .. import records, valence
from . import track1, track2, tools


def _load(path):
    from chilmesh import CHILmesh

    mesh = CHILmesh.read_from_fort14(str(path), compute_layers=False, compute_adjacencies=False)
    pts = np.asarray(mesh.points)
    depths = pts[:, 2] if pts.shape[1] > 2 else np.zeros(len(pts))
    return pts[:, :2], np.asarray(mesh.connectivity_list), depths


def run_mesh(name, out_dir, *, dt=1.0, courant_max=1.0, crs="cartesian", units="m", timeout_s=3600, cache="unknown"):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fetch = valence.fetch_mesh_file(name, out_dir / "meshes")
    meta = dict(name=name, version=fetch.get("sha256", "unknown")[:12], sha256=fetch.get("sha256", "unknown"),
                provenance="Valence registry", crs=crs, units=units, n_nodes=0, n_elems=0)
    recs = []
    if fetch["status"] != "ok":
        for track in (1, 2):
            for s in sorted(tools.CHILMESH_SMOOTHERS):
                recs.append(track1._failed(f"chilmesh:{s}", "PolyForm-NC", "subprocess", meta, track,
                                           fetch["reason"], [], status="not-run", cache=cache))
            for t, cfg in tools.EXTERNAL.items():
                recs.append(track1._failed(t, cfg["license"], cfg["run_mode"], meta, track, fetch["reason"], [],
                                           status="not-run", cache=cache))
    else:
        pts, conn, depths = _load(fetch["path"])
        for track, mod in ((1, track1), (2, track2)):
            for s in sorted(tools.CHILMESH_SMOOTHERS):
                extra = (lambda before, after: {"vertices_deleted": records.metric(int(len(before) - len(after)))}) if track == 2 else None
                rec = track1.run_chilmesh(s, pts, conn, depths, mesh=meta, crs=crs, dt=dt,
                                          courant_max=courant_max, track=track, frozen=track == 1,
                                          timeout_s=timeout_s, cache=cache, extra_metrics_fn=extra)
                if track == 2 and rec["status"] != "failed":
                    rec["reason"] = "; ".join(x for x in (rec["reason"], track2.CHILMESH_NOTE) if x)
                recs.append(rec)
            for t in tools.EXTERNAL:
                recs.append(mod.run_external(t, pts, conn, depths, mesh=meta, crs=crs, dt=dt,
                                             courant_max=courant_max, workdir=str(out_dir / f"work_{t}_{track}"),
                                             fort14_path=fetch["path"], timeout_s=timeout_s, cache=cache))
    for i, rec in enumerate(recs):
        records.write_record(rec, out_dir / f"w1_{name}_t{rec['track']}_{rec['tool'].replace(':', '-')}_{i}.json")
    return recs
