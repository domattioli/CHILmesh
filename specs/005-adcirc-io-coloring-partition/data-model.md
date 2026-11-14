<!-- provenance: author=domattioli model=claude-opus-5-5 effort=medium date=2026-10-12 skill=speckit-pipeline repo=CHILmesh session=4e33ac09-e474-4ebc-a2dd-027a9e0a0442 -->
# Data Model: #289 B1-B3

| Entity | Fields | Rules |
|---|---|---|
| Nodestring | `node_ids: int[]` (ordered) | >=2 nodes; IDs exist in mesh |
| ElementColoring | `colors: int[n_elems]`, `n_colors: int` | no two node-sharing elements equal; `n_colors <= max_degree+1` |
| PartitionLabels | `labels: int[n_elems]`, `k: int` | 0 <= label < k; 1 < k <= n_elems |
| PartitionReport | `sizes: int[k]`, `imbalance: float` (max/mean), `cut_edges: int`, `pieces: int[k]` | `pieces[i] >= 1` for non-empty part |
| OutputSeries | `times: float[nt]`, `iterations: int[nt]`, `values: float[nt, n]` or `[nt, n, 2]`, `kind: scalar|vector`, `location: nodes|stations`, `dropped_records: int` | NaN where dry/fill; n matches mesh when mesh given; maxele/maxvel second set (time of max) is record index 1 |
| StormTrack | records of `time, lat, lon, vmax, pmin, isotach, radii[4], rmw, extra: str[]` (untyped columns, record order) | time ascending; missing as None; `extra` written back unchanged |
| OWISnapshot | `time, ilat, ilon, dx, dy, swlat, swlon, values` (pressure) or `u, v` | shape `(ilat, ilon)` |
| BoundaryForcingSeries | `time_increment, values: float[nt, nbn]` | nbn equals count of matching boundary nodes when mesh given |
| TidalLoadField (fort.24) | per constituent `alpha, frequency, nodal_factor, name, amp[n], phase[n]` (no equilibrium argument; see research.md R8) | n equals mesh nodes |
| ManifestEntry | `name, format, stage, url, lfs: bool, sha256, size` | sha256 64 hex |
