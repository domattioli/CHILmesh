<!-- provenance: author=domattioli model=claude-opus-5-5 effort=medium date=2026-10-12 skill=speckit-pipeline repo=CHILmesh session=4e33ac09-e474-4ebc-a2dd-027a9e0a0442 -->
# Public API contract (additive)

```python
CHILmesh.write_to_2dm(filename, *, nodestrings=True) -> None   # _write_2dm stays as alias
CHILmesh.read_from_2dm(filename, ...)                           # now fills boundary segments from NS cards
chilmesh.color_elements(mesh) -> ElementColoring
chilmesh.check_coloring(mesh, colors) -> bool
chilmesh.partition_report(mesh, labels) -> PartitionReport
chilmesh.partition_mesh(mesh, k, method="spectral"|"metis")    # ONLY if exploration fails
chilmesh.read_adcirc_output(path, mesh=None, step=None) -> OutputSeries        # fort.61/62/63/64, maxele, maxvel
chilmesh.read_adcirc_output_nc(path, mesh=None, step=None, variable=None) -> OutputSeries  # needs extra 'netcdf'; fort.61/62/63/64 .nc, maxele, maxvel
chilmesh.chilplotting.plot_field(mesh, values, *, ax=None, cmap="viridis", colorbar=True)
chilmesh.read_fort22_besttrack(path, nws=None) -> StormTrack ; write_fort22_besttrack(track, path, nws=8|20)
chilmesh.read_owi(path) -> list[OWISnapshot] ; write_owi(snapshots, path)
chilmesh.read_fort19 / write_fort19, read_fort20 / write_fort20, read_fort24 / write_fort24   # read_* take mesh=None; count check when given. read_fort19/read_fort20 also take n_nodes=None (no mesh needed)
```

Errors: `ValueError` on node-count mismatch (message names both counts), unsupported NWS, k out of range. `ImportError` naming the extra when netCDF4 or pymetis missing. `UserWarning` on truncated output naming the dropped count.

Script: `python scripts/check_adcirc_testsuite.py [--stage B1|B2|B3|all] [--cache DIR]` → exit 0 only when every pinned file for the selected stages downloads, hashes and parses.
