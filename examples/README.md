# CHILmesh examples

Runnable scripts demonstrating common CHILmesh tasks. Each script is
self-contained — clone the repo (or `pip install chilmesh`) and run.

| Script | What it shows |
|---|---|
| [`01_quickstart.py`](01_quickstart.py) | Load a built-in mesh, print stats, save a plot |
| [`02_fort14_roundtrip.py`](02_fort14_roundtrip.py) | Write a mesh to `.fort.14` and read it back |
| [`03_smoothing.py`](03_smoothing.py) | Perturb interior vertices and relax them with the angle-based smoother |
| [`04_spatial_queries.py`](04_spatial_queries.py) | Phase 5 spatial-indexing API: `find_element`, `find_elements_in_radius`, `nearest_vertices` |

## Running

```bash
# from repo root
python examples/01_quickstart.py
```

Examples use only the public API surface (`chilmesh.CHILmesh`,
`chilmesh.examples`, `chilmesh.write_fort14`). They run against the
bundled fixtures in `chilmesh.data`, so no external mesh files are needed.
