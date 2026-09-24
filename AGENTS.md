# CHILmesh

CHILmesh is a Python library for generating, manipulating, smoothing, and analyzing 2D triangular, quadrilateral, and mixed-element meshes. It supports ADCIRC `fort.14`, SMS `.2dm`, related ADCIRC formats, concentric mesh layerization, mesh-quality analysis, and spatial queries. Python is the reference implementation. Optional C++ and frozen Rust backends provide acceleration and must remain output-equivalent to Python.

The original MATLAB implementation remains in `src/@CHILmesh/CHILmesh.m` but is not actively developed.

## Project hard rules

- Keep public method signatures stable. Breaking changes require deprecation warnings.
- Preserve lossless `fort.14` round trips, including boundary metadata.
- Preserve layer-peel behavior and layer invariants.
- Keep triangular, quadrilateral, and mixed-element support.
- Keep Python and compiled backend outputs equivalent. Run backend-equivalence tests after backend changes.
- Require type hints on public APIs.
- Validate file input and user input at boundaries.
- Document every adjacency invariant.

## Repository layout

```text
CHILmesh/
├── src/chilmesh/              # Python package and reference algorithms
│   ├── CHILmesh.py            # Mesh class, adjacency build, and layer peel
│   ├── fort13_io.py           # ADCIRC fort.13 I/O
│   ├── fort14_io.py           # ADCIRC fort.14 I/O
│   ├── fort15_io.py           # ADCIRC fort.15 I/O
│   ├── mesh_topology.py       # Topology helpers
│   ├── backends/              # Python wrappers for compiled backends
│   └── data/                  # Built-in mesh fixtures
├── src/chilmesh_cpp/          # Optional C++ half-edge backend
├── src/chilmesh_core/         # Frozen Rust backend
├── src/@CHILmesh/             # Original MATLAB implementation
├── tests/                     # Pytest suite
├── examples/                  # Runnable examples
├── scripts/                   # Build, benchmark, release, and hook scripts
├── docs/                      # Architecture, API, format, and benchmark docs
├── .planning/                 # Project plans, audits, and decisions
└── pyproject.toml             # Package metadata and dependencies
```

Key implementation points:

- `_build_adjacencies()` in `src/chilmesh/CHILmesh.py` constructs mesh topology.
- `_peel()` performs concentric layerization. It is not medial-axis extraction.
- `layers` contains `OE`, `IE`, `OV`, `IV`, and `bEdgeIDs` entries per layer.
- `compute_layers=True` requires adjacencies. With `compute_layers=False`, callers may independently request adjacencies and spatial indices.

## Install and run

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"

# CLI
chilmesh info mesh.fort.14
python -m chilmesh --help

# Optional C++ backend from source
pip install ./src/chilmesh_cpp
# or
bash scripts/build_cpp.sh
```

## Test commands

```bash
# Full suite
pytest -v

# Fast suite
pytest -m "not slow" -v

# Specific invariant
pytest tests/test_invariants.py::test_layers_disjoint_cover -v

# Coverage
pytest --cov=src/chilmesh --cov-report=html tests/

# Backend parity
pytest tests/test_backend_equivalence.py -v
```

See `tests/TESTING.md` for markers, backend setup, and external MATLAB parity tests.

## Code and data conventions

- Support Python 3.10 and newer.
- Use NumPy arrays for dense mesh data and dictionaries of sets for sparse vertex adjacencies.
- Keep internal refactors behind the existing public API.
- Comments explain why a choice exists. Do not narrate obvious code.
- Use NumPy-style docstrings for public APIs.
- Add a reproducer for every bug fix.
- Add feature-specific tests for new behavior.
- Run the existing suite unchanged for refactors.
- Include benchmarks and parity tests for algorithmic or backend optimizations.

### Adjacency invariants

| Structure | Representation | Invariant |
|---|---|---|
| `Elem2Vert` | `ndarray[n_elems, 3|4]` | Element vertex IDs are valid. Mixed triangles use the established padding convention. |
| `Edge2Vert` | `ndarray[n_edges, 2]` | Endpoints are stored in canonical `(min, max)` order. |
| `Elem2Edge` | `ndarray[n_elems, 3|4]` | Edge IDs follow element-major, slot-minor first-encounter ordering. |
| `Edge2Elem` | `ndarray[n_edges, 2]` | Boundary edges use `-1` for the missing adjacent element. |
| `Vert2Edge` | `dict[int, set[int]]` | Every vertex has an entry containing all incident edge IDs. |
| `Vert2Elem` | `dict[int, set[int]]` | Every vertex has an entry containing all incident element IDs. |
| `EdgeMap` | hash map | Canonical endpoint pairs map to edge IDs in constant expected time. |

Topology mutations must rebuild or explicitly invalidate affected adjacencies and layers. Do not return stale topology.

### Built-in test fixtures

`tests/conftest.py` exposes five fixtures: `annulus`, `donut`, `block_o`, `structured`, and `quad_2x2`. Tests that mutate a cached mesh must copy it first. Parametrize across all applicable element types. Use `TRI_FIXTURE_NAMES` for triangle-only behavior.

## Branch workflow

- The default working branch is `development` because `origin/development` exists.
- Releases flow through a pull request from `development` to `main`.
- Never push directly to `main`.
- Never force-push.
- CI must pass before merge.

## Related repositories

| Repository | Relationship |
|---|---|
| MADMESHR | Downstream research project that uses CHILmesh. |
| ADMESH | Mesh adaptation framework that uses CHILmesh. |
| ADMESH-Domains | Domain handling used with ADMESH and CHILmesh. |
| Valence | Registry used by reference benchmark workloads. |

## Project references

- `.planning/constitution.md`: repository governance and project principles.
- `.planning/project_plan.md`: roadmap and milestone status.
- `.planning/MODERNIZATION_LESSONS_LEARNED.md`: optimization decisions and tradeoffs.
- `.planning/TEST-AUDIT.md`: test coverage audit.
- `.planning/HOOKS-AUDIT.md`: hook audit and historical findings.
- `docs/ARCHITECTURE.md`: graph engine, adjacency tables, and layerization boundary.
- `docs/ADJACENCY_STRUCTURES.md`: adjacency details.
- `docs/RUST_EVALUATION.md`: rationale for freezing the Rust backend.

## Governance

This repository is a downstream consumer of `domattioli/DomI`. Universal git, coding-dispatch, secrets, session-lifecycle, and communication rules live in DomI `.claude/policies/`.
The `.domi-pin` drift check runs at session start through `scripts/instructions_on_start.sh`.
Spec-kit artifacts for CHILmesh live in DomI `specs/consumers/CHILmesh/`, never in a local `.specify/` directory.
