# Changelog

All notable changes to this project will be documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and the project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### 🔧 Fixed

- **Windows portability of fort.14 / 2dm I/O** (#121, partial). The three
  `open()` sites in `chilmesh.CHILmesh` now pass explicit
  `encoding='utf-8'` and `newline=` arguments. Reads use
  `newline=''` (universal-newline translation) so CRLF-terminated files
  authored on Windows roundtrip identically to LF-terminated files
  authored on POSIX; writes force `newline='\n'` so the produced bytes
  are byte-identical across all platforms. Adds
  `tests/test_io_portability.py` to lock the invariant. The `windows-latest`
  CI runner and OS-portable `build-and-smoke` venv path remain open
  under #121.

## [0.4.1] — 2026-05-18 (Consumer-Readiness + Zenodo Patch)

Patch release primarily for Zenodo first-archive and consumer-facing
documentation polish. No public-API breaking changes.

### ✨ Added

- `plot_quality_histogram(bins, cmap, auto_norm, ax)` — 1-D companion to
  `plot_quality`. Bars are coloured by colormap midpoint so the
  distribution reads as a horizontal slice of the 2-D quality plot.
  `auto_norm=True` compresses the colormap to the observed quality range
  for narrow distributions; default `False` keeps the [0, 1] norm to
  match `plot_quality` colours pixel-for-pixel.
- `scripts/generate_wnat_showcase.py` — stacks `plot_quality` over
  `plot_quality_histogram` to regenerate `output/wnat_hagen_showcase.png`.
  Falls back to the largest bundled fixture if WNAT_Hagen.14 is not on
  disk; honours `--mesh` / `WNAT_HAGEN_PATH` env var.
- `CITATION.cff` at repo root — GitHub renders the "Cite this repository"
  button. Includes author / affiliation / version metadata; identifiers
  block ready to fill once Zenodo mints DOIs.
- `conda/meta.yaml` + `conda/README.md` — conda-forge recipe scaffold +
  submission workflow.
- `.planning/ZENODO-SETUP.md` — five-minute one-time auth runbook plus
  per-release maintenance steps.
- README Gallery: figures wrapped in centred `<p>` blocks with
  `<sub><em>` caption styling, figure numbers, and tighter
  functionality-oriented captions.
- README Performance: compact v0.2.0 → v0.4.0 comparison table back at
  the top; full breakdown stays in `docs/BENCHMARK.md`.
- README badges: Zenodo DOI badge added (repo id 693749657).
- README "Citation": 1-2 line context paragraph above the bibtex block
  explaining the QuADMESH+ → CHILmesh lineage.
- README Installation: snippets for pip, [uv](https://docs.astral.sh/uv/),
  conda-forge (pending acceptance), and from-source.

### 🛠 Build & CI

- `pyproject.toml`: `pytest-xdist` added to `[dev]` extras.
- `.github/workflows/python-package.yml`: PR-only fast matrix
  (`ubuntu-latest × py3.11`, `pytest -n auto -m "not slow"`); push-to-main
  runs full matrix (Ubuntu+macOS × py3.10/3.11/3.12). Pip cache via
  `actions/setup-python cache: pip`. Wheel build + twine + smoke install
  split into `build-and-smoke` job gated on
  `github.event_name != 'pull_request'`. Expected PR cycle 10 min → 2 min.
- `scripts/generate_3row_admesh.py`: layer colorbar uses
  `BoundaryNorm` + integer ticks (was continuous 0–1).

### 🧹 Internal

- Repo-root cleanup: 5 audit/plan docs + `docs/introspections/`
  relocated to `.planning/`.
- README trimmed 288 → 224 lines; install command above the fold; deep
  content lives under `docs/`.
- `.specify/memory/constitution.md` consolidated as canonical
  governance; `.planning/constitution.md` retained as redirect stub.

### 📚 Documentation

- "Mesh generation" wording removed from README tagline and
  `docs/API.md` overview — replaced with "processing, smoothing, and
  analysis" to accurately scope the library. Advancing-front-related
  API references retain the "generation" terminology where it applies.
- Smoothing section rewritten to enumerate the three algorithms
  (Balendran direct FEM, Zhou-Shimada angle-based, ADMESH Spring-Based
  Truss Smoother)
  in a comparison table with selection guidance.
- `CITATION.cff` keyword `mesh-generation` → `mesh-processing` /
  `mesh-smoothing`.

### 🏷 Issue label hygiene

- Issues #101, #105 normalised onto the `priority:*` convention
  (was `low-priority` / `high-priority` legacy form).
- Issues #94, #93, #110, #111 received explicit `priority:*` labels.
- Issue #92 (Phase 5 spatial indexing) closed as completed — shipped
  in v0.4.0.

---

## [0.4.0] — 2026-05-18 (Spatial Indexing + Layer-Paths + CLI + Audit Polish)

Builds on the v0.3.0 vectorisation foundation with first-class spatial
queries, a shell-friendly CLI, an outer-vertex traversal port, and
broad consumer-readiness polish.

### ✨ Added

#### Command-line interface (issue #117)
- `chilmesh` CLI with four subcommands wrapping the existing public API
  - `chilmesh info <mesh>` — vert/elem/edge/layer counts + skew-quality stats
  - `chilmesh convert <in> <out>` — format conversion (.14/.fort.14/.2dm in, .14/.fort.14 out)
  - `chilmesh smooth <in> -o <out> --method {angle-based,fem,laplacian} --iter N` — in-place relaxation; reports median/min quality delta
  - `chilmesh plot <in> -o <fig> [--layers|--quality]` — static PNG/PDF/SVG by suffix
- New module `src/chilmesh/cli.py` + `__main__.py` so `python -m chilmesh ...` also works
- `[project.scripts] chilmesh = "chilmesh.cli:main"` entry point in `pyproject.toml`
- `tests/test_cli.py` exercises each subcommand via `subprocess.run([sys.executable, "-m", "chilmesh", ...])` on the bundled `annulus` fixture
- No new runtime dependencies — stdlib `argparse` over the existing public API

#### Layer-paths feature (PR #118)
- `chilmesh.layer_paths` module — port of MATLAB `PathsOnOV` vertex-ordering heuristic for outer-vertex traversal
- `tests/test_layer_paths.py` — regression suite (207 LOC)
- Scoped subgraph build: `O(L*m)` → `O(m)` per layer
- Demo visualization: `output/annulus_layer1_paths.png`

#### Phase 5 — Spatial query API (PR #115)
- `find_element(point)` — locate the element containing a 2D point via centroid kd-tree narrowing + barycentric check; returns `-1` outside the mesh
- `find_elements_in_radius(point, radius)` — ball-tree query returning all elements within `radius` of a point
- `nearest_vertices(point, k)` — k-nearest-neighbour vertex lookup
- `insert_vertex(point)` — vertex insertion primitive (initial mesh-mutation surface; see #94 for the full Phase 5 mutation API)
- Lazy-built kd-tree on element centroids; bounded `< 0.5s` build on WNAT_Hagen, `< 50μs` per `find_element` call

#### Examples directory
- `examples/01_quickstart.py`, `02_fort14_roundtrip.py`, `03_smoothing.py`, `04_spatial_queries.py` — runnable consumer scripts against the bundled fixtures; no external mesh files required

### 🧪 Testing

#### Test-suite holistic audit (#110, see `.planning/TEST-AUDIT.md`)
- **F1** — `tests/test_angle_based_smoother.py` (21 tests) covering the Zhou-Shimada angle-based smoother; `CHILmesh.py` line coverage 73% → 89%
- **F2** — `tests/test_plot_utils.py` (104 tests) promoted from smoke tests to real assertions (artist counts, axis limits, titles, colormaps)
- **F4** — `TRI_FIXTURE_NAMES` centralised in `conftest.py`
- **F5** — `tests/test_copy.py` (65 tests) locks the deep-copy contract the fixture cache relies on
- **F6** — assertion messages added to numerical comparisons in 5 test files
- **F7** — `tests/test_2dm_reader.py` migrated to pytest's `tmp_path` fixture
- **F9** — `slow` marker registered in `pyproject.toml`
- **F11** — README test badge points at the correct workflow
- **F14** — defensive-branch `# pragma: no cover` annotations on truly-unreachable guards in `_skeletonize` and `pinch_points`; behavioural lock-in tests added for `pinch_points` semantics
- Total: 582 tests passing / 13 skipped (was 439 / 9 in the original audit); `CHILmesh.py` line coverage 89% → 90%, total 88% → 89%

### 📚 Documentation

- README trimmed 288 → 224 lines: install command now visible on first laptop screen; deep content (BENCHMARK, API) lives under `docs/`. All three showcase images preserved in a dedicated Gallery section
- README ToC added (#106); `Contributing` section resolves previously-stale ToC link
- `examples/README.md` indexes the runnable scripts
- `docs/API.md`: Visualization section added; `plot_boundary` / `plot_interior_edges` documented
- `TESTING.md` added: pytest invocation patterns, fixture table, debugging tips
- Constitution consolidated into `.specify/memory/constitution.md` (#107); `.planning/constitution.md` and `.specify/speckit-constitution.md` retained as redirect stubs

### 🏛 Governance

- CHILmesh Constitution v1.0 → v1.1: Decision Authority table, Breaking Changes policy, Release Process, Feature-Specific Principles (I/O, Mesh Smoothing with DOMsmooth fallback, Skeletonization), External Upstream DomI section

### 🧹 Internal

- `matplotlib >= 3.10` compat: `axis_chilmesh` aspect-ratio assertions accept both `'equal'` and numeric `1.0`
- Repo-root cleanup: 5 audit/plan docs and `docs/introspections/` relocated into `.planning/` for a cleaner GitHub landing-page impression; no packaging change (sdist already pruned `.planning/`)
- DomI sync contract: `.domi-pin` and `scripts/instructions_on_start.sh` enforce upstream-skill drift detection at session start (#81, #112, #113)
- HOOKS-AUDIT.md inventoried hook surface and identified upstream-relevant gaps (#112)

---

## [0.3.0] — 2026-05-10 (Vectorization & Benchmark Validation)

Tag `v0.3.0` (commit `cf8f399`, PR #91) shipped the core-operation vectorisation pass and the reproducible benchmark harness.

### ✨ Added

#### New public methods
- `boundary_node_indices()` — returns unique boundary vertex indices (PR #86, closes README API gap)
- `plot_boundary()` — convenience plot for boundary edges only
- `plot_interior_edges()` — convenience plot for interior edges only

#### Reproducible benchmarking
- `scripts/benchmark_wnat_hagen.py` — measures init/quality/query performance on WNAT_Hagen reference mesh; markdown report + optional `--json` for CI archival (issue #55)

### ⚡ Performance

#### Core operation vectorization (issue #75)
- `signed_area`: O(n²) loop → fully vectorised shoelace formula on `(n,3,2)` / `(n,4,2)` arrays
- `interior_angles`: outer loop replaced with numpy gather `self.points[verts, :2]`
- `elem_quality`: O(n²) `[elem_id in tri_elems for elem_id in elem_ids]` → boolean mask from connectivity shape
- `_plot_polys`: per-element Python loop → numpy fancy indexing fed to `PolyCollection`
- `plot_point` centroid path: loop → vectorised `points[conn,:2].mean(axis=1)`

#### Measured impact (WNAT_Hagen, 52,774 verts · 98,365 elems)
- Total workflow: 14.3 s → **3.33 s** (4.3× faster than v0.2.0; ~4,000× vs the v0.1.1 baseline)
- Quality analysis: 6.6 s → **0.07 s** (94× faster than v0.2.0)
- `Vert2Edge` lookup: 0.7 μs → **0.17 μs** per call

These numbers carry forward unchanged into v0.4.0 — Phase 5 spatial queries and layer-paths are additive features, not workflow-level speedups.

### 🛠 Packaging

- `MANIFEST.in` added: PyPI sdist excludes `.claude/`, `.planning/`, `.specify/`, `specs/`, `tests/`, `scripts/`, `.github/`, `.domi-pin`. Wheel layout unchanged.

### 📚 Documentation

- `docs/BENCHMARK.md`: appended May 2026 re-run section preserving original April 2026 baselines

## [0.2.1] — TBD (Bug Fix Release)

### 🐛 Fixed

#### Skeletonization Layer Separation Invariant (#74)
- Fixed layer separation invariant violations in `_skeletonize()` method
- Replaced iterative-peel algorithm with MATLAB-equivalent BFS distance-based skeletonization
- Layer counts now match original MATLAB `meshLayers` function output
- Layer assignment is now correct: no elements from layer k share vertices with elements from layer m where |k-m| ≥ 2
- Comprehensive regression tests added: `test_skeletonization_invariant.py`, `test_skeletonization_matlab_parity.py`
- **Impact**: Layer counts may change for existing meshes (this is correct behavior; previous counts were wrong)
- **Backward Compatibility**: API remains unchanged; layer data structure format unchanged

### ✨ Added

#### Enhanced Test Coverage
- `test_skeletonization_invariant.py`: Layer separation invariant validation across all fixtures
- `test_skeletonization_matlab_parity.py`: MATLAB reference layer count validation
- `test_skeletonization_matlab_parity_external.py`: External ADMESH-Domains catalog parity validation
- `test_smoothing.py`: Comprehensive FEM smoother tests for triangles, quads, and mixed-element meshes

## [0.2.0] — 2026-04-27 (Modernization Release)

This is a major performance and architecture modernization release. **All existing APIs remain backward compatible** (no code changes needed) while internal data structures are redesigned for 150x+ performance improvement on large meshes.

### 🎯 Breaking Changes (Internal Only)

**User-facing API: 100% backward compatible**
- All v0.1.1 code works without modification
- Public method signatures unchanged
- Existing test suite passes as-is

**Internal changes (if you accessed private structures):**
- `adjacencies['Vert2Edge']`: Changed from `List[List[int]]` to `Dict[int, Set[int]]`
- `adjacencies['Vert2Elem']`: Changed from `List[List[int]]` to `Dict[int, Set[int]]`
- ✅ Use public methods instead: `mesh.get_vertex_edges(v)`, `mesh.get_vertex_elements(v)`

See [DOWNSTREAM_MIGRATION_GUIDE.md](docs/DOWNSTREAM_MIGRATION_GUIDE.md) for migration details.

### ✨ Added

#### Phase 1: EdgeMap O(1) Lookup (150x speedup)
- Hash-based edge discovery: O(n log n) vs previous O(n²)
- O(1) edge lookups throughout adjacency building
- Canonical edge form prevents duplicates
- Zero performance regression on small meshes

#### Phase 2: Adjacency Modernization
- `Dict[int, Set[int]]` adjacency structures
- Automatic validation and invariant checking
- Public API: `get_vertex_edges()`, `get_vertex_elements()`
- Type hints for static analysis support

#### Phase 3: Bridge Infrastructure (Domain-Specific APIs)

**MeshAdapterForMADMESHR** (mesh adaptation research)
- `get_element_neighbors(elem_id)` - Find adjacent elements
- `get_element_quality_neighborhood(elem_id)` - Local quality assessment  
- `get_refinement_region(elem_ids)` - Identify refinement zones

**MeshAdapterForADMESH** (adaptive mesh refinement)
- `get_mesh_quality_report()` - Comprehensive quality metrics
- `get_element_angles_summary()` - Angle-based analysis

**MeshAdapterForADMESHDomains** (domain decomposition)
- `get_domain_boundaries()` - Extract boundary vertices
- `get_mesh_connectivity_info()` - Connectivity statistics

#### Phase 4: Advancing-Front API (MADMESHR)

```python
# Get boundary for element placement
boundary = mesh.advancing_front_boundary_edges()

# Add elements incrementally during mesh generation
new_elem_id = mesh.add_advancing_front_element([v1, v2, v3], 'tri')

# Remove residual closure when boundary shrinks
mesh.remove_boundary_loop([edge1, edge2, edge3, edge4])

# Identify bottlenecks for domain splitting
pinches = mesh.pinch_points(width_threshold=0.3)
```

#### Documentation & Release Infrastructure
- `API.md` - Complete 25+ method reference with stability guarantees
- `BENCHMARK.md` - Performance comparison v0.1.1 vs v0.2.0 with real-world impact
- `DOWNSTREAM_MIGRATION_GUIDE.md` - Integration guide for MADMESHR, ADMESH, ADMESH-Domains
- `docs/CHILmesh_Access_Interface.md` - Stability guarantees and usage patterns
- Stable API contract (CAI) through v1.0

### 🚀 Performance

#### Initialization (WNAT_Hagen: 52,774 verts, 98,365 elems)

| Operation | v0.1.1 | v0.2.0 | Speedup |
|-----------|--------|--------|---------|
| Fast init (no layers) | ~3,200s | **3.9s** | **822x** |
| Full init (30 layers) | ~5,400s | **7.7s** | **701x** |
| Quality analysis | ~4,800s | **6.6s** | **727x** |
| **Total** | ~13,400s | **14.3s** | **937x** |

#### Query Performance

| Operation | v0.1.1 | v0.2.0 | Speedup |
|-----------|--------|--------|---------|
| Adjacency lookup | ~2,000μs | **4.0μs** | **500x** |
| Vertex neighbors | ~3,500μs | **0.7μs** | **5,000x** |
| Element neighbors | ~4,500μs | **4.4μs** | **1,022x** |

#### Real-World Impact (MADMESHR)
- Mesh adaptation workflow: 3,800s → 14.7s (**259x speedup**)
- Interactive development instead of hourly batch operations

### 🧪 Test Coverage

- **239 existing tests** (all passing, zero regressions)
- **26 new advancing-front tests** (comprehensive workflow validation)
- **265 total tests**
- Coverage: adjacencies, performance regression, bridge adapters, advancing-front, fort.14 roundtrip

### 📚 Stability Guarantees

All methods in `API.md` guaranteed stable through v1.0:
- ✅ Method signatures won't change
- ✅ Return types won't change
- ✅ Behavior won't change
- ✅ Clear deprecation policy (major version + 2-week notice)

### 🔧 Technical Implementation

**Phase 1 (Issues #11–16):** EdgeMap O(1) lookup, 18 regression tests  
**Phase 2 (Issues #17–23):** Adjacency modernization with validation  
**Phase 3 (Issues #24–31):** Bridge adapters + 44 integration tests  
**Phase 4 (Issues #32+):** Advancing-front API, documentation, release  

### 🔄 Backward Compatibility

✅ Drop-in replacement for v0.1.1  
✅ No code changes needed for downstream projects  
✅ All public APIs identical  
⚠️ Internal structures changed (use public APIs)  
🚀 150x+ performance improvement  

### 📦 Dependencies

No changes from v0.1.1:
- Python 3.10+
- numpy, scipy, matplotlib

### 🏷️ Version Support

| Version | Status | Support |
|---------|--------|---------|
| 0.2.0 | **Current** | Active |
| 0.1.1 | Legacy | Bug fixes until 0.3.0 |

[0.2.0]: https://github.com/domattioli/CHILmesh/releases/tag/v0.2.0

## [0.1.1] — 2026-04-07

This is an honest hotfix release. 0.1.0 shipped with three latent
data-corruption bugs and a README quickstart that no pip user could run.
0.1.1 fixes the bugs, adds regression tests, and ships the example
fixtures inside the wheel so the quickstart actually works.

### Added
- `chilmesh.examples` module with `annulus()`, `donut()`, `block_o()`,
  `structured()` factories that return ready-to-use `CHILmesh` instances
  loaded from bundled `.fort.14` fixtures.
- `chilmesh/data/` package shipping all four fixtures via
  `[tool.setuptools.package-data]`, accessible from any install.
- `chilmesh/utils/__init__.py` so `chilmesh.utils` is a real package
  and `plot_utils.py` always lands in the wheel (B6).
- Regression tests for every B1-B6 bug, parametrized over the four
  fixtures where applicable; layer-disjoint-cover invariants replacing
  the old value-pinned annulus assertions.
- `pyproject.toml` keywords, classifiers, `[project.urls]`, and a
  `[project.optional-dependencies] dev` group.
- GitHub Actions CI workflow (`.github/workflows/ci.yml`) running on
  Python 3.10/3.11/3.12 × ubuntu/macos.
- This `CHANGELOG.md`.

### Fixed
- **B1**: `CHILmesh.write_to_fort14` recursed into itself with the wrong
  arity. Now delegates to the module-level `write_fort14`. Roundtrip
  identity test catches any future regression.
- **B3**: `_elem_type` no longer treats vertex id `0` as a sentinel for
  "missing 4th vertex". Padded triangles are detected by row-pattern.
- **B4**: `_ensure_ccw_orientation` no longer applies the quad
  permutation `[0,3,2,1]` to padded-triangle rows in mixed-element
  meshes. Each CW row is classified before being flipped.
- **B5**: `interior_angles` quad branch guards against zero-length edges
  with the same `+1e-12` epsilon already used by the triangle branch and
  takes the real part of the result, mirroring MATLAB's `real(acosd(...))`.
- `read_from_fort14` accepts node indices written in float form
  (e.g. `1.000000`) so legacy fixtures load.

### Changed
- Bumped `requires-python` to `>=3.10` (drops EOL 3.8/3.9).
- `requirements.txt` reduced from a 60-line MADMESHR ML environment
  freeze to the three actual runtime dependencies (`numpy`, `scipy`,
  `matplotlib`).
- README quickstart rewritten around `chilmesh.examples.annulus()`;
  the broken `/kaggle/working/...` path is gone, the install command
  now includes `-r`, the stale "pip installation" TODO is removed,
  an alpha-status banner has been added, and the MADMESHR "Future Work"
  badge has been reframed as a "See also / Downstream projects"
  section that explicitly notes CHILmesh is agnostic to MADMESHR.

### Investigated, no fix needed
- **B2** (`signed_area` quad shoelace): the audit reported the formula
  used 4 cross terms instead of 8. On inspection the existing 4-row
  formula expands to the correct 8 terms, and the new random-quad
  cross-check test against an independent shoelace passes. Tests
  retained as a hardening guard.

### Known issues / deferred to 0.1.2
- The O(n²) Python loops in `_build_elem2edge` / `_build_edge2elem`
  make Block_O (~5,200 elements) slow to load. Vectorisation lands
  in 0.1.2.
- Coverage gating, Hypothesis property tests, mixed-element smoothing
  characterisation, Windows + Python 3.13 CI matrix.

[0.1.1]: https://github.com/domattioli/CHILmesh/releases/tag/v0.1.1
