# Changelog

All notable changes to this project will be documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and the project adheres to [Semantic Versioning](https://semver.org/).

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
