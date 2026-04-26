# Tasks: ADMESH-Domains I/O Parity

**Phase 2 output for `/speckit-tasks`**
**Branch**: `planning-optimize_modernize` | **Date**: 2026-04-26 | **Plan**: [plan.md](plan.md)

---

## Overview

**Feature**: CHILmesh ↔ ADMESH-Domains I/O Parity
**Total Tasks**: 42 | **Phases**: 6 | **Suggested MVP**: Phases 1–2 (Setup + Foundational fix)

**Task Dependencies**:
```
Phase 1 (Setup) → Phase 2 (Foundational) → Phase 3 (US1) → Phase 4 (US2) → Phase 5 (US3) → Phase 6 (Polish)
                        ↓
            Unblocks all user stories
```

**Parallel Opportunities**:
- Tasks T003–T007 (Tests for US1): can run in parallel once T002 completes
- Tasks T008–T012 (Metadata tasks for US2): can run in parallel with T007
- Tasks T013–T018 (Entry point for US1): can run in parallel with US2 tasks

---

## Phase 1: Setup

*Create test fixtures and register in project metadata.*

- [ ] T001 Create synthetic quad test fixture at `src/chilmesh/data/quad_2x2.fort.14` (9 nodes, 4 quads, exact format in `contracts/public-api.md`)
- [ ] T002 Register `quad_2x2()` factory in `src/chilmesh/examples.py` and add `"quad_2x2"` to `conftest.py` `FIXTURE_NAMES`

---

## Phase 2: Foundational — Group A (Quad I/O)

*Critical fix: reader and writer must handle quad/mixed elements. This unblocks all user stories.*

### Understanding Current State
- [ ] T003 Inspect `read_from_fort14()` at `src/chilmesh/CHILmesh.py:662–699` and identify the `ValueError` on line 692–693
- [ ] T004 Inspect `write_fort14()` at `src/chilmesh/CHILmesh.py:969–995` and identify hardcoded `"3"` on line 991

### Reader Implementation
- [ ] T005 [P] Modify `read_from_fort14()` signature to add `compute_layers: bool = True` kwarg
- [ ] T006 [P] Implement two-pass read logic in `read_from_fort14()`: scan for max `num_nodes`, allocate `elements` array with correct width (3 or 4)
- [ ] T007 [P] Apply padded-triangle convention in `read_from_fort14()`: triangles stored as `[v0, v1, v2, v0]` in 4-col arrays
- [ ] T008 [P] Add error handling in `read_from_fort14()`: raise `ValueError(f"Unsupported element type: element {id} has {n} nodes (only 3 or 4 supported).")` for 5+ nodes
- [ ] T009 [P] Forward `compute_layers` kwarg from `read_from_fort14()` to `CHILmesh.__init__()`

### Writer Implementation
- [ ] T010 [P] Modify `write_fort14()` to detect quad vs triangle per row: check if 4-col and `row[3] != row[0]`
- [ ] T011 [P] Write `"i 4 n1 n2 n3 n4"` for quads and `"i 3 n1 n2 n3"` for triangles in `write_fort14()`
- [ ] T012 [P] Preserve existing 3-column-only array path in `write_fort14()` (backward compat)

### Integration & Tests
- [ ] T013 Create `tests/test_fort14_quad_roundtrip.py`: load `quad_2x2`, assert `n_verts==9`, `n_elems==4`, `type=="Quadrilateral"`, write to tmpfile, reload, assert counts match
- [ ] T014 [P] Add quad roundtrip tests for all quad-supporting fixtures (Block_O, annulus if mixed, etc.)
- [ ] T015 Ensure `quad_2x2` passes all parametrised tests in `conftest.py` (e.g., `test_ccw`, `test_elem_type`)

---

## Phase 3: User Story 1 — Load Real Coastal Mesh (P1)

*Requirement*: Researchers can load any ADMESH-Domains mesh (including quad/mixed elements) into CHILmesh and analyse it.
*Depends on*: Phase 2 (quad I/O fixed)
**Independent Test Criteria**:
- Load 1 quad and 1 mixed-element `.fort.14` from ADMESH-Domains catalog
- Verify `n_verts`, `n_elems`, `type` are correct
- Call `elem_quality()` and `interior_angles()` — no crashes, return valid results
- Save back to `.fort.14` → reload → verify roundtrip fidelity

### Entry Point Implementation
- [ ] T016 [P] Add `CHILmesh.from_admesh_domain()` classmethod signature in `src/chilmesh/CHILmesh.py` (after `read_from_fort14()`)
- [ ] T017 [P] Implement duck-typed record access in `from_admesh_domain()`: use `getattr(record, "filename", None)` and `getattr(record, "type", None)`
- [ ] T018 [P] Implement routing logic in `from_admesh_domain()`: `"SMS_2DM"` → `read_from_2dm()`, else → `read_from_fort14()`
- [ ] T019 [P] Add file-not-found guard in `from_admesh_domain()`: check `Path(filename).exists()`, raise `FileNotFoundError("File not found: {path}. If using ADMESH-Domains, call mesh_record.load() first.")`
- [ ] T020 [P] Handle unknown `record.type` in `from_admesh_domain()`: emit `warnings.warn(f"Unrecognised mesh type '{t}', falling back to ADCIRC reader.")` then route to `read_from_fort14()`
- [ ] T021 [P] Forward `compute_layers` kwarg through `from_admesh_domain()` to the chosen reader

### Testing US1
- [ ] T022 Create `tests/test_from_admesh_domain.py` with mock ADMESH-Domains record (SimpleNamespace with `filename`, `type`)
- [ ] T023 [P] Test ADCIRC routing: `record.type="ADCIRC"` → loads via `read_from_fort14()`
- [ ] T024 [P] Test ADCIRC_GRD routing: `record.type="ADCIRC_GRD"` → routes to `read_from_fort14()`
- [ ] T025 [P] Test unknown type: `record.type="UNKNOWN"` → emits warning, falls back to `read_from_fort14()`
- [ ] T026 [P] Test missing file: `record.filename=/nonexistent` → raises `FileNotFoundError` with guidance message
- [ ] T027 [P] Test `compute_layers=False` forwarding: mesh initialised with no layers
- [ ] T028 Create integration test: load a real quad/mixed mesh from ADMESH-Domains (if available locally), run `elem_quality()` and `interior_angles()`, verify no crashes

---

## Phase 4: User Story 2 — Bulk Metadata Queries (P1)

*Requirement*: Developers can load 10 catalog meshes with layers disabled, each in <2s, and call `admesh_metadata()` to get accurate node/element/type/bbox.
*Depends on*: Phase 2 (quad reader), Phase 3 (entry point)
**Independent Test Criteria**:
- Load each fixture with `compute_layers=False` → init completes in <2s
- Call `admesh_metadata()` on each → `node_count` and `element_count` match file
- `element_type` is one of the three canonical strings
- `bounding_box` keys exist and have sensible values

### Fast Initialisation Implementation
- [ ] T029 [P] Add `compute_layers: bool = True` kwarg to `CHILmesh.__init__()` signature at `src/chilmesh/CHILmesh.py:59`
- [ ] T030 [P] Add `compute_layers` parameter to `_initialize_mesh()` at `src/chilmesh/CHILmesh.py:101`
- [ ] T031 [P] Implement conditional skip in `_initialize_mesh()`: when `compute_layers=False`, skip `_build_adjacencies()` and `_skeletonize()`
- [ ] T032 [P] Ensure `self.n_layers = 0` and `self.adjacencies = {}` when `compute_layers=False`

### Element Type Classification
- [ ] T033 [P] Set `self.type` in `_initialize_mesh()` after `_ensure_ccw_orientation()` by calling `_elem_type()` and deriving canonical string
- [ ] T034 [P] Implement `self.type` classification logic: all-tri → `"Triangular"`, all-quad → `"Quadrilateral"`, mixed → `"Mixed-Element"`

### Metadata Method Implementation
- [ ] T035 [P] Implement `CHILmesh.admesh_metadata()` method in `src/chilmesh/CHILmesh.py` (after `get_layer()`)
- [ ] T036 [P] Populate metadata dict with `node_count` (from `n_verts`), `element_count` (from `n_elems`), `element_type` (from `self.type`)
- [ ] T037 [P] Compute and populate `bounding_box` dict in `admesh_metadata()`: `min_x/max_x/min_y/max_y` from `points[:, 0]` and `points[:, 1]` extrema
- [ ] T038 [P] Ensure `admesh_metadata()` is safe to call with `compute_layers=False` (no dependencies on layers)

### Testing US2
- [ ] T039 Create `tests/test_admesh_metadata.py`
- [ ] T040 [P] Test fast init: load each fixture with `compute_layers=False`, assert init time <2s (use `time.perf_counter`), assert `n_layers==0`
- [ ] T041 [P] Test metadata accuracy: call `admesh_metadata()` on each fixture; assert `node_count == n_verts`, `element_count == n_elems`, `element_type` is canonical string
- [ ] T042 [P] Test bounding_box: verify `min_x <= max_x`, `min_y <= max_y`, extrema match actual point coordinates
- [ ] T043 [P] Test metadata with `compute_layers=False`: confirm metadata still available (no layers required)

---

## Phase 5: User Story 3 — Validate Metadata Before Contributing (P3)

*Requirement*: Contributors can load a mesh, call `admesh_metadata()`, and compare returned values against manifest entries to catch mismatches.
*Depends on*: Phase 2 (quad reader), Phase 4 (metadata method)
**Independent Test Criteria**:
- Load a mesh with wrong `node_count` in metadata (manually crafted test case)
- Call `admesh_metadata()` → `node_count` mismatch is immediately visible
- Metadata dict is complete (all 4 schema fields present)

### Guard Improvements
- [ ] T044 Improve `get_layer()` error message in `src/chilmesh/CHILmesh.py:649`: prepend guard `if self.n_layers == 0: raise RuntimeError("Layers not computed. Re-initialise with compute_layers=True.")`

### Testing US3
- [ ] T045 Create `tests/test_metadata_validation.py`
- [ ] T046 [P] Test metadata completeness: call `admesh_metadata()`, assert dict has all 4 keys (`node_count`, `element_count`, `element_type`, `bounding_box`)
- [ ] T047 [P] Test `get_layer()` guard: initialise with `compute_layers=False`, call `get_layer(0)`, assert `RuntimeError` with guidance message
- [ ] T048 [P] Test contributor workflow: manually construct a mesh record with wrong `node_count`, load via `from_admesh_domain()`, call `admesh_metadata()`, verify mismatch detectable

---

## Phase 6: Polish & Cross-Cutting Concerns

*Complete remaining features and error handling.*

### SMS .2dm Reader (Group E)
- [ ] T049 [P] Implement `CHILmesh.read_from_2dm()` staticmethod in `src/chilmesh/CHILmesh.py` (after `read_from_fort14`)
- [ ] T050 [P] Single-pass `.2dm` parser: keyword dispatch for `ND`, `E3T`, `E4Q`, comment/skip lines
- [ ] T051 [P] Handle 1-based node IDs in `.2dm` file; subtract 1 for 0-based storage
- [ ] T052 [P] Apply padded-triangle convention in `read_from_2dm()` for mixed-element meshes
- [ ] T053 [P] Add file-not-found and unsupported-element error handling in `read_from_2dm()`
- [ ] T054 [P] Forward `compute_layers` kwarg from `read_from_2dm()` to `CHILmesh.__init__()`

### 2dm Reader Tests
- [ ] T055 Create `tests/test_2dm_reader.py`
- [ ] T056 [P] Write minimal `.2dm` strings to tmpfiles; test triangle-only, quad-only, mixed-element parsing
- [ ] T057 [P] Assert `n_verts`, `n_elems`, `type` correctness
- [ ] T058 [P] Test file-not-found and unsupported-element errors

### Documentation
- [ ] T059 Update docstrings in `read_from_fort14()`, `write_fort14()`, `__init__()` to document quad/mixed support and `compute_layers` kwarg
- [ ] T060 Add docstring to `from_admesh_domain()` noting duck-typed interface and no `admesh-domains` import required
- [ ] T061 Add docstring to `admesh_metadata()` documenting return dict schema and Cartesian coordinate convention
- [ ] T062 Add docstring to `read_from_2dm()` documenting SMS format support and padded-triangle convention

### Final Validation
- [ ] T063 Run full test suite (`pytest tests/` -v): all existing tests pass (zero regressions — SC-005)
- [ ] T064 [P] Run parametrised fixture tests on `quad_2x2` (ccw, elem_type, interior_angles, elem_quality)
- [ ] T065 [P] Verify all 4 original fixtures still work with new code (`compute_layers=True` path unchanged)

---

## Implementation Strategy

**MVP Scope** (Phases 1–2 + core US1):
- T001–T015: Setup + quad I/O fix
- T016–T028: Entry point + US1 integration test
- **Result**: Researchers can load any ADMESH-Domains mesh and analyse it.
- **Effort**: ~20 tasks, ~40% of total

**Phase-2 Extension** (add US2):
- T029–T043: Fast init + metadata
- **Result**: Catalog queries fast (<2s), metadata accurate (SC-004)
- **Incremental effort**: ~15 tasks, ~35% of total

**Full Scope** (all phases):
- T044–T065: US3 + polish + 2dm reader
- **Result**: Complete ADMESH-Domains integration, SMS support, full error handling

**Parallelisation**:
- Setup (T001–T002): Serial prerequisite
- Foundational (T003–T015): T005–T009 [P] (reader), T010–T012 [P] (writer) can run in parallel
- US1 (T016–T028): T016–T021 [P] can run in parallel with foundational
- US2 (T029–T043): T029–T038 [P] can run in parallel; T039–T043 [P] tests can run in parallel once implementations complete
- US3 (T044–T048): Serialized (depends on US2 completion)
- Polish (T049–T065): T049–T058 [P] 2dm reader, T059–T062 [P] docs, T063–T065 [P] tests can run in parallel

**Estimated velocity** (if parallelised):
- Critical path: Phases 1–2 (Setup → Quad fix) → Phase 3 (Entry point) → Phase 4 (Metadata) → Phase 5 (Guard) → Phase 6 (Polish)
- Serial: ~15 phases × ~3 days/phase ≈ 45 days (conservative)
- With parallelisation: ~25 days (50% faster, mainly in Groups B–E and testing phases)
