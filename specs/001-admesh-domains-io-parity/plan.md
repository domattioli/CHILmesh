# Implementation Plan: ADMESH-Domains I/O Parity

**Branch**: `planning-optimize_modernize` | **Date**: 2026-04-26 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `specs/001-admesh-domains-io-parity/spec.md`

## Summary

CHILmesh must load any mesh from the ADMESH-Domains catalog, save it back without data loss, and expose metadata fields compatible with the ADMESH-Domains schema. Two critical bugs block this today: the `.fort.14` reader crashes on quad elements, and the writer also crashes on quad elements. This plan fixes both atomically, adds fast initialisation, a metadata method, an ADMESH-Domains entry point, and an SMS `.2dm` reader — enabling CHILmesh to serve as the computational engine behind the catalog.

## Technical Context

**Language/Version**: Python 3.10, 3.11, 3.12 (all three tested in CI)
**Primary Dependencies**: numpy>=1.23, scipy>=1.10, matplotlib>=3.6
**Storage**: On-disk mesh files (`.fort.14`, `.14`, `.2dm`); no database
**Testing**: pytest>=7, pytest-cov
**Target Platform**: All platforms (pure Python + numpy/scipy)
**Project Type**: Library (`pip install chilmesh`)
**Performance Goals**: `compute_layers=False` init in <2s for meshes up to ~15,000 elements (current: ~30s for Block_O ~5,200 elements)
**Constraints**: Zero regressions on existing 4 triangular fixtures; no new mandatory dependencies
**Scale/Scope**: Single library module (`src/chilmesh/CHILmesh.py`); 1,000 LOC change estimate

## Constitution Check

*Constitution is a blank template — no project-specific gates defined. Standard quality gates apply:*

- [x] No new mandatory dependencies introduced
- [x] All existing tests continue to pass (zero regression requirement — SC-005)
- [x] Each change is independently testable
- [x] Reader and writer bugs fixed atomically (C-1 from clarify phase)

## Project Structure

### Documentation (this feature)

```text
specs/001-admesh-domains-io-parity/
├── plan.md              ← this file
├── research.md          ← Phase 0 output
├── data-model.md        ← Phase 1 output
├── contracts/
│   └── public-api.md    ← Phase 1 output
├── checklists/
│   └── requirements.md  ← quality checklist (all pass)
├── spec.md
└── tasks.md             ← /speckit-tasks output (not yet created)
```

### Source Code (repository root)

```text
src/chilmesh/
├── CHILmesh.py          ← primary change target
└── data/
    └── quad_2x2.fort.14 ← new synthetic fixture (4 quads, 9 nodes)

tests/
├── conftest.py          ← add quad fixture to FIXTURE_NAMES
├── test_fort14_quad_roundtrip.py   ← new: FR-001, FR-002, SC-002
├── test_fast_init.py               ← new: FR-005, SC-003
├── test_admesh_metadata.py         ← new: FR-004, SC-004
├── test_from_admesh_domain.py      ← new: FR-006, FR-007
└── test_2dm_reader.py              ← new: FR-003
```

## Implementation Groups

### Group A — Quad I/O (P0, atomic, blocks everything else)
Maps to: FR-001, FR-002, SC-001, SC-002

**A-1: Add synthetic quad fixture**
- Create `src/chilmesh/data/quad_2x2.fort.14` (exact content in `contracts/public-api.md`)
- Register `"quad_2x2"` example in `src/chilmesh/examples.py` (`quad_2x2()` factory)
- Register in `pyproject.toml` package-data (`.fort.14` glob already covers it)

**A-2: Fix `read_from_fort14()` for quad/mixed elements**
- File: `src/chilmesh/CHILmesh.py:688–699`
- Two-pass: scan element block to find max `num_nodes`; allocate `elements` with correct column count
- Apply padded-triangle convention for mixed meshes (col 4 = col 1 for triangles)
- Raise `ValueError(f"Unsupported element type: element {id} has {n} nodes (only 3 or 4 supported).")` for 5+ nodes
- Add `compute_layers: bool = True` kwarg and forward to `CHILmesh.__init__`

**A-3: Fix `write_fort14()` for quad/mixed elements**
- File: `src/chilmesh/CHILmesh.py:969–995`
- Per-row: detect quad vs triangle (4-col array + `row[3] != row[0]`)
- Write `"i 4 n1 n2 n3 n4"` for quads, `"i 3 n1 n2 n3"` for triangles
- Keep existing triangle-only path for 3-column arrays

**A-4: Tests**
- `tests/test_fort14_quad_roundtrip.py`: load `quad_2x2`, assert `n_verts==9`, `n_elems==4`, `type=="Quadrilateral"`, write to temp file, reload, assert same counts
- Add `"quad_2x2"` to `conftest.py` `FIXTURE_NAMES` so it participates in all parametrised tests

---

### Group B — Fast Initialisation (P1)
Maps to: FR-005, SC-003
Depends on: Group A (quad fixture must load correctly first)

**B-1: Add `compute_layers` kwarg to `__init__` and `_initialize_mesh`**
- File: `src/chilmesh/CHILmesh.py:59, 101`
- `compute_layers: bool = True` parameter; when `False`, skip `_build_adjacencies()` and `_skeletonize()`
- `self.n_layers` stays `0`; `self.adjacencies` stays `{}`

**B-2: Improve `get_layer()` guard**
- File: `src/chilmesh/CHILmesh.py:649`
- Prepend check: `if self.n_layers == 0: raise RuntimeError("Layers not computed. Re-initialise with compute_layers=True.")`

**B-3: Set `self.type` in `_initialize_mesh()`**
- File: `src/chilmesh/CHILmesh.py:101`
- After CCW orientation: call `_elem_type()`, derive string, set `self.type`
- Logic: all-tri → `"Triangular"`, all-quad → `"Quadrilateral"`, mixed → `"Mixed-Element"`

**B-4: Tests**
- `tests/test_fast_init.py`: load each fixture with `compute_layers=False`, assert init <2s (use `time.perf_counter`), assert `n_layers==0`, assert `get_layer(0)` raises `RuntimeError`

---

### Group C — Metadata Method (P1)
Maps to: FR-004, SC-004
Depends on: Group B (needs `self.type` set)

**C-1: Implement `admesh_metadata()`**
- File: `src/chilmesh/CHILmesh.py` (add after `get_layer`)
- Returns dict with `node_count`, `element_count`, `element_type`, `bounding_box` (see `contracts/public-api.md`)
- `element_type` from `self.type` (set in B-3)
- `bounding_box` from `points[:, 0].min()` / `.max()`, `points[:, 1].min()` / `.max()`
- Safe to call when `compute_layers=False`

**C-2: Tests**
- `tests/test_admesh_metadata.py`: call `admesh_metadata()` on each fixture; assert `node_count == n_verts`, `element_count == n_elems`, `bounding_box` keys present, `element_type` is canonical string
- Test with `compute_layers=False` to confirm metadata still available

---

### Group D — ADMESH-Domains Entry Point (P1)
Maps to: FR-006, FR-007 (file-not-found case)
Depends on: Groups A, B, C

**D-1: Implement `CHILmesh.from_admesh_domain()`**
- File: `src/chilmesh/CHILmesh.py` (add as classmethod after `read_from_fort14`)
- Check `Path(filename).exists()`; raise `FileNotFoundError` with guidance message
- Route by `record.type`: `"SMS_2DM"` → `read_from_2dm(filename, compute_layers=compute_layers)` (from Group E); else → `read_from_fort14(filename, compute_layers=compute_layers)`
- Unknown type: `warnings.warn(f"Unrecognised mesh type '{t}', falling back to ADCIRC reader.")`

**D-2: Tests**
- `tests/test_from_admesh_domain.py`:
  - Mock record (simple dataclass or SimpleNamespace) with `filename`, `type`
  - ADCIRC routing test: record.type="ADCIRC", loads correctly
  - ADCIRC_GRD routing test: routes to fort14 reader
  - SMS_2DM routing test: routes to 2dm reader
  - Missing file test: raises `FileNotFoundError` with guidance
  - Unknown type test: emits warning, still loads
  - `compute_layers=False` forwarded correctly

---

### Group E — SMS `.2dm` Reader (P2)
Maps to: FR-003
Depends on: Group A (padded-triangle convention established)

**E-1: Implement `CHILmesh.read_from_2dm()`**
- File: `src/chilmesh/CHILmesh.py` (add as staticmethod after `read_from_fort14`)
- Single-pass keyword dispatch: `ND` → node, `E3T` → triangle, `E4Q` → quad, `#` / `MESH2D` → skip
- 1-based node IDs in file; subtract 1 for 0-based storage
- Mixed meshes → 4-column padded-triangle connectivity
- `compute_layers` kwarg forwarded to `CHILmesh.__init__`
- Raise `ValueError` for unsupported element keywords (5-node or unknown)

**E-2: Tests**
- `tests/test_2dm_reader.py`: write a minimal `.2dm` string to a tmp file; assert correct `n_verts`, `n_elems`, `type`; test triangle-only, quad-only, mixed; test FileNotFoundError; test unsupported element ValueError

---

## Delivery Order

```
A (quad I/O fix) → B (fast init) → C (metadata) → D (entry point) → E (2dm reader)
```

Groups A–D are required for SC-001 through SC-004. Group E (2dm) completes FR-003 and can be delivered independently after D.

## Complexity Tracking

No constitution violations. No new abstractions beyond the feature spec.
