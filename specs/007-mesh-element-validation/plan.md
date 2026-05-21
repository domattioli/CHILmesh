# Implementation Plan: Mesh Element Validity Test Suite

**Branch**: `claude/mesh-quad-triangle-spec-IIvKw` | **Date**: 2026-05-21 | **Spec**: [`spec.md`](./spec.md)
**Input**: Feature specification from `/specs/007-mesh-element-validation/spec.md`

## Summary

Ship a regression-gate test suite that verifies a `CHILmesh` is quad-dominant with triangles confined to layer 0 (boundary), all elements planar 2D and non-self-intersecting, no element-element interior overlap, and no edge-edge crossings between non-adjacent elements. Degenerate quads (collinear, zero-area, padded-triangle, duplicate vertex) remain accepted and are recorded as informational notes.

Scope is **test-suite-only** at this stage (Clarify Q6). All code lives under `tests/_validity/`; no `src/chilmesh/` change. A discussion issue (#142) tracks long-term promotion to `chilmesh.validate`.

Technical approach:

- Single entry point `validate_mesh_elements(mesh) -> MeshValidityReport` in `tests/_validity/validator.py`.
- Two-tier predicates: cheap per-element checks (arity, layer membership, self-intersection) in pure numpy; pair-wise checks (interior overlap, edge crossing) behind a uniform-grid broadphase to keep block_o under budget.
- Synthetic negative fixtures built by post-mutating `connectivity_list` after `CHILmesh.__init__` (the constructor's degeneracy fallback would otherwise repair planted violations).
- Pytest UX: one test per fixture asserting `report.ok`; failure message aggregates all violation categories.

## Technical Context

**Language/Version**: Python 3.10+
**Primary Dependencies**: numpy (already a CHILmesh dep), pytest (already a dev dep). No new deps.
**Storage**: N/A (in-memory only).
**Testing**: pytest, parametrized over the four built-in fixtures (`annulus`, `donut`, `block_o`, `structured`) + five synthetic negative fixtures.
**Target Platform**: Linux / macOS / Windows (per CHILmesh CI matrix; currently Linux py3.11 only on this PR).
**Project Type**: library (CHILmesh) — extending its test surface, not its public API.
**Performance Goals**: Total runtime ≤ 60 s across the four built-in fixtures (block_o ≤ 45 s; others ≤ 5 s each). Broadphase MUST reduce narrowphase pair calls by ≥100× vs `O(n²)` on block_o (SC-007).
**Constraints**: No new public API surface (Clarify Q6). No new deps. Constitution principle V — coordinate-system agnostic (tolerances are bbox-relative).
**Scale/Scope**: ~5,200 elements on block_o (largest fixture). One new test module, one helper package, ≤ ~600 lines of new code total.

## Constitution Check

| Principle | Gate | Verdict |
|-----------|------|---------|
| I — Library-first | New module independently testable; clear surface | ✅ `tests/_validity/` is self-contained, importable, depends only on `CHILmesh` public surface plus `_skeletonize()` (existing) |
| II — Mesh immutability | Validator MUST NOT mutate the mesh | ✅ FR-001 forbids mutation; FR-007 mutation of layer cache is via existing `_skeletonize()` which CHILmesh already caches |
| III — Test-first + regression-gated | Tests written before implementation | ✅ Spec defines acceptance scenarios; T-tasks order tests before validator code |
| IV — Geometric correctness over performance | No floating-point shortcuts that mask bowtie | ✅ FR-009 uses segment-intersection predicate, not signed-area proxy (Clarify Q3) |
| V — Coordinate-system agnostic | No hard-coded coord scale | ✅ FR-013 tolerance is bbox-relative (Clarify Q7) |
| VI — Format pluralism | N/A — operates on parsed mesh | ✅ |
| VII — Public API stability | No public API changes | ✅ Confined to `tests/_validity/` (Clarify Q6) |

**No violations. No Complexity Tracking entries needed.**

## Phase 0: Research

**Output:** `research.md` — design decisions captured. See that doc for full detail. Key bullets:

- **Bowtie predicate.** Use the standard "two diagonals of the polygon" formulation: for vertices `v0, v1, v2, v3`, the quad is non-self-intersecting iff edges `(v0,v1)` and `(v2,v3)` do NOT properly cross AND edges `(v1,v2)` and `(v3,v0)` do NOT properly cross. Robust orientation predicate via numpy 2D cross-product sign with tol-aware tie-breaking.
- **Point-in-polygon.** Winding-number test (handles non-convex quads; convex ray-casting also works but winding is one branch simpler for quads + tris uniformly). Robust on degenerate quads because we treat boundary as "outside" (strict interior).
- **Broadphase.** Uniform grid hash sized by mean element bbox diagonal × 2. For block_o (~5,200 elems, square-ish domain), cell side ≈ √(area/n_elems) × 2 yields ≈ 1,200 cells, ~4-5 elements per cell average. Narrowphase pair calls drop from `O(n²)≈2.7e7` to `O(n * avg_per_cell)≈2.6e4`, ≥1000× reduction (SC-007 ≥100× passes with margin).
- **Edge sharing detection.** Build a `frozenset({va, vb})`-keyed map from vertex pair → element list (one pass over `connectivity_list`). Two elements are "edge-adjacent" iff they share an entry; "vertex-adjacent" iff they share ≥1 vertex but no edge entry. Pairs that are edge-adjacent are exempt from the edge-crossing test (their shared edge trivially "crosses" otherwise).
- **Padded-triangle classification.** A row `[a, b, c, d]` is a triangle iff `d == -1` OR `d ∈ {a, b, c}`. Otherwise quad.
- **Tolerance flow.** Spec FR-013: compute `bbox_diag = ||(max_x-min_x, max_y-min_y)||₂`; `tol_effective = max(1e-15, 1e-12 * bbox_diag)`. Caller `tol` override propagates everywhere.
- **Synthetic fixture construction.** `CHILmesh.__init__` triggers degeneracy fallback that would repair bowties. Strategy: build a *valid* mesh, then post-mutate `mesh.connectivity_list[row_id] = bad_row` and clear cached adjacency / layers so subsequent `validate_mesh_elements` recomputes against the corrupted state. Helper `corrupt_to(mesh, row_id, new_row)` lives in `tests/_validity/fixtures.py`.

## Phase 1: Design

**Outputs:**

- `data-model.md` — types for `Violation`, `InformationalNote`, `MeshValidityReport`.
- `contracts/validator.py` (signature only) — frozen `validate_mesh_elements` signature so tests can be written before implementation.
- `quickstart.md` — short reproducible example a maintainer can paste into a notebook to dry-run the validator.

### Module layout (under `tests/_validity/`)

```text
tests/
├── _validity/
│   ├── __init__.py          # exports validate_mesh_elements, types
│   ├── types.py             # MeshValidityReport, Violation, InformationalNote (frozen dataclasses)
│   ├── predicates.py        # segment_proper_cross, point_in_polygon, classify_element, bbox_diag
│   ├── broadphase.py        # UniformGridIndex (build + query bbox-overlap pairs)
│   ├── validator.py         # validate_mesh_elements main entry; orchestrates predicates
│   ├── fixtures.py          # bowtie_quad_mesh, interior_triangle_mesh, etc. + corrupt_to helper
│   └── README.md            # short doc pointing at spec.md
└── test_mesh_element_validity.py  # pytest entry; parametrized over (built-in + synthetic) fixtures
```

### Surface (frozen for test writing)

```python
# tests/_validity/validator.py
from chilmesh import CHILmesh
from tests._validity.types import MeshValidityReport

def validate_mesh_elements(
    mesh: CHILmesh,
    *,
    tol: float | None = None,
) -> MeshValidityReport:
    """Verify mesh element validity per spec 007.

    Args:
        mesh: CHILmesh instance.
        tol: Override for the absolute tolerance. If None, uses
             1e-12 * bbox_diag (floored at 1e-15) per FR-013.

    Returns:
        MeshValidityReport with ok, violations, notes, n_elems_checked, runtime_s.
    """
```

```python
# tests/_validity/types.py
from dataclasses import dataclass, field

@dataclass(frozen=True)
class Violation:
    category: str
    element_ids: tuple[int, ...]
    edge_ids: tuple[int, ...] | None
    detail: str

@dataclass(frozen=True)
class InformationalNote:
    category: str
    element_ids: tuple[int, ...]
    detail: str

@dataclass(frozen=True)
class MeshValidityReport:
    ok: bool
    violations: tuple[Violation, ...]
    notes: tuple[InformationalNote, ...]
    n_elems_checked: int
    runtime_s: float
```

### Project structure (repository root)

```text
src/chilmesh/         # UNCHANGED
specs/007-mesh-element-validation/
├── spec.md
├── plan.md              # this file
├── research.md          # phase 0 output
├── data-model.md        # phase 1 output
├── quickstart.md        # phase 1 output
├── contracts/
│   └── validator.py     # frozen signature
├── tasks.md             # phase 2 output (/speckit.tasks)
└── analyze.md           # phase 2.5 output (/speckit.analyze)
tests/
├── _validity/           # NEW — helper package
└── test_mesh_element_validity.py   # NEW
```

**Structure Decision**: Option 1 (single-project library) — no src/ change. New code is confined to a `tests/_validity/` helper package + one pytest module.

## Phase 2 → Phase 5 milestones

- **Phase 2 (Setup)**: scaffold `tests/_validity/` directory + `__init__.py`, write `types.py`. ~30 LOC.
- **Phase 3 (Foundational)**: predicates (`predicates.py`), broadphase (`broadphase.py`), bbox helpers. ~200 LOC. Unit-tested via private cases inside each module.
- **Phase 4 (US1 — element-type composition)**: `classify_element`, FR-003/004/005/006/007 checks. ~120 LOC. Parametrize against built-in fixtures + `interior_triangle_mesh` + `pentagon_mesh`.
- **Phase 5 (US2 — geometric validity)**: FR-008/009/011/012. ~180 LOC. Parametrize + add `bowtie_quad_mesh`, `overlapping_quads_mesh`, `edge_crossing_mesh`.
- **Phase 6 (US3 — degenerate-quad first-class)**: ensure `DEGENERATE_QUAD_DUPLICATE_VERTEX` and `DEGENERATE_ZERO_AREA` paths emit notes not violations; cross-test against existing `tests/test_degeneracy.py`. ~40 LOC.
- **Phase 7 (US4 — pytest UX + budget)**: `test_mesh_element_validity.py` parametrization; failure-message aggregation (FR-016/017); SC-004 runtime assertion via `pytest --durations=10` + an explicit `assert report.runtime_s < budget`. ~80 LOC.
- **Phase 8 (Polish)**: docstrings, `tests/_validity/README.md` pointing at spec, `CHANGELOG.md` "Internal" entry (not public API).

## Complexity Tracking

No constitution violations. Section intentionally empty.

## Risk register

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| block_o runtime > 45 s | medium | Profile broadphase early (Phase 3 exit gate); fall back to `scipy.spatial.cKDTree` if uniform grid loses cache locality (scipy already in dev deps) |
| Constructor degeneracy fallback repairs synthetic bowties | medium | Post-mutate connectivity AFTER construct + cache-bust (`mesh._invalidate_cache()` or manual clear); document in `fixtures.py` |
| `_skeletonize()` cost on first call adds ~14 s on block_o (FR-007) | low | One-time; cached on `mesh.layers`. Acceptable because validator is not on the hot path |
| `mesh.boundary_vertices()` semantics differ from expected (e.g., returns int IDs vs ndarray) | low | Verified in Phase 0 by reading `tests/test_invariants.py:74` — it returns iterable of int. Wrap in `set(int(x) for x in ...)` defensively |
| Floating-point ties in segment-cross predicate flap near tol boundary | medium | Use a symmetric robust orientation predicate (sign of 2D cross) and treat `|cross| <= tol_effective` as collinear (no crossing). Documented in `predicates.py` |

## Re-check after Phase 1

- All FRs (FR-001 → FR-017) have a home in the module layout above. ✅
- All SC items (SC-001 → SC-008) map to a test case in `test_mesh_element_validity.py` or the issue-tracker entry (#142). ✅
- Constitution principles still hold (no public API growth, no perf regression risk on existing code). ✅
