# Implementation Plan: ADMESH Warm-Start Truss Optimization

**Branch**: `main` (CHILmesh single-branch policy; planning-optimize_modernize is in a divergent state) | **Date**: 2026-05-02 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/005-admesh-warm-start-truss/spec.md`

## Summary

Build a generic adapter that feeds an existing triangulation (points + triangles + boundary) into ADMESH's distmesh truss optimizer while preserving the input boundary bit-exactly. The adapter exposes two layered entry points: a low-level array-based form and a high-level CHILmesh wrapper. The implementation **vendors the inner truss loop** of `admesh.distmesh.distmesh2d` (exactly as written upstream at the pinned commit) into a private CHILmesh helper, with the only modification being that the initial point distribution comes from the caller's existing mesh instead of `_initial_distribution` + `_rejection_method`. This lets us bypass ADMESH's broken `routine.py` (missing `MeshOutput`) entirely while still using ADMESH's actual truss algorithm — every line of the truss loop is byte-identical to upstream. Once ADMESH ships an `initial_points` parameter or a `distmesh2d_warmstart` entry point (filed as issue ADMESH-B), the vendored copy is removed and we call ADMESH directly.

The user-facing demo is the existing `generate_4row_admesh.py` restructured per Q2=d: Row 1 = raw annulus, Row 2 = warm-start of Row 1, Row 3 = FEM smoother applied to Row 2, Row 4 = right-isoceles smoother applied to Row 2. The fresh-ADMESH-from-bbox row of the existing pipeline is dropped.

## Technical Context

**Language/Version**: Python 3.10+ (matches CHILmesh and ADMESH minimums)
**Primary Dependencies**: numpy, scipy.spatial.Delaunay (re-triangulation inside truss loop), matplotlib (rendering only), chilmesh (this repo), admesh (sibling repo, used at pinned commit `05bc68fc81060f7d710b8f4abb2cc382f85df33f` per FR-015 / Q1=a — only `admesh.distmesh.distmesh2d` is actually called; `admesh.routine` is bypassed because of the broken `MeshOutput` import)
**Storage**: None (function returns in-memory arrays / CHILmesh instances). Demo PNG written to `tests/output/annulus_quickstart.png` — same path as before so README image links don't break.
**Testing**: pytest. Adapter unit tests exercise both API forms × multiple input sources × multiple domains (per SC-007/SC-008). Demo script self-tests via fail-loud assertions V_BND, V_BND_PROP, V_QI, V_CONN, V_CHAIN before saving the PNG.
**Target Platform**: Cross-platform Python (Linux/macOS/Windows) with numpy ≥ 1.20, scipy ≥ 1.7. No platform-specific code.
**Project Type**: Library + standalone visualization script. New module added to `src/chilmesh/`; demo script lives at repo root.
**Performance Goals**: Adapter completes in < 30 s wall time on a 580-element annulus mesh on a developer laptop (per SC-003). Re-triangulation cadence dominates; default `niter=500` from ADMESH is the cap.
**Constraints**: Bit-exact boundary preservation is the contract — verified by `np.array_equal` on the boundary subset, not approximate equality (per FR-008, SC-002). Must not modify input mesh (per FR-010). Must use deterministic RNG (`seed=0` default) so demo PNG is reproducible.
**Scale/Scope**: Adapter handles meshes from ~100 to ~10,000 elements in V1. Larger meshes (Block_O at 50k) are out of scope for V1 — they may take minutes due to per-iteration Delaunay rebuild. SC-003 wall-time benchmark is for the annulus fixture only.

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

CHILmesh's `.specify/memory/constitution.md` is an unfilled template — no concrete principles ratified. The binding governance document is `CLAUDE.md`, which mandates the single-branch policy. Per the precedence rule (CLAUDE.md > system reminders), this plan is committed on `main` because:

1. The README + image work (spec 004 lineage) was pushed to `main` directly when `planning-optimize_modernize` had divergence issues earlier this session.
2. The remote `planning-optimize_modernize` was force-updated to a different lineage (issue #4 FEM-quad work) and no longer contains spec 004 history.
3. Putting spec 005 on `main` keeps it adjacent to the other in-tree specs that still live there (001, 002, 002-fem-smoother-quads, and now 005).

**Status**: PASS (no constitution gates; branch policy navigated as documented).

## Project Structure

### Documentation (this feature)

```text
specs/005-admesh-warm-start-truss/
├── spec.md                      # /speckit-specify output (clarified)
├── plan.md                      # This file (/speckit-plan output)
├── research.md                  # Phase 0: distmesh2d API, MeshOutput status, pfix preservation
├── data-model.md                # Phase 1: entities, runtime data flow
├── quickstart.md                # Phase 1: how to call the adapter, three worked examples
├── contracts/
│   ├── api-contract.md          # Phase 1: function signatures, validation, errors
│   └── visualization-output.md  # Phase 1: new 4-row PNG structure (replaces spec 004's contract)
└── tasks.md                     # Phase 2 output (/speckit-tasks — NOT created here)
```

### Source Code (repository root)

```text
CHILmesh/
├── src/
│   └── chilmesh/
│       ├── admesh_warmstart.py       # NEW: the adapter module (FR-001a + FR-001b)
│       ├── _vendor_admesh_truss.py   # NEW: vendored inner truss loop (private; remove when ADMESH-B lands)
│       ├── CHILmesh.py               # UNCHANGED
│       ├── examples.py               # UNCHANGED (annulus, donut fixtures)
│       └── __init__.py               # MODIFIED: re-export `optimize_with_admesh_truss` and array form
├── tests/
│   ├── test_admesh_warmstart.py      # NEW: unit tests (FR-001 through FR-018, SC-001 through SC-008)
│   └── output/
│       └── annulus_quickstart.png    # MODIFIED: regenerated with new 4-row layout (Row 1 raw, Row 2 warm-start, Row 3 FEM-of-row2, Row 4 right-iso-of-row2)
├── generate_4row_admesh.py           # MODIFIED: restructured per Q2=d
├── README.md                         # MODIFIED: caption updated to describe new pipeline rows
└── specs/005-admesh-warm-start-truss/  # NEW: this spec/plan/contracts directory
```

**Structure Decision**: Single library module. `optimize_with_admesh_truss` and `optimize_with_admesh_truss_arrays` are added to `src/chilmesh/admesh_warmstart.py`, re-exported from `chilmesh/__init__.py` so callers can write `chilmesh.optimize_with_admesh_truss(...)`. The vendored truss-loop helper lives in `_vendor_admesh_truss.py` (underscore-prefixed = private) and is removed once ADMESH-B is fixed upstream.

## Architecture Overview

```text
┌─────────────────────────────────────────────────────────────────────────┐
│                        Caller (user / demo script)                       │
└──────────────────────────────────┬──────────────────────────────────────┘
                                   │
              ┌────────────────────┴────────────────────┐
              │                                          │
              ▼ (CHILmesh form, FR-001b)                 ▼ (raw arrays form, FR-001a)
   ┌──────────────────────────┐               ┌─────────────────────────────────┐
   │ optimize_with_admesh_    │               │ optimize_with_admesh_truss_     │
   │   truss(mesh, sdf,       │  thin wrapper │   arrays(points, triangles,     │
   │   size_fn, **kwargs) ────┼──────────────►│   sdf, size_fn,                 │
   │   -> CHILmesh            │               │   boundary_indices=None,        │
   └──────────────────────────┘               │   **kwargs)                     │
                                              │   -> (points, triangles)        │
                                              └────────────┬────────────────────┘
                                                           │
                                          ┌────────────────┴────────────────┐
                                          │  Validation layer (FR-005-7,16) │
                                          │  - boundary on SDF zero set     │
                                          │  - triangle-only                │
                                          │  - positive areas               │
                                          └────────────────┬────────────────┘
                                                           │
                                          ┌────────────────┴────────────────┐
                                          │  _vendor_admesh_truss.          │
                                          │  distmesh2d_warmstart(...)      │
                                          │  - byte-identical truss loop    │
                                          │  - skips _initial_distribution  │
                                          │  - skips _rejection_method      │
                                          │  - uses caller's initial points │
                                          └────────────────┬────────────────┘
                                                           │
                                          ┌────────────────┴────────────────┐
                                          │  Non-degradation guard (FR-011) │
                                          │  if median_q(out) < median_q(in)│
                                          │     return input + RuntimeWarning│
                                          └────────────────┬────────────────┘
                                                           │
                                                           ▼
                                                  optimized (points, triangles)
```

Key architectural decisions:

- **Two-tier API**: The high-level form is purely a wrapper that extracts `(points, triangles, boundary_indices)` and rewraps. All real logic — validation, truss invocation, non-degradation guard — lives in the array form. This satisfies FR-001 + FR-016 (input-source-agnostic).
- **Vendored truss loop**: The inner loop is copied byte-for-byte from `admesh.distmesh.distmesh2d` at the pinned SHA. The only changes are (1) the initial point distribution comes from the caller, (2) `pfix` is the caller's boundary subset (instead of an optional argument). The vendored module includes a header comment with the source SHA and a TODO referencing ADMESH-B for removal.
- **Non-degradation as a wrapper layer**: FR-011 (return input on regression) is a thin guard around the truss call, not part of the truss itself. This keeps the vendor module pure and lets us evolve the regression policy without touching the truss math.
- **Validation upfront**: All input validation runs *before* the truss loop. By the time we enter the truss, the input is guaranteed valid; the truss never has to handle malformed inputs.

## Constraints & Assumptions

| Item | Decision | Source |
|------|----------|--------|
| ADMESH version pin | `05bc68fc81060f7d710b8f4abb2cc382f85df33f` (current `main` HEAD as of 2026-05-02) | FR-015, Q1=a |
| MeshOutput workaround | Bypass `admesh.routine` entirely; call `admesh.distmesh.distmesh2d` directly OR use vendored truss loop | research.md R1 |
| Boundary identification (CHILmesh form) | `boundary_edges()` → `Edge2Vert` lookup → `np.unique` on flattened indices | FR-002 |
| Boundary identification (array form) | Caller passes `boundary_indices` OR adapter infers via "edges in exactly one triangle" | FR-002 |
| pfix preservation contract | ADMESH's `Ftot[:nfix] = 0.0` makes pfix bit-exact preserved (verified by reading source); spec depends on this | research.md R2 |
| Output PNG path | `tests/output/annulus_quickstart.png` (unchanged — README links keep working) | FR-012 |
| RNG seed default | `seed=0` (matches ADMESH default; ensures deterministic demo PNG) | FR-009 |
| Non-degradation policy | Return input + `RuntimeWarning` (not raise, not silent worse output) | FR-011, Q4=b |
| Demo benchmark | annulus median quality ≥ 0.60 | SC-001, Q5=b |
| Second test domain | donut (CHILmesh fixture) with caller-supplied SDF | SC-008, Q3=d |

## Phase Plan

The plan executes in three phases. **Phase 0 (Research)** is documented in `research.md` and is already complete — findings inline below. **Phase 1 (Design)** produces `data-model.md`, `quickstart.md`, and the two contracts. **Phase 2 (Tasks)** is generated by `/speckit-tasks` — NOT in this file.

### Phase 0: Research (DONE — see research.md)

Key findings (full detail in `research.md`):

- **R1**: `admesh.routine` imports `MeshOutput` from `admesh.distmesh` but `MeshOutput` is not defined there on `main` HEAD (`05bc68f`). The fix exists on `daily-issue-fixing` branch (`f21e340`) but has not landed on `main`. We bypass `routine` entirely; `admesh.distmesh.distmesh2d` is directly importable and complete.
- **R2**: `distmesh2d`'s `pfix` mechanism stores fixed points at indices `0..nfix-1` and applies `Ftot[:nfix] = 0.0` every iteration — bit-exact preservation is guaranteed by source inspection. Re-triangulation via scipy `Delaunay` does not reorder these indices.
- **R3**: `distmesh2d`'s public signature accepts `(fd, fh, h0, bbox, pfix, *, dptol, ttol, Fscale, deltat, geps_factor, niter, seed)`. All truss-relevant tunables are exposed; we forward them through our adapter as kwargs.
- **R4**: There is **no** way to inject custom initial points into `distmesh2d` without source modification — it always calls `_initial_distribution` then `_rejection_method`. The vendor approach is therefore necessary; there is no monkey-patch hook.
- **R5**: The "right-isoceles smoother" (`admesh.quad_prep.smooth_for_quadrangulation`) is in a different module that does NOT depend on `routine.py` — it imports cleanly from the pinned commit. No change needed for Row 4.

### Phase 1: Design & Contracts

Outputs:

1. **`data-model.md`** — Runtime entities: `WarmStartInput`, `BoundaryVertexSet`, `InteriorVertexSet`, `TrussSolverConfig`, `WarmStartOutput`, plus the four-row figure entities. (No persistent data model.)
2. **`quickstart.md`** — Three worked examples per FR-018: (a) bundled annulus, (b) bundled donut, (c) raw-arrays from a non-CHILmesh source.
3. **`contracts/api-contract.md`** — The function signatures (FR-001a/b), validation rules (FR-005, 006, 007), error catalog, kwargs forwarded to `distmesh2d`, return value shape.
4. **`contracts/visualization-output.md`** — The new 4-row PNG: row layout, column layout (mesh / layers / quality), colormaps (parula / cool_r — preserved from spec 004), fail-loud assertions (V_BND, V_BND_PROP, V_QI, V_CONN, V_CHAIN), title strings.

### Phase 2: Tasks (NOT in this file)

Generated by `/speckit-tasks` from this plan. Will produce a sequenced task list grouping work by user story (US1 boundary preservation = MVP; US2 4-row demo; US3 generic input). Each task names exact file paths and acceptance criteria.

## Dependencies & Risks

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| ADMESH `main` ships further breaking changes | medium | medium | Pin to `05bc68f` exactly; refuse to upgrade without a manual review |
| `pfix` bit-exact preservation isn't actually byte-stable across scipy `Delaunay` versions | low | high | Add a regression test that runs the adapter and compares `output_points[boundary_indices]` to `input_points[boundary_indices]` via `np.array_equal` — if scipy ever introduces drift, this fails loud |
| Vendored truss loop drifts from upstream | medium (if upstream evolves) | low (we're pinned anyway) | Header comment in `_vendor_admesh_truss.py` lists the source SHA; an explicit task in `/speckit-tasks` adds a "diff against upstream pinned commit" check to CI |
| Annulus quality regresses on warm-start (rare but possible if size_fn is hostile) | low | low | Demo uses a known-good constant size_fn; non-degradation guard (FR-011) catches it for arbitrary callers |
| New 4-row layout breaks existing README image link | low | high | PNG path stays at `tests/output/annulus_quickstart.png`; only the *contents* change |
| Donut SDF supplied by tests is non-trivial to write | low | low | The donut is two concentric annuli; SDF is `max(max(r-R_outer, R_inner-r), max(r-R_inner_hole, R_inner_hole_outer-r))` or similar. Write it once in test fixtures and reuse |
| Spec 004's V1-V7 assertions in `generate_4row_admesh.py` no longer match the new pipeline | high | medium | Plan task: REPLACE V1-V7 with V_BND, V_BND_PROP, V_QI, V_CONN, V_CHAIN; comment in script explains the migration |

## Open Questions Pushed to /speckit-tasks or Implementation

These are decisions the planner consciously defers:

- **Concrete demo size function for warm-start**: the user said "constant size function" implicitly via the spec, but the original `generate_4row_admesh.py` used a both-boundary refinement size function. Should warm-start (Row 2) use constant or graded? Default plan: constant for Row 2 (so the demo is purely about quality improvement, not size grading); the graded function from spec 004 is removed because that row is dropped.
- **Whether to keep `generate_4row_admesh.py` filename or rename to `generate_4row_warmstart.py`**: Plan says keep the filename to avoid invalidating any external bookmarks; the script's *narrative* changes but the filename and PNG path do not. Implementation can revisit if the team prefers a clean rename.
- **Whether to add a new fixture `examples.warm_start_demo()`**: out of scope for V1.

## Complexity Tracking

The constitution check is PASS (no gates), so no complexity violations to justify. The architectural choice that warrants explicit tracking:

- **Vendored truss loop** is a deliberate temporary measure pending ADMESH-B. It's complexity-positive in the short term (we maintain a copy) but complexity-negative long term once ADMESH-B lands and we delete the file. The vendor module is ~50 lines (the truss inner loop plus the warm-start preamble); the cost is bounded.

## Acceptance for `/speckit-plan`

This plan is ready when:

- [x] Summary captures the technical approach in one paragraph
- [x] Technical Context filled with concrete versions / dependencies / pins
- [x] Constitution Check completed (PASS)
- [x] Project Structure shows files added / modified
- [x] Architecture diagram shows the call graph and key decisions
- [x] Constraints & Assumptions table aligns with clarify answers Q1-Q5
- [x] Phase 0 research findings inline (full detail in research.md)
- [x] Phase 1 deliverables enumerated
- [x] Risk register populated with concrete mitigations
- [x] Open questions explicitly pushed to /speckit-tasks or implementation

Status: **Ready for `/speckit-tasks`**.
