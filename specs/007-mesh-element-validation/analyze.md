# Analyze: spec ↔ plan ↔ tasks consistency

**Date**: 2026-05-21
**Spec**: `spec.md` (17 FRs, 8 SCs, 4 user stories, 4 clarifications Q1-Q9)
**Plan**: `plan.md`
**Tasks**: `tasks.md` (T001-T063)

## Coverage matrix

### Functional Requirements → Plan section → Task

| FR | Plan section | Task(s) | Status |
|----|--------------|---------|--------|
| FR-001 | Module layout / Surface | T001, T003, T021 | ✅ |
| FR-002 | Data model | T003 | ✅ |
| FR-003 | Phase 4 US1 | T021 | ✅ |
| FR-004 | Phase 4 US1 | T021 | ✅ |
| FR-005 | Phase 4 US1 | T021 | ✅ |
| FR-006 | Phase 4 US1 | T021 | ✅ |
| FR-007 | Phase 4 US1 + risk register | T021 | ✅ |
| FR-008 | Phase 5 US2 | T031, T032 | ✅ |
| FR-009 | Phase 5 US2 + research bowtie | T032 | ✅ |
| FR-010 | Phase 5 US2 (no-op note) | T032 | ✅ trivially |
| FR-011 | Phase 5 US2 | T033 | ✅ |
| FR-012 | Phase 5 US2 + research broadphase | T013, T033, T053 | ✅ |
| FR-013 | Tech context + Constitution V | T012 | ✅ |
| FR-014 | Phase 6 US3 | T041 | ✅ |
| FR-015 | Phase 1 design + research fixture construction | T014 | ✅ |
| FR-016 | Phase 7 US4 | T052 | ✅ |
| FR-017 | Phase 7 US4 | T052 | ✅ |

### Success Criteria → Task

| SC | Task | Status |
|----|------|--------|
| SC-001 | T022, T034, T042 | ✅ |
| SC-002 | T022, T034, T040 | ✅ |
| SC-003 | T001, T003 | ✅ |
| SC-004 | T051 | ✅ |
| SC-005 | T062 | ✅ |
| SC-006 | T050 | ✅ |
| SC-007 | research.md + T053 | ✅ (verified empirically post-impl) |
| SC-008 | done — issue #142 already open | ✅ |

### User stories → Tasks

| US | Priority | Tasks | Independent test |
|----|----------|-------|------------------|
| US1 element-type | P1 | T020-T022 | T020 parametrization |
| US2 geometric | P1 | T030-T034 | T030 parametrization |
| US3 degenerate-accepted | P1 | T040-T042 | T040 explicit tests |
| US4 pytest UX | P2 | T050-T053 | T050 + T051 |

### Edge cases → Coverage

| Edge case | Where covered |
|-----------|---------------|
| Triangle-as-padded-quad | T021 (classify_element padded form), T040 |
| Coincident vertices (distinct IDs) | T040 + T032 (bowtie predicate handles via robust orient) |
| Numerical near-bowtie | Predicate tol semantics (research.md) + T010 |
| Floating-point coincident edges | Edge-sharing map uses vert-ID equality (research.md); coincident-coord-but-distinct-ID treated as not shared (acceptable — would be `INTERIOR_OVERLAP` if it matters) |
| Boundary tri (all on boundary) | T021 FR-006 |
| Boundary tri (1 on boundary) | T021 FR-006 |
| Interior tri (0 on boundary) | T021 → INTERIOR_TRIANGLE_FORBIDDEN |
| Disconnected components | Validator iterates all elements; no component-level branching needed |
| Open-ocean boundaries | Treated as boundary by `mesh.boundary_vertices()` semantics |
| Empty mesh | T021 short-circuit + EMPTY_MESH note |
| Single-element mesh | T033 broadphase yields zero pairs; per-element rules still run |

## Constitution gates re-checked post-plan

| Principle | Status |
|-----------|--------|
| I library-first | ✅ self-contained `tests/_validity/` |
| II immutability | ✅ no mesh mutation; `_skeletonize()` cache write is existing behavior |
| III test-first | ✅ every implementation task preceded by a test task |
| IV correctness > perf | ✅ bowtie via segment-cross, not area sign |
| V coord-agnostic | ✅ bbox-relative tol |
| VI format pluralism | ✅ N/A |
| VII API stability | ✅ no public API change |

## Risks flagged in plan, status

| Risk | Mitigation task |
|------|----------------|
| block_o runtime > 45 s | T053 fallback to KDTree |
| Constructor degeneracy fallback repairs synthetic bowtie | T014 corrupt_to + cache-bust |
| `_skeletonize()` adds ~14 s on block_o (FR-007) | Cached; one-time. Counted toward T051 budget |
| `boundary_vertices()` shape mismatch | T021 wraps in `set(int(x) for x in ...)` |
| FP ties in segment-cross | T010 covers near-tol cases |

## Inconsistencies / gaps found and resolved

1. **`pentagon_mesh` cannot use real `CHILmesh`** (4-column connectivity_list rules out 5-vertex elements). Resolution: data-model.md introduces a `_FakeMesh` test double for this single negative fixture. Plan and tasks updated to call it out (T014, fixtures.py).
2. **`LAYERS_AUTO_TRIGGERED` note category was implied by FR-007 but not in the spec's violation-category table.** Resolution: T060 polish task adds it.
3. **`_FakeMesh` necessity introduces a tiny duck-type contract**: `connectivity_list`, `points`, `n_elems`, `boundary_vertices`, `layers`, `boundary_edges`. The validator MUST NOT call any other attribute on the mesh, or the synthetic fixture breaks. Plan section "Module layout" implicitly enforces this; flagged here for T021 / T033 implementer.
4. **FR-013 tolerance flow** is documented in research.md and plan.md but the spec text could be tighter on the floor value `1e-15`. Acceptable as-is.

## Approval to proceed to implement

All FRs and SCs trace to tasks. No constitution violations. No unaddressed risks. Spec is internally consistent.

**Ready for Phase 3 (`/speckit.implement`).**
