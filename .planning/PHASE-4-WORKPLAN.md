# Phase 4: Implementation Workplan

**Status:** In Progress  
**Branch:** `planning-optimize_modernize`  
**Target:** 0.2.0 release with full integration

---

## Task Breakdown

### 1. MADMESHR Advancing-Front API
- [ ] **P4-01:** Implement `advancing_front_boundary_edges()` — returns boundary edge IDs; queries Edge2Elem (e2 == -1)
- [ ] **P4-02:** Implement `add_advancing_front_element(vertices, elem_type)` — incremental add; update adjacencies in-place; return new elem ID
- [ ] **P4-03:** Implement `remove_boundary_loop(edge_ids)` — remove residual boundary closure; update adjacencies/layers
- [ ] **P4-04:** Test advancing-front scenario (annulus fill center) — 10+ elements via API; verify validity; test Fort.14 export/reimport

### 2. Documentation

- [ ] **P4-05:** Create `API.md` — public method reference (25+ methods); old vs. new (0.2.0); deprecation warnings; usage examples
- [ ] **P4-06:** Update `CHANGELOG.md` — "0.2.0" section with breaking changes; list new methods; ref MIGRATION_GUIDE.md
- [ ] **P4-07:** Update governance docs — CLAUDE.md: "Modernization Task"; constitution.md: "Graph Representation Governance"; PROJECT_PLAN.md: "0.2.0 Modernization Release"

### 3. Version & Release

- [ ] **P4-08:** Bump version to 0.2.0 — setup.py, pyproject.toml, __init__.py (if needed)
- [ ] **P4-09:** Create lessons learned doc — what worked (spec-kit, benchmarks, research); what to improve; architectural decisions + rationale

### 4. Verification

- [ ] **P4-10:** ADMESH-Domains bulk-load benchmark — 5+ meshes; <500ms each (compute_layers=False); document results
- [ ] **P4-11:** All tests passing — 239 existing + 4+ advancing-front; zero regressions
- [ ] **P4-12:** Tag release candidate — v0.2.0-rc1; push to remote; verify stable

---

## Acceptance Criteria

- [x] Phases 1-3 complete (prerequisite)
- [ ] 4 MADMESHR API methods implemented
- [ ] 4+ advancing-front integration tests
- [ ] API.md complete (all public methods)
- [ ] CHANGELOG.md with 0.2.0 section
- [ ] Version bumped to 0.2.0
- [ ] CLAUDE.md, constitution.md, PROJECT_PLAN.md updated
- [ ] Lessons learned doc written
- [ ] ADMESH-Domains bulk-load verified (<500ms/mesh)
- [ ] All tests passing (239+)
- [ ] v0.2.0-rc1 tagged and pushed

---

## Dependencies

- All Phase 1-3 issues closed
- WNAT_Hagen benchmark complete (Issue #55)
- DOWNSTREAM_MIGRATION_GUIDE.md exists
- Bridge adapters functional

---

## Timeline

Estimated: 2-3 days (docs + coordination, not heavy coding)
