# Phase 4: Implementation Workplan

**Status:** In Progress  
**Branch:** `planning-optimize_modernize`  
**Target:** 0.2.0 release with full integration

---

## Task Breakdown

### 1. MADMESHR Advancing-Front API
- [ ] **P4-01:** Implement `advancing_front_boundary_edges()` method
  - Returns edge IDs on current mesh boundary
  - Queries Edge2Elem for boundary (e2 == -1)
  
- [ ] **P4-02:** Implement `add_advancing_front_element(vertices, elem_type)` method
  - Add element to mesh incrementally
  - Update adjacencies in-place
  - Return new element ID
  
- [ ] **P4-03:** Implement `remove_boundary_loop(edge_ids)` method
  - Remove residual boundary closure
  - Update adjacencies/layers

- [ ] **P4-04:** Test advancing-front scenario (annulus fill center)
  - Add 10+ elements via advancing-front API
  - Verify mesh validity and consistency
  - Test Fort.14 export/reimport

### 2. Documentation

- [ ] **P4-05:** Create `API.md`
  - Public method reference (all 25+ methods)
  - Old vs. new (0.2.0) comparison
  - Deprecation warnings
  - Usage examples
  
- [ ] **P4-06:** Update `CHANGELOG.md`
  - New "0.2.0" section with breaking changes
  - List all new methods
  - Reference MIGRATION_GUIDE.md
  
- [ ] **P4-07:** Update governance docs
  - CLAUDE.md: Add "Modernization Task" section
  - constitution.md: Add "Graph Representation Governance"
  - PROJECT_PLAN.md: Add "0.2.0 Modernization Release"

### 3. Version & Release

- [ ] **P4-08:** Bump version to 0.2.0
  - setup.py
  - pyproject.toml
  - __init__.py (if needed)
  
- [ ] **P4-09:** Create lessons learned document
  - What worked (spec-kit, benchmarks, research docs)
  - What could improve
  - Architectural decisions + rationale

### 4. Verification

- [ ] **P4-10:** ADMESH-Domains bulk-load benchmark
  - Load 5+ sample meshes from catalog
  - Verify each <500ms (with compute_layers=False)
  - Document results
  
- [ ] **P4-11:** All tests passing
  - 239 existing tests
  - 4+ new advancing-front integration tests
  - Zero regressions

- [ ] **P4-12:** Tag release candidate
  - Tag: v0.2.0-rc1
  - Push to remote
  - Verify tag is stable

---

## Acceptance Criteria

- [x] Phases 1-3 complete (prerequisite)
- [ ] 4 MADMESHR API methods implemented
- [ ] 4+ advancing-front integration tests
- [ ] API.md complete with all public methods
- [ ] CHANGELOG.md with 0.2.0 section
- [ ] Version bumped to 0.2.0
- [ ] CLAUDE.md, constitution.md, PROJECT_PLAN.md updated
- [ ] Lessons learned document written
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

Estimated: 2-3 days (documentation + coordination, not heavy coding)
