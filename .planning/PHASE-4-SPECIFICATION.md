# Phase 4: Downstream Integration & Release Specification

**Branch:** `planning-optimize_modernize`  
**Date:** 2026-04-26  
**Phase:** 4 – Integration & Release (0.2.0)  
**Methodology:** Spec-Kit  
**Prerequisite:** Phases 1–3 complete

---

## Executive Summary

Phase 4 integrates CHILmesh's modernized API with MADMESHR, ADMESH-Domains, hypothetical ADMESH. Culminates in v0.2.0 with breaking-change signaling, migration guide, architectural decisions documented.

**Deliverables:**
- MADMESHR integration (advancing-front API, domain splitting)
- ADMESH-Domains bulk-load verification
- v0.2.0 release with migration guide
- Governance docs updated (CLAUDE.md, constitution.md, PROJECT_PLAN.md, API.md)
- Lessons learned doc

**Success:** MADMESHR uses CHILmesh for advancing-front generation; 0.2.0 released with clear signaling.

---

## Specification

### 4.1 Goals

**Primary:** Ship v0.2.0 as breaking-change release with clear migration path.

**Secondary:** Validate MADMESHR and ADMESH-Domains leverage new CHILmesh APIs.

### 4.2 In-Scope Integration

#### **4.2.1 MADMESHR Advancing-Front API**

**What MADMESHR needs:**
1. Place quads/triangles incrementally on advancing-front boundary
2. Detect pinch points (bottlenecks); split domain
3. Remove boundary loop elements; pack into result

**CHILmesh Support:**
```python
class CHILmesh:
    def advancing_front_boundary_edges(self) -> List[int]:
        """Return edge IDs on current mesh boundary."""
        return self.adjacencies['Edge2Elem'].apply(
            lambda e: e[1] == -1  # Boundary edge has -1 in right elem
        )
    
    def add_advancing_front_element(
        self, vertices: List[int], elem_type: str
    ) -> int:
        """Add element to mesh during advancing-front generation."""
        # Wrapper around add_element() with advancing-front semantics
        return self.add_element(vertices, elem_type)
    
    def remove_boundary_loop(self, edge_ids: List[int]) -> None:
        """Remove boundary loop (residual closure) when shrinks to ≤4 verts."""
        for edge_id in edge_ids:
            # Remove adjacent elements
            pass
    
    def pinch_points(self, width_threshold: float) -> List[int]:
        """From Phase 3: identify bottlenecks."""
        return self._pinch_points(width_threshold)
```

**Integration Path:**
- MADMESHR calls `mesh.add_advancing_front_element(verts, type)` in loop
- CHILmesh state stays consistent (adjacencies updated, layers stale until rebuild)
- MADMESHR queries `advancing_front_boundary_edges()` to find next placement

**Test Case:** Simple advancing-front scenario (annulus → fill center) using MADMESHR API.

#### **4.2.2 ADMESH-Domains Bulk Load Optimization**

**What ADMESH-Domains needs:**
- Fast mesh loading (O(n log n) vs. O(n²))
- Metadata queries (size, quality, type)
- Lazy init (skip skeletonization unless requested)

**CHILmesh Support (Already From Phase 1–3):**
```python
class CHILmesh:
    def __init__(self, ..., compute_layers: bool = False):
        """If compute_layers=False, skip skeletonization (saves time for ADMESH-Domains)."""
        self._initialize_mesh(compute_layers=compute_layers)
    
    # Phase 3 optimizations make _identify_edges() fast
    # ADMESH-Domains can load, query metadata, export without layer computation
```

**Verification:** Benchmark ADMESH-Domains loading 10 sample meshes; confirm <500ms each.

#### **4.2.3 ADMESH (Hypothetical) Readiness**

**Status:** GitHub 404 — future availability assumed. Prepare APIs without direct integration.

**Likely needs (from thesis):**
- Edge-swapping (improve quality)
- Node repositioning (smoothing)
- Refinement/coarsening (local modification)

**CHILmesh Readiness (Phase 2):**
- `add_element()`, `remove_element()` → refinement/coarsening
- `swap_edge()` → edge-swapping
- Existing smoothing methods (FEM, angle-based) → node repositioning

**No additional work needed; Phase 2–3 APIs sufficient.**

### 4.3 Documentation Updates

1. **API.md** (New) — public methods: deprecated vs. new (0.2.0+); examples; deprecation warnings
2. **MIGRATION_GUIDE.md** (New) — 0.1.x → 0.2.0 breaking changes; old vs. new code examples; FAQ
3. **CLAUDE.md** (Update) — add "Modernization Task" section; link MIGRATION_GUIDE.md
4. **constitution.md** (Update) — add "Graph Representation Governance"; document Decisions A1–A5
5. **PROJECT_PLAN.md** (Update) — add "0.2.0 Modernization Release" section; roadmap updated
6. **CHANGELOG.md** (Update) — "0.2.0" entry with breaking changes; list new methods; ref MIGRATION_GUIDE.md

### 4.4 Governance Documents

Create `MODERNIZATION_LESSONS_LEARNED.md`:
- What worked (spec-kit, benchmarks, research docs)
- What didn't
- Recommendations for future modernization
- Architectural decisions + rationale

### 4.5 Test Coverage

**Integration Tests:**
- [ ] MADMESHR advancing-front scenario (add 10 elements, verify valid)
- [ ] ADMESH-Domains bulk load (5 meshes, <500ms each)
- [ ] Fort.14 roundtrip after dynamic ops (modify, export, reimport, identical)
- [ ] All Phase 1–3 tests passing

**Release Tests:**
- [ ] Version bumped to 0.2.0 in setup.py, pyproject.toml
- [ ] CHANGELOG.md updated with breaking-change section
- [ ] MIGRATION_GUIDE.md exists and complete
- [ ] API.md documents all public methods

### 4.6 Acceptance Criteria

- [ ] MADMESHR integration API defined and tested
- [ ] ADMESH-Domains bulk-load benchmarked (<500ms/mesh)
- [ ] Version bumped to 0.2.0
- [ ] CHANGELOG.md, MIGRATION_GUIDE.md, API.md created
- [ ] CLAUDE.md, constitution.md, PROJECT_PLAN.md updated
- [ ] Lessons learned doc written
- [ ] All tests passing (0 regressions)
- [ ] RC1 tagged and tested

### 4.7 Out of Scope

- ❌ Publishing to PyPI (separate release process)
- ❌ Porting MATLAB original to Python (done in 0.1.x)
- ❌ New features beyond downstream integration
- ❌ CI/CD setup (separate audit follow-up)

---

## Timeline

**2–3 days** (coordination + docs; APIs from Phase 3 ready)

---

## Final Gate (Release)

Before 0.2.0 ship:

1. ✓ Issues #35–38 (Phase 1 research) closed
2. ✓ Phase 2 complete + all tests passing
3. ✓ Phase 3 complete + benchmarks pass
4. ✓ MADMESHR integration validated
5. ✓ ADMESH-Domains bulk-load verified
6. ✓ Docs complete (API.md, MIGRATION_GUIDE.md, etc.)
7. ✓ Version bumped to 0.2.0
8. ✓ Lessons learned approved

Then: Tag RC1, get final approval, publish to PyPI.

---

## Handoff to Future Maintenance

After 0.2.0:
- Break into Issues #40+ for maintenance/feature work
- Standing task: monitor MADMESHR integration, support evolving needs
- Quarterly review: ADMESH updates? MADMESHR feedback?

