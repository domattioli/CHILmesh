# Phase 4: Downstream Integration & Release Specification

**Branch:** `planning-optimize_modernize`  
**Date:** 2026-04-26  
**Phase:** 4 – Integration & Release (0.2.0)  
**Methodology:** Spec-Kit  
**Prerequisite:** Phases 1–3 complete

---

## Executive Summary

Phase 4 integrates CHILmesh's modernized API with MADMESHR, ADMESH-Domains, and hypothetical ADMESH. The phase culminates in releasing version 0.2.0 with breaking-change signaling, migration guide, and documented architectural decisions.

**Deliverables:**
- MADMESHR integration (advancing-front API, domain splitting)
- ADMESH-Domains bulk-load optimization verification
- Version 0.2.0 release with migration guide
- Governance docs updated (CLAUDE.md, constitution.md, PROJECT_PLAN.md, API.md)
- Lessons learned document

**Success:** MADMESHR can use CHILmesh for advancing-front generation; 0.2.0 released with clear signaling.

---

## Specification

### 4.1 Goals

**Primary:** Ship version 0.2.0 as breaking-change release with clear migration path.

**Secondary:** Validate that MADMESHR and ADMESH-Domains can leverage new CHILmesh APIs effectively.

### 4.2 In-Scope Integration

#### **4.2.1 MADMESHR Advancing-Front API**

**What MADMESHR needs:**
1. Place quads/triangles incrementally on advancing-front boundary
2. Detect pinch points (bottlenecks) and split domain
3. Remove boundary loop elements and pack into result

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
- CHILmesh internal state stays consistent (adjacencies updated, layers stale until rebuild)
- MADMESHR can query `advancing_front_boundary_edges()` to find where to place next element

**Test Case:** Simple advancing-front scenario (annulus → fill center) using MADMESHR API.

#### **4.2.2 ADMESH-Domains Bulk Load Optimization**

**What ADMESH-Domains needs:**
- Fast mesh loading (O(n log n) vs. O(n²))
- Metadata queries (size, quality, type)
- Lazy initialization (don't compute skeletonization unless requested)

**CHILmesh Support (Already From Phase 1–3):**
```python
class CHILmesh:
    def __init__(self, ..., compute_layers: bool = False):
        """If compute_layers=False, skip skeletonization (saves time for ADMESH-Domains)."""
        self._initialize_mesh(compute_layers=compute_layers)
    
    # Phase 3 optimizations make _identify_edges() fast
    # ADMESH-Domains can load, query metadata, export without layer computation
```

**Verification:** Benchmark ADMESH-Domains loading 10 sample meshes; confirm <500ms per mesh.

#### **4.2.3 ADMESH (Hypothetical) Readiness**

**Status:** GitHub 404 — assuming future availability. Prepare APIs without direct integration.

**Likely needs (from thesis):**
- Edge-swapping (improve mesh quality)
- Node repositioning (smoothing)
- Refinement/coarsening (local modification)

**CHILmesh Readiness (From Phase 2):**
- `add_element()`, `remove_element()` → supports refinement/coarsening
- `swap_edge()` (Phase 2 or 3) → edge-swapping support
- Existing smoothing methods (FEM, angle-based) → node repositioning

**No additional work needed for Phase 4; APIs from Phases 2–3 are sufficient.**

### 4.3 Documentation Updates

**Create/Update:**

1. **API.md** (New)
   - Public methods: old (deprecated) vs. new (0.2.0+)
   - Examples: adding elements, querying adjacencies, detecting components
   - Deprecation warnings (which old methods to avoid)

2. **MIGRATION_GUIDE.md** (New)
   - 0.1.x → 0.2.0 breaking changes
   - Code examples: old way vs. new way
   - FAQ: "Why did you change the adjacency dict?"

3. **CLAUDE.md** (Update)
   - Add "Modernization Task" section (copy from GOVERNANCE_UPDATES.md Part 3a)
   - Link to MIGRATION_GUIDE.md

4. **constitution.md** (Update)
   - Add "Graph Representation & Data Structure Governance" section (from GOVERNANCE_UPDATES.md Part 3b)
   - Document Decisions A1–A5

5. **PROJECT_PLAN.md** (Update)
   - Add "0.2.0 Modernization Release" section (from GOVERNANCE_UPDATES.md Part 3c)
   - Roadmap updated, Phase 4 complete

6. **CHANGELOG.md** (Update)
   - New "0.2.0" entry with breaking-change section
   - List all new methods
   - Reference MIGRATION_GUIDE.md

### 4.4 Governance Documents

**Commit Lessons Learned Document** (final deliverable):

Create `MODERNIZATION_LESSONS_LEARNED.md` summarizing:
- What worked (spec-kit methodology, benchmark framework, research docs)
- What didn't (if anything)
- Recommendations for future modernization efforts
- Architectural decisions made + rationale

See Part 5 of this spec.

### 4.5 Test Coverage

**Integration Tests:**
- [ ] MADMESHR advancing-front scenario (add 10 elements, verify mesh valid)
- [ ] ADMESH-Domains bulk load (load 5 meshes, time <500ms each)
- [ ] Fort.14 roundtrip after dynamic ops (modify, export, reimport, identical)
- [ ] All Phase 1–3 tests passing

**Release Tests:**
- [ ] Version bumped to 0.2.0 in setup.py, pyproject.toml
- [ ] CHANGELOG.md updated with breaking-change section
- [ ] MIGRATION_GUIDE.md exists and is complete
- [ ] API.md documents all public methods

### 4.6 Acceptance Criteria

- [ ] MADMESHR integration API defined and tested
- [ ] ADMESH-Domains bulk-load benchmarked (<500ms per mesh)
- [ ] Version bumped to 0.2.0
- [ ] CHANGELOG.md, MIGRATION_GUIDE.md, API.md created
- [ ] CLAUDE.md, constitution.md, PROJECT_PLAN.md updated
- [ ] Lessons learned document written
- [ ] All tests passing (0 regressions)
- [ ] Release candidate (RC1) tagged and tested

### 4.7 Out of Scope (Phase 4)

- ❌ Publishing to PyPI (separate release process)
- ❌ Porting MATLAB original to Python (already done in 0.1.x)
- ❌ New feature requests beyond downstream integration
- ❌ CI/CD setup (separate audit follow-up)

---

## Timeline

**2–3 days** (most work is coordination + documentation; APIs from Phase 3 are ready)

---

## Final Gate (Release)

Before 0.2.0 ship:

1. ✓ All issues #35–38 (Phase 1 research) closed
2. ✓ Phase 2 implementation complete + all tests passing
3. ✓ Phase 3 optimization complete + benchmarks meet targets
4. ✓ MADMESHR integration validated
5. ✓ ADMESH-Domains bulk-load optimization verified
6. ✓ Documentation complete (API.md, MIGRATION_GUIDE.md, etc.)
7. ✓ Version bumped to 0.2.0
8. ✓ Lessons learned document approved

Then: Tag release candidate (RC1), get final approval, publish to PyPI.

---

## Handoff to Future Maintenance

After 0.2.0 release:
- Break into Issues #40+ for maintenance/feature work
- Standing task: "Monitor MADMESHR integration, support evolving needs"
- Quarterly review: Any ADMESH updates? MADMESHR feedback?

