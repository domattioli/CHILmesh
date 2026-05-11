# CHILmesh Project Plan & Roadmap

**Current Version:** 0.1.1 (Alpha)
**Planning Date:** 2026-04-26
**Horizon:** 12 months (through Q1 2027)

---

## Executive Summary

CHILmesh transitions from bug-fix phase (0.1.1) to data structure modernization + bridge architecture (0.2.0):

1. **Immediate priorities** (Weeks 1-4): GitHub issues, stakeholder coordination
2. **Development phases** (Months 2-6): Three phases of architecture improvements
3. **Stabilization** (Months 7-9): Testing, performance validation, docs
4. **Release & communication** (Months 10-12): v0.2.0 release, downstream integration

**Key Success Criteria:**
- ✅ 1.5×+ performance improvement on large meshes
- ✅ Zero breaking changes to public API
- ✅ Clear bridge interfaces for MADMESHR/ADMESH/ADMESH-Domains
- ✅ 100% test pass rate on all platforms
- ✅ Complete docs of architecture decisions

---

## Phase 0: Planning & Coordination (CURRENT — Weeks 1-4)

**Deliverables:**
- ✅ `PLANNING_DATA_STRUCTURE_MODERNIZATION.md`
- ✅ `CLAUDE.md`
- ✅ `constitution.md`
- ✅ `project_plan.md` (this doc)
- 🔄 GitHub issues (P1-01 through P3-04, 16 total)
- 🔄 Communication with MADMESHR/ADMESH/ADMESH-Domains authors

**Timeline:**
- Week 1-2: Complete planning docs ← **You are here**
- Week 3: Publish to main, create GitHub issues
- Week 4: Stakeholder feedback; adjust timeline if needed

---

## Phase 1: Hash Map Edge Lookup (Weeks 5-8)

**Focus:** Eliminate O(n²) bottleneck in adjacency building  
**Effort:** 8-12 hours | **Risk:** Low

| ID | Task | Hours | Dependencies |
|----|----|-------|--------------|
| P1-01 | Create EdgeMap class | 2 | None |
| P1-02 | Refactor _identify_edges() | 2 | P1-01 |
| P1-03 | Optimize _build_elem2edge() | 2 | P1-02 |
| P1-04 | Optimize _build_edge2elem() | 2 | P1-02 |
| P1-05 | Store EdgeMap in adjacencies | 1 | P1-03, P1-04 |
| P1-06 | Performance regression tests | 3 | P1-05 |

**Success:** all tests pass; Block_O build <45s; ≥1.5× improvement on large fixtures; no new dependencies

**Outputs:** `src/chilmesh/mesh_topology.py`; updated `_build_adjacencies()`; benchmark report

---

## Phase 2: Adjacency Modernization (Weeks 9-14)

**Focus:** Modernize sparse adjacency structures (Vert2Edge, Vert2Elem)  
**Effort:** 12-16 hours | **Risk:** Medium

| ID | Task | Hours | Dependencies |
|----|----|-------|--------------|
| P2-01 | Migrate Vert2Edge to dict | 2 | P1-06 |
| P2-02 | Migrate Vert2Elem to dict | 2 | P1-06 |
| P2-03 | Add explicit accessor methods | 2 | P2-01, P2-02 |
| P2-04 | Update traversal patterns | 4 | P2-03 |
| P2-05 | Type hints & validation | 3 | P2-04 |
| P2-06 | Adjacency documentation | 3 | P2-05 |

**Success:** all tests pass; Vert2Edge/Vert2Elem as `Dict[int, Set[int]]`; type hints complete; no list-of-lists iteration; performance ≥ Phase 1

**Outputs:** refactored `_build_vert2edge()` and `_build_vert2elem()`; new `get_vertex_edges()`, `get_vertex_elements()`; comprehensive adjacency docs

---

## Phase 3: Bridge Infrastructure (Weeks 15-26)

**Focus:** Clean bridge APIs for downstream projects  
**Effort:** 20-24 hours | **Risk:** High (depends on downstream coordination)

| ID | Task | Hours | Dependencies |
|----|----|-------|--------------|
| P3-01 | Define CAI (CHILmesh Access Interface) | 4 | P2-06 |
| P3-02 | Create bridge adapters | 10 | P3-01 |
| P3-03 | Integration tests | 6 | P3-02 |
| P3-04 | Downstream migration guide | 4 | P3-03 |

**Success:** CAI published; bridge adapters for MADMESHR/ADMESH/ADMESH-Domains; integration tests passing; no external dependency on internals

**Outputs:** `docs/CHILmesh_Access_Interface.md`; `src/chilmesh/bridge.py`; `tests/test_bridge_integration.py`; `docs/DOWNSTREAM_MIGRATION_GUIDE.md`

**Note:** requires early feedback from downstream authors; may need timeline adjustment

---

## Phase 4: Stabilization & Testing (Weeks 27-36)

**Focus:** Comprehensive testing, performance validation, final hardening  
**Effort:** 16-20 hours | **Risk:** Low

**Activities:**
1. Full regression testing — Python 3.10/3.11/3.12; Ubuntu/macOS/Windows; fort.14 roundtrip all fixtures
2. Performance benchmarking — baseline before Phase 1, after Phase 1, after Phase 2; document in release notes
3. Memory profiling — verify no >10% increase; no leaks in new structures
4. Documentation review — update README; verify all examples; add architecture docs
5. Community feedback — share RC with downstream; gather feedback; address before final release

**Deliverables:** test coverage report (>90% target); perf benchmark report; memory profiling results; updated README; resolved feedback from downstream

---

## Phase 5: Release & Communication (Weeks 37-40)

**Effort:** 4-6 hours

**Release Checklist:**
- [ ] All tests pass on all platforms
- [ ] Performance benchmarks meet targets
- [ ] CHANGELOG.md updated
- [ ] Version bumped to 0.2.0 in pyproject.toml
- [ ] GitHub release notes prepared
- [ ] PyPI build validated with `twine check`
- [ ] Migration guide ready
- [ ] Announce in CHIL Lab channels
- [ ] Tag v0.2.0 in git
- [ ] Publish to PyPI

**Communication:** GitHub Release (detailed notes); email to MADMESHR/ADMESH/ADMESH-Domains (integration guide); lab announcement

---

## Post-Release (Months 10-12)

**v0.2.1-0.2.x:** bug fixes from community feedback; docs improvements; minor perf enhancements

**v0.3.0 exploration:** gather downstream requirements; evaluate spatial indexing; mesh mutation APIs; profile memory on largest real-world meshes

**v1.0.0 vision:** API freeze assessment; GPU acceleration for smoothing/quality; distributed mesh handling; external format support (Gmsh, VTK)

---

## Resource Allocation

- **Lead Developer:** Claude Code (AI-assisted)
- **Code Review:** Dominik Mattioli (scientific oversight)
- **Testing:** Automated via pytest CI/CD
- **No new external dependencies planned for 0.2.0** (numpy, scipy, matplotlib only)

---

## Risk Management

| Risk | Impact | Mitigation |
|------|--------|-----------|
| Skeletonization breaks | Critical | Comprehensive tests, careful refactoring |
| Downstream incompatibility | Critical | Early coordination, integration tests |
| Performance regression | High | Benchmark every step, revert if slower |
| Complex merge conflicts | High | Clear branch strategy, small focused PRs |
| Scope creep | Medium | Strict phase boundaries, issue templates |
| Timeline slip | Medium | Realistic estimates, track velocity |

---

## Success Metrics

**Quantitative:**
- ✅ 1.5×+ improvement on adjacency building
- ✅ 100% test pass rate on all platforms
- ✅ Zero breaking changes to public API
- ✅ 95%+ of public methods have docstrings
- ✅ 3+ downstream projects successfully integrating

**Qualitative:**
- ✅ Clearer adjacency structures (easier to understand)
- ✅ Fewer complexity-related open issues
- ✅ Positive feedback from downstream authors

---

## Phase Summary Chart

```
Phase 0 (Weeks 1-4): Planning ━━━━━━━━━━━
   ├── Planning docs ✅
   ├── GitHub issues 🔄
   └── Stakeholder coordination 🔄

Phase 1 (Weeks 5-8): Edge Lookup ━━━━━━
   └── 8-12 hours, Low risk

Phase 2 (Weeks 9-14): Adjacencies ━━━━━━━━
   └── 12-16 hours, Medium risk

Phase 3 (Weeks 15-26): Bridge ━━━━━━━━━━━━━━
   └── 20-24 hours, High risk

Phase 4 (Weeks 27-36): Stabilization ━━━━━━━━━
   └── 16-20 hours, Low risk

Phase 5 (Weeks 37-40): Release ━━━━
   └── 4-6 hours, Low risk
```

---

## File Structure After Modernization

```
CHILmesh/
├── src/chilmesh/
│   ├── CHILmesh.py              (refactored: Phase 1-2)
│   ├── mesh_topology.py         (new: Phase 1)
│   ├── bridge.py                (new: Phase 3)
│   ├── utils/
│   │   ├── plot_utils.py        (unchanged)
│   │   └── __init__.py
│   └── data/                    (unchanged)
├── docs/
│   ├── CHILmesh_Access_Interface.md         (new: Phase 3)
│   └── DOWNSTREAM_MIGRATION_GUIDE.md        (new: Phase 3)
├── tests/
│   ├── test_bridge_integration.py           (new: Phase 3)
│   ├── test_performance_edge_building.py    (new: Phase 1)
│   └── ... (existing tests)
├── PLANNING_DATA_STRUCTURE_MODERNIZATION.md (new: Phase 0)
├── CLAUDE.md                                (new: Phase 0)
├── constitution.md                          (new: Phase 0)
├── project_plan.md                          (new: Phase 0)
└── ... (existing docs)
```

---

## Glossary

| Term | Definition |
|------|-----------|
| CAI | CHILmesh Access Interface (stable bridge APIs) |
| DCEL | Doubly-Connected Edge List (half-edge data structure) |
| MADMESHR | Downstream research project for mesh adaptation |
| ADMESH | Mesh adaptation framework |
| ADMESH-Domains | Domain handling for ADMESH |
| O(n²) | Quadratic time complexity |
| O(n log n) | Linearithmic time complexity |
| Skeletonization | Medial axis extraction via boundary peeling |
| Fort.14 | ADCIRC mesh file format |

---

**Document Version:** 1.0  
**Last Updated:** 2026-04-26  
**Next Review:** End of Phase 0 (Week 4)
