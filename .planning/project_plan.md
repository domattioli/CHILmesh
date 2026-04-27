# CHILmesh Project Plan & Roadmap

**Current Version:** 0.1.1 (Alpha)
**Planning Date:** 2026-04-26
**Horizon:** 12 months (through Q1 2027)

---

## Executive Summary

CHILmesh is transitioning from a bug-fix phase (0.1.1) to a **data structure modernization and bridge architecture phase** (0.2.0). This plan outlines:

1. **Immediate priorities** (Weeks 1-4): GitHub issues, stakeholder coordination
2. **Development phases** (Months 2-6): Three phases of architecture improvements
3. **Stabilization** (Months 7-9): Testing, performance validation, documentation
4. **Release & communication** (Months 10-12): v0.2.0 release, integration with downstream projects

**Key Success Criteria:**
- ✅ 1.5x+ performance improvement on large meshes
- ✅ Zero breaking changes to public API
- ✅ Clear bridge interfaces for MADMESHR/ADMESH/ADMESH-Domains
- ✅ 100% test pass rate on all platforms
- ✅ Complete documentation of architecture decisions

---

## Phase 0: Planning & Coordination (CURRENT - Weeks 1-4)

### Objectives
1. Create detailed specification for architectural changes
2. Coordinate with downstream project authors
3. Set up GitHub issues for tracking
4. Establish governance principles

### Deliverables
- ✅ `PLANNING_DATA_STRUCTURE_MODERNIZATION.md` (comprehensive spec)
- ✅ `CLAUDE.md` (developer guide for AI assistance)
- ✅ `constitution.md` (governance principles)
- ✅ `project_plan.md` (this document)
- 🔄 GitHub issues (P1-01 through P3-04, 16 total)
- 🔄 Communication with MADMESHR/ADMESH/ADMESH-Domains authors

### Timeline
- **Week 1-2**: Complete planning documents ← **You are here**
- **Week 3**: Publish to main branch, create GitHub issues
- **Week 4**: Stakeholder feedback, adjust timeline if needed

---

## Phase 1: Hash Map Edge Lookup (Weeks 5-8)

**Focus:** Eliminate O(n²) bottleneck in adjacency building
**Estimated Effort:** 8-12 hours
**Risk Level:** Low

### Tasks
| ID | Task | Owner | Hours | Dependencies |
|----|----|-------|-------|--------------|
| P1-01 | Create EdgeMap class | TBD | 2 | None |
| P1-02 | Refactor _identify_edges() | TBD | 2 | P1-01 |
| P1-03 | Optimize _build_elem2edge() | TBD | 2 | P1-02 |
| P1-04 | Optimize _build_edge2elem() | TBD | 2 | P1-02 |
| P1-05 | Store EdgeMap in adjacencies | TBD | 1 | P1-03, P1-04 |
| P1-06 | Performance regression tests | TBD | 3 | P1-05 |

### Success Criteria
- All existing tests pass
- Block_O build time: <45s (currently ~30-60s range)
- Benchmarks show ≥1.5x improvement on large fixtures
- No new dependencies

### Outputs
- `src/chilmesh/mesh_topology.py` (new EdgeMap class)
- Updated `_build_adjacencies()` using EdgeMap
- Performance benchmark report
- All tests passing

---

## Phase 2: Adjacency Modernization (Weeks 9-14)

**Focus:** Modernize sparse adjacency structures (Vert2Edge, Vert2Elem)
**Estimated Effort:** 12-16 hours
**Risk Level:** Medium

### Tasks
| ID | Task | Owner | Hours | Dependencies |
|----|----|-------|-------|--------------|
| P2-01 | Migrate Vert2Edge to dict | TBD | 2 | P1-06 |
| P2-02 | Migrate Vert2Elem to dict | TBD | 2 | P1-06 |
| P2-03 | Add explicit accessor methods | TBD | 2 | P2-01, P2-02 |
| P2-04 | Update traversal patterns | TBD | 4 | P2-03 |
| P2-05 | Type hints & validation | TBD | 3 | P2-04 |
| P2-06 | Adjacency documentation | TBD | 3 | P2-05 |

### Success Criteria
- All tests still passing
- Vert2Edge and Vert2Elem explicitly documented as Dict[int, Set[int]]
- Type hints complete for adjacency methods
- Traversal patterns updated (no list-of-lists iteration)
- Performance at least as good as Phase 1

### Outputs
- Refactored `_build_vert2edge()` and `_build_vert2elem()`
- New public methods: `get_vertex_edges()`, `get_vertex_elements()`
- Comprehensive adjacency structure documentation
- Type hint coverage report

---

## Phase 3: Bridge Infrastructure (Weeks 15-26)

**Focus:** Define and implement clean bridge APIs for downstream projects
**Estimated Effort:** 20-24 hours
**Risk Level:** High (depends on downstream project coordination)

### Tasks
| ID | Task | Owner | Hours | Dependencies |
|----|----|-------|-------|--------------|
| P3-01 | Define CAI (CHILmesh Access Interface) | TBD | 4 | P2-06 |
| P3-02 | Create bridge adapters | TBD | 10 | P3-01 |
| P3-03 | Integration tests | TBD | 6 | P3-02 |
| P3-04 | Downstream migration guide | TBD | 4 | P3-03 |

### Success Criteria
- CAI document published and stable
- Bridge adapters for MADMESHR, ADMESH, ADMESH-Domains
- Integration tests covering common usage patterns
- Clear migration path for downstream projects
- Zero external dependency on internal implementation

### Outputs
- `docs/CHILmesh_Access_Interface.md` (CAI specification)
- `src/chilmesh/bridge.py` (adapter implementations)
- `tests/test_bridge_integration.py` (integration tests)
- `docs/DOWNSTREAM_MIGRATION_GUIDE.md` (migration documentation)

### Coordination Notes
- Requires early feedback from MADMESHR/ADMESH/ADMESH-Domains authors
- May need to adjust timeline based on downstream project needs
- Should include calls/meetings to clarify integration patterns

---

## Phase 4: Stabilization & Testing (Weeks 27-36)

**Focus:** Comprehensive testing, performance validation, final hardening
**Estimated Effort:** 16-20 hours
**Risk Level:** Low (mostly verification)

### Activities
1. **Full Regression Testing**
   - Run all tests on Python 3.10, 3.11, 3.12
   - Test on Ubuntu, macOS, Windows
   - Verify fort.14 roundtrip on all fixtures

2. **Performance Benchmarking**
   - Baseline before Phase 1, after Phase 1, after Phase 2
   - Compare to v0.1.1 on all fixtures
   - Document improvements in release notes

3. **Memory Profiling**
   - Check memory usage hasn't increased >10%
   - Identify any leaks in new structures

4. **Documentation Review**
   - Update README with new capabilities
   - Verify all examples still work
   - Add architecture documentation

5. **Community Feedback**
   - Share release candidate with MADMESHR/ADMESH/ADMESH-Domains
   - Gather feedback on bridge APIs
   - Address issues before final release

### Deliverables
- Test coverage report (target: >90%)
- Performance benchmark report
- Memory profiling results
- Updated README and documentation
- Resolved feedback from downstream projects

---

## Phase 5: Release & Communication (Weeks 37-40)

**Focus:** Final validation, release, and ecosystem communication
**Estimated Effort:** 4-6 hours

### Release Checklist
- [ ] All tests pass on all platforms
- [ ] Performance benchmarks meet targets
- [ ] CHANGELOG.md updated with all changes
- [ ] Version bumped to 0.2.0 in pyproject.toml
- [ ] GitHub release notes prepared
- [ ] PyPI build validated with `twine check`
- [ ] Migration guide ready for downstream projects
- [ ] Announce in CHIL Lab channels
- [ ] Tag v0.2.0 in git
- [ ] Publish to PyPI

### Communication
- GitHub Release: Detailed notes on Phase 1-3 improvements
- Email to MADMESHR/ADMESH/ADMESH-Domains authors: Integration guide
- Lab announcement: What changed and why
- Update MATLAB documentation (if applicable)

---

## Post-Release: Future Planning (Months 10-12)

### Immediate (v0.2.1-0.2.x patches)
- Bug fixes from community feedback
- Documentation improvements
- Minor performance enhancements

### Medium-term (v0.3.0 exploration)
- Gather requirements from downstream projects
- Evaluate spatial indexing needs
- Consider mesh mutation APIs (add/remove nodes/edges)
- Profile memory usage on largest real-world meshes

### Long-term (v1.0.0 vision)
- Assess when API can be frozen
- Evaluate GPU acceleration for smoothing/quality
- Consider distributed mesh handling
- Plan for external mesh format support (Gmsh, VTK)

---

## Resource Allocation

### Team
- **Lead Developer:** Claude Code (AI-assisted development)
- **Code Review:** Dominik Mattioli (scientific oversight)
- **Testing:** Automated via pytest CI/CD
- **Documentation:** Embedded in code (docstrings) + standalone guides

### Tools & Infrastructure
- GitHub Issues/Projects for task tracking
- GitHub Actions for CI/CD (Python 3.10, 3.11, 3.12)
- PyPI for package distribution
- pytest for test automation

### External Dependencies
- **numpy** (required): dense numerical operations
- **scipy** (required): sparse matrices, spatial structures
- **matplotlib** (required): visualization
- **pytest** (dev): testing framework

**No new external dependencies planned for 0.2.0**

---

## Risk Management

### High-Risk Items

| Risk | Impact | Mitigation |
|------|--------|-----------|
| Skeletonization breaks | Critical | Comprehensive tests, careful refactoring |
| Downstream incompatibility | Critical | Early coordination, integration tests, clear APIs |
| Performance regression | High | Benchmark every step, revert if slower |
| Complex merge conflicts | High | Clear branch strategy, small focused PRs |

### Medium-Risk Items

| Risk | Impact | Mitigation |
|------|--------|-----------|
| Scope creep | Medium | Strict phase boundaries, issue templates |
| Timeline slip | Medium | Realistic estimates, track velocity |
| Technical debt | Medium | Design reviews, architecture documentation |
| Breaking changes leak | Medium | Automated compatibility tests, code review |

### Low-Risk Items
- Small bug fixes during development
- Documentation updates
- Test improvements
- Type hint additions

---

## Success Metrics

### Quantitative
- ✅ **Performance**: 1.5x+ improvement on adjacency building
- ✅ **Coverage**: 100% test pass rate on all platforms
- ✅ **Compatibility**: Zero breaking changes to public API
- ✅ **Documentation**: 95%+ of public methods have docstrings
- ✅ **Integration**: 3+ downstream projects successfully integrating

### Qualitative
- ✅ **Code clarity**: Easier to understand adjacency structures
- ✅ **Maintainability**: Fewer open issues related to complexity
- ✅ **Stakeholder satisfaction**: Positive feedback from downstream authors
- ✅ **Scientific integrity**: All changes properly documented and tested
- ✅ **Community engagement**: Active GitHub discussions, clear communication

---

## Dependency on External Knowledge

### Known Dependencies
- **MADMESHR internals**: Need to understand what APIs it expects
- **ADMESH requirements**: Mesh adaptation strategies it uses
- **ADMESH-Domains patterns**: Domain handling and tagging

### Action Items
1. **Week 2-3**: Reach out to all three project authors
2. **Week 3-4**: Schedule calls to clarify integration needs
3. **Week 5+**: Incorporate feedback into Phase 3 design

---

## Governance & Decision-Making

**Principle:** This plan is governed by `constitution.md`

### Decision Authority
- **Minor changes** (task details, test improvements): Any active developer
- **Phase scope adjustments**: Maintainer + consensus
- **Timeline changes**: Maintainer + stakeholder communication
- **API additions**: Design review + downstream coordination

### Escalation
1. Start in GitHub issue with context
2. Request feedback from stakeholders
3. Escalate to maintainer if consensus not reached
4. Document decision and rationale in issue

---

## Appendix A: Glossary of Terms

| Term | Definition |
|------|-----------|
| CAI | CHILmesh Access Interface (stable bridge APIs) |
| DCEL | Doubly-Connected Edge List (half-edge data structure) |
| MADMESHR | Downstream research project for mesh adaptation |
| ADMESH | Mesh adaptation framework |
| ADMESH-Domains | Domain handling for ADMESH |
| O(n²) | Quadratic time complexity (slow on large n) |
| O(n log n) | Linearithmic time complexity (efficient) |
| Skeletonization | Medial axis extraction via boundary peeling |
| Fort.14 | ADCIRC mesh file format |

---

## Appendix B: Phase Summary Chart

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

## Appendix C: File Structure After Modernization

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

**Document Version:** 1.0
**Last Updated:** 2026-04-26
**Next Review:** End of Phase 0 (Week 4)
