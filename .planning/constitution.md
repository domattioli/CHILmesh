# CHILmesh Governing Constitution

**Effective Date:** 2026-04-26
**Version:** 1.0
**Scope:** Development governance, decision-making, and evolution principles

---

## 1. PREAMBLE

Establishes principles to:
- Ensure sustainable, high-quality development
- Maintain clear decision-making processes
- Preserve scientific integrity and reproducibility
- Enable productive collaboration with downstream projects
- Document principles guiding long-term evolution

Not about implementation details — about values and processes.

---

## 2. CORE PRINCIPLES

### 2.1 Scientific Integrity
- Algorithms rooted in published literature or novel research
- Reproducibility non-negotiable: changes must be testable, benchmarkable
- Performance claims backed by measurements
- Breaking changes documented with rationale

### 2.2 Backward Compatibility
- Public APIs stable until major version bump (v1.0 → v2.0)
- Internal refactoring hidden behind consistent public interfaces
- Deprecations require warning period: mark old methods, support new in parallel
- Tests define contracts: if tests pass, users' code works

### 2.3 Transparency
- Design decisions documented before implementation
- Trade-offs explicit: performance vs memory, complexity vs maintainability
- Rationale preserved in commit messages
- Limitations acknowledged

### 2.4 Pragmatism
- Incremental improvement beats comprehensive redesign
- Measure before optimizing
- One tool per job; avoid over-engineering
- Minimize dependencies: each external library is long-term commitment

### 2.5 Collaborative Integration
- Downstream projects are first-class stakeholders
- Bridge layers explicit and documented
- Integration testing catches incompatibilities early
- Communication with MADMESHR/ADMESH/ADMESH-Domains is proactive

---

## 3. DEVELOPMENT GOVERNANCE

### 3.1 Decision Authority

#### Minor (Implementation Details)
**Authority:** Any active developer with commit access  
**Examples:** variable naming, internal refactoring, test additions, docs  
**Process:** commit with clear message; tests must pass

#### Medium (API Additions, Internal Structure)
**Authority:** Maintainer (Dominik Mattioli) + consensus  
**Examples:** new public methods, adjacency structure changes, optional dependencies, deprecations  
**Process:**
1. Create GitHub issue with `design` label
2. Document: what, why, tradeoffs
3. Solicit feedback (minimum 3 days)
4. Seek consensus; maintainer decides if blocked
5. Record decision before implementation

#### Major (Version Strategy, Mission Scope)
**Authority:** Maintainer + downstream project authors  
**Examples:** v0.2.0 roadmap, major API redesign, large new dependencies  
**Process:**
1. Draft in GitHub issue with `strategic` label
2. Propose to MADMESHR/ADMESH/ADMESH-Domains authors
3. Gather downstream requirements
4. Plan in detail (spec doc)
5. Announce in CHANGELOG and project_plan.md
6. Execute in phases with clear milestones

### 3.2 Breaking Changes

**Definition:** Any change requiring downstream code modifications (API signature, return type, behavior changes)

**Policy:**
1. Avoid until v1.0
2. Document extensively: rationale in GitHub issue; migration guide in CHANGELOG; deprecation warning (minimum 1 version)
3. Coordinate with downstream — notify MADMESHR/ADMESH/ADMESH-Domains
4. Major version bump signals breaking changes

### 3.3 Code Review Standards

**Required for all commits:**
1. Semantic correctness
2. Test coverage (new behavior has tests)
3. Documentation (public APIs have docstrings)
4. Backward compatibility (old tests pass)
5. Performance benchmark if algorithmic change

**When in doubt:**
- Algorithm from paper? → Cite it
- Touches skeletonization? → Comprehensive testing
- Affects I/O? → All fixtures must round-trip
- Slow on large meshes? → Profile and report

---

## 4. QUALITY STANDARDS

### 4.1 Testing Requirements

**Baseline:** all commits pass existing suite; parametrize over all four fixtures; new features require tests; no regressions without documented justification

**Performance Baseline:**
- Annulus adjacency build: <1ms
- Donut adjacency build: <10ms
- Structured adjacency build: <20ms
- Block_O: <45s (currently ~30s, target to improve)
- Skeletonization: <2× adjacency build time
- Fort.14 I/O: <1s all fixtures

**Correctness Invariants:**
- ✅ Signed area > 0 (CCW orientation)
- ✅ All vertices indexed [0, n_verts)
- ✅ All elements indexed [0, n_elems)
- ✅ Connectivity valid (no out-of-range)
- ✅ Skeletonization: disjoint cover of elements
- ✅ Fort.14 roundtrip: mesh → write → read ≡ mesh (geometry equality)

### 4.2 Code Style

**Language:** Python 3.10+  
**Import style:** standard library → numpy/scipy → local modules  
**Type hints:** required for public APIs  
**Comments:** document *why*, not *what*; reference papers/issues for non-obvious algorithms; invariants explicit (docstring or assert)

```python
def skeletonize(self) -> None:
    """
    Extract medial axis via boundary peeling.
    
    Implements layer-based decomposition: iteratively identify boundary elements,
    then their inner neighbors, creating concentric layers outside inward.
    
    Invariant: self.layers['OE'][i] ∩ self.layers['IE'][i] = ∅
               (outer and inner elements disjoint per layer)
    """
```

### 4.3 Documentation Standards

**Public APIs:** docstring with purpose/parameters/returns; type hints; example usage

**Algorithms:** paper reference; complexity analysis; limitations or special cases

**Data structures:** explicit invariants; access patterns; modification rules

```python
"""
Returns Vert2Elem adjacency as explicit dict.
    
Invariant: vert2elem[v] is a Set[int] of all elements incident to vertex v.
           For each elem ∈ vert2elem[v], vertex v appears in connectivity_list[elem].
           
Access: O(1) lookup, O(k) iteration (k = vertex degree)
Update: When adding element e with vertices [v1, v2, v3, v4]:
        for v in [v1, v2, v3, v4]:
            if v not in vert2elem: vert2elem[v] = set()
            vert2elem[v].add(e)
"""
```

---

## 5. RELEASE PROCESS

### 5.1 Versioning

**Format:** MAJOR.MINOR.PATCH (e.g., 0.1.1)
- **MAJOR:** Reserved for v1.0.0 when API stabilizes
- **MINOR:** New features, internal refactoring, performance improvements
- **PATCH:** Bug fixes, tests, documentation

### 5.2 Release Checklist

- [ ] All tests pass on Python 3.10, 3.11, 3.12 × Ubuntu, macOS
- [ ] CHANGELOG.md updated with user-visible changes
- [ ] Version updated in pyproject.toml
- [ ] README examples verified
- [ ] Performance benchmarks recorded
- [ ] GitHub release notes written with migration guide (if needed)
- [ ] PyPI package built + checked with `twine check`

### 5.3 Deprecation Policy

Timeline:
1. Announce in GitHub issue (1 week minimum)
2. Mark with `@warnings.warn(..., DeprecationWarning)`
3. Document in CHANGELOG as "Deprecated in 0.x.y"
4. Support old API for minimum 2 releases
5. Remove in next MINOR version bump

```python
def old_method(self):
    warnings.warn(
        "old_method() is deprecated and will be removed in 0.3.0. "
        "Use new_method() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return self.new_method()
```

---

## 6. COLLABORATION WITH DOWNSTREAM PROJECTS

### 6.1 CHILmesh Commits To

1. **Stability**: public APIs stable until v1.0
2. **Communication**: advance notice of major changes (2+ weeks)
3. **Bridge layer**: clear, documented integration interface
4. **Testing**: integration tests verify downstream compatibility
5. **Support**: help with API questions and migration

### 6.2 Expectations from Downstream

- Use only documented public APIs
- Provide feedback on limitations or missing features
- Coordinate major changes affecting CHILmesh
- Report bugs with reproducible examples
- Join design discussions on strategic changes

### 6.3 Integration Patterns

**Recommended:**
```python
# Downstream uses CHILmesh via stable API
from chilmesh import CHILmesh

mesh = CHILmesh.read_from_fort14(path)
layers = mesh.get_layers()  # Stable, documented API
quality = mesh.elem_quality()
mesh.smooth_mesh(method='fem', acknowledge_change=True)
```

**Discouraged:**
```python
# Directly accessing internal structures (may break)
for edge_id in mesh.adjacencies['Vert2Edge'][v]:
    # This structure changed in 0.2.0!
```

---

## 7. LONG-TERM VISION

### 7.1 Phase Milestones

**0.1.x (Current):** bug fixes; clear docs; bridge architecture design

**0.2.0 (Proposed):** data structure modernization; public bridge APIs; O(n log n) edge building

**0.3.0–0.9.x (Future):** advanced features as needed; spatial indexing; parallel operations

**1.0.0 (Eventually):** API stabilization; full backward compatibility guarantee; mature ecosystem

### 7.2 Success Metrics

- ✅ All tests pass consistently
- ✅ Zero regressions from downstream projects
- ✅ Performance improves monotonically per release
- ✅ Code clarity improves (fewer bugs, easier to understand)
- ✅ Documentation completeness 95%+
- ✅ Community engagement (GitHub issues/discussions active)

### 7.3 Technical Debt Management

**Acceptable:** short-term hacks with clear FIXME; intentional simplifications; legacy code paths flagged for refactor

**Unacceptable:** undocumented behavior changes; silent performance regressions; breaking changes without deprecation; commented-out code

---

## 8. DISPUTE RESOLUTION

### 8.1 Disagreement Process

1. Document both positions in GitHub issue
2. Present evidence: performance data, scientific rationale, downstream impact
3. Seek consensus; iterate until agreement or compromise
4. Escalate to maintainer if consensus not reached within 1 week
5. Maintainer decides with explicit rationale

### 8.2 Reverting Changes

Revert if: tests fail unexpectedly; performance regresses >10% on large fixtures; breaking change not caught; scientific error discovered

**Process:** identify regression → create issue with reproducible example → revert commit → post-mortem → re-implement correctly

---

## 9. AMENDMENTS

**Change Process:**
1. Propose in GitHub issue labeled `constitution`
2. Justify why current principles insufficient
3. Draft new principle/policy
4. Solicit feedback (minimum 1 week)
5. Maintainer approves or suggests refinement
6. Merge to main branch

**Vote Required For:** core principle changes; backward compatibility policy; versioning scheme

**Maintainer Decision For:** governance structure; process improvements; clarifications

---

## 10. DOCUMENT CONTROL

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2026-04-26 | Dominik Mattioli | Initial constitution (Planning Phase) |

---

## External Upstream: DomI

`domattioli/DomI` manages foundational skills and policy used by this repo.
`.domi-pin` committed; session start auto-checks drift via
`scripts/instructions_on_start.sh`. Hard stop on drift; `/sync-from-domi` unblocks.

CHILmesh-specific rules (branch policy, API stability) take precedence over
DomI universal defaults where they conflict. This precedence flows from session-start
read order: `.planning/constitution.md` is read before `.claude/CLAUDE.md`, and
both override DomI universal defaults.

---

## Appendix A: Principles in Action

### Example 1: Fixing Bug in Skeletonization

**Correct:** create issue with reproducible example → write regression test (fails) → analyze root cause → fix with explanation → verify test passes → document in CHANGELOG → communicate to downstream

**Incorrect:** quietly "fix" without tests; change behavior without documentation; break tests and suppress them

### Example 2: Performance Optimization

**Correct:** profile to confirm bottleneck → design new approach → implement with tests verifying same output → benchmark before/after all fixtures → document what changed and why → preserve API → report improvement in release notes

**Incorrect:** change adjacency structure without updating all using code; improve speed but break test semantics; introduce non-determinism

### Example 3: Adding Downstream Integration

**Correct:** ask MADMESHR authors what API they want → document interface → implement with clear semantics → add integration tests → communicate API guarantees → provide migration guide

**Incorrect:** guess without asking; expose internal structures as public API; surprise downstream with incompatibility

---

**This constitution is a living document. It reflects current values and will evolve as CHILmesh matures.**
