# CHILmesh Governing Constitution

**Effective Date:** 2026-04-26
**Version:** 1.0
**Scope:** Development governance, decision-making, and evolution principles

---

## 1. PREAMBLE

We, the developers of CHILmesh, establish this constitution to:
- Ensure sustainable, high-quality development
- Maintain clear decision-making processes
- Preserve scientific integrity and reproducibility
- Enable productive collaboration with downstream projects
- Document principles that guide long-term evolution

This constitution does not dictate implementation details, but rather the values and processes that inform them.

---

## 2. CORE PRINCIPLES

### 2.1 Scientific Integrity
- **All algorithms are rooted in published literature** or novel research
- **Reproducibility is non-negotiable**: changes must be testable, benchmarkable
- **Performance claims are backed by measurements**, not intuition
- **Breaking changes are documented** with rationale, not silently introduced

### 2.2 Backward Compatibility
- **Public APIs are stable** until major version bump (v1.0 → v2.0)
- **Internal refactoring is hidden** behind consistent public interfaces
- **Deprecations require warning period**: mark old methods, support new ones in parallel
- **Tests define contracts**: if tests pass, users' code works

### 2.3 Transparency
- **Design decisions are documented** before implementation
- **Trade-offs are explicit**: performance vs memory, complexity vs maintainability
- **Rationale is preserved** in commit messages, not discarded
- **Limitations are acknowledged**: no feature pretending to be universal

### 2.4 Pragmatism
- **Perfect is the enemy of good**: incremental improvement beats comprehensive redesign
- **Measure before optimizing**: profile, don't guess
- **One tool per job**: avoid over-engineering or premature abstraction
- **Minimize dependencies**: each external library is a long-term commitment

### 2.5 Collaborative Integration
- **Downstream projects are first-class stakeholders**
- **Bridge layers are explicit and documented**
- **Integration testing catches incompatibilities early**
- **Communication with MADMESHR/ADMESH/ADMESH-Domains is proactive**

---

## 3. DEVELOPMENT GOVERNANCE

### 3.1 Decision Authority

#### Minor Decisions (Implementation Details)
**Authority:** Any active developer with commit access
**Examples:**
- Variable naming, code organization
- Internal refactoring (same semantics, different implementation)
- Test additions
- Documentation improvements

**Process:**
- Commit with clear message
- Reference related issue if applicable
- Tests must pass

#### Medium Decisions (API Additions, Internal Structure)
**Authority:** Project maintainer (Dominik Mattioli) + consensus with active developers
**Examples:**
- New public method signatures
- Changing adjacency structure representations
- Adding optional dependencies
- Deprecating old features

**Process:**
1. Create GitHub issue with `design` label
2. Document in comment: what, why, tradeoffs
3. Solicit feedback (minimum 3 days)
4. Seek consensus; if blocked, maintainer decides
5. Record decision in issue before implementation

#### Major Decisions (Version Strategy, Mission Scope)
**Authority:** Maintainer + consultation with downstream project authors
**Examples:**
- v0.2.0 roadmap and timeline
- Major API redesign
- Accepting large new dependencies
- Changing mesh element types supported

**Process:**
1. Draft in GitHub issue with `strategic` label
2. Propose to MADMESHR/ADMESH/ADMESH-Domains authors
3. Gather requirements from downstream projects
4. Plan in detail (e.g., PLANNING_DATA_STRUCTURE_MODERNIZATION.md)
5. Announce in CHANGELOG and project_plan.md
6. Execute in phases with clear milestones

### 3.2 Breaking Changes

**Definition:** Any change requiring downstream code modifications
- API signature changes
- Return type changes
- Behavior changes in existing methods

**Policy:**
1. **Avoid until v1.0**: minimize breaking changes in 0.x releases
2. **Document extensively** if necessary:
   - Rationale in GitHub issue
   - Migration guide in CHANGELOG
   - Deprecation warning in old code (minimum 1 version)
3. **Coordinate with downstream**: notify MADMESHR/ADMESH/ADMESH-Domains authors
4. **Major version bump** signals breaking changes (0.1.1 → 1.0.0)

### 3.3 Code Review Standards

**Required for all commits:**
1. **Semantic correctness**: code does what it claims
2. **Test coverage**: new behavior has corresponding tests
3. **Documentation**: public APIs have docstrings
4. **Backward compatibility**: old tests still pass
5. **Performance**: benchmark reported if algorithmic change

**When in doubt, ask:**
- Is this a scientific paper's algorithm? → Cite it
- Does this change touch skeletonization? → Comprehensive testing
- Does this affect I/O? → All fixtures must round-trip
- Is this slow on large meshes? → Profile and report

---

## 4. QUALITY STANDARDS

### 4.1 Testing Requirements

**Baseline:**
- All commits must pass existing test suite
- Parametrize over all four built-in fixtures
- New features require new tests (test-driven preferred)
- No test regressions without documented justification

**Performance Baseline:**
- Adjacency building:
  - Annulus: < 1ms
  - Donut: < 10ms
  - Structured: < 20ms
  - Block_O: < 45s (currently ~30s, target is to improve)
- Skeletonization: < 2x adjacency build time
- Fort.14 I/O: < 1s for all fixtures

**Correctness Invariants:**
- ✅ Signed area > 0 (counter-clockwise orientation)
- ✅ All vertices indexed [0, n_verts)
- ✅ All elements indexed [0, n_elems)
- ✅ Connectivity valid (no out-of-range indices)
- ✅ Skeletonization: disjoint cover of elements
- ✅ Fort.14 roundtrip: mesh → write → read ≡ mesh (geometry equality)

### 4.2 Code Style

**Language:** Python 3.10+
**Import style:** Standard library → numpy/scipy → local modules
**Type hints:** Required for public APIs
**Comments:**
- Document *why*, not *what*
- Avoid comments that restate code
- Reference papers/issues for non-obvious algorithms
- Invariants should be explicit (docstring or assert)

**Example:**
```python
def skeletonize(self) -> None:
    """
    Extract medial axis via boundary peeling.
    
    Implements the layer-based decomposition from Mattioli's thesis:
    iteratively identify boundary elements, then their inner neighbors,
    creating concentric layers from outside inward.
    
    Invariant: self.layers['OE'][i] ∩ self.layers['IE'][i] = ∅
               (outer and inner elements disjoint per layer)
    """
```

### 4.3 Documentation Standards

**Required for all public APIs:**
- Docstring with purpose, parameters, returns
- Type hints (e.g., `def foo(x: ndarray) -> Set[int]`)
- Example usage in docstring or test

**Required for algorithms:**
- Paper reference (if applicable)
- Complexity analysis (e.g., O(n log n))
- Limitations or special cases

**Required for data structures:**
- Explicit invariants (what must always be true)
- Access patterns (how to query efficiently)
- Modification rules (how to maintain invariants when changing)

**Example:**
```python
# In docstring:
"""
Returns Vert2Elem adjacency as explicit dict.
    
Invariant: vert2elem[v] is a Set[int] containing all elements incident to vertex v.
           For each elem ∈ vert2elem[v], vertex v must appear in 
           self.connectivity_list[elem].
           
Access: O(1) lookup, O(k) iteration (k = degree of vertex)
Update: When adding element e with vertices [v1, v2, v3, v4]:
        for v in [v1, v2, v3, v4]:
            if v not in vert2elem: vert2elem[v] = set()
            vert2elem[v].add(e)
"""
```

---

## 5. RELEASE PROCESS

### 5.1 Versioning Scheme

**Format:** MAJOR.MINOR.PATCH (e.g., 0.1.1)
- **MAJOR:** Reserved for v1.0.0 when API stabilizes
- **MINOR:** New features, internal refactoring, performance improvements
- **PATCH:** Bug fixes, test improvements, documentation

### 5.2 Release Checklist

Before releasing version X.Y.Z:
- [ ] All tests pass on Python 3.10, 3.11, 3.12
- [ ] All tests pass on Ubuntu, macOS
- [ ] CHANGELOG.md is updated with user-visible changes
- [ ] Version number updated in pyproject.toml
- [ ] README examples are verified to run and produce expected output
- [ ] Performance benchmarks recorded (adjacency, skeletonization, fort.14)
- [ ] GitHub release notes written with migration guide (if needed)
- [ ] PyPI package built and checked with `twine check`

### 5.3 Deprecation Policy

**Deprecation Timeline:**
1. Announce in GitHub issue (1 week minimum discussion)
2. Mark with `@warnings.warn(..., DeprecationWarning)` in code
3. Document in CHANGELOG as "Deprecated in 0.x.y"
4. Support old API for minimum 2 releases
5. Remove in next MINOR version bump

**Example:**
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

### 6.1 Commitment

CHILmesh commits to:
1. **Stability**: Public APIs stable until v1.0
2. **Communication**: Advance notice of major changes (2+ weeks)
3. **Bridge layer**: Clear, documented interface for integration
4. **Testing**: Integration tests verify downstream compatibility
5. **Support**: Help with API questions and migration issues

### 6.2 Expectations from Downstream Projects

- Use only documented public APIs
- Provide feedback on limitations or missing features
- Coordinate major changes affecting CHILmesh
- Report bugs with reproducible examples
- Join design discussions on strategic changes

### 6.3 Integration Patterns

**Recommended:**
```python
# Downstream project uses CHILmesh via stable API
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

**0.1.x (Current)**
- Bug fixes and hardening
- Clear documentation of existing functionality
- Bridge architecture design

**0.2.0 (Proposed)**
- Data structure modernization (Phase 1-3 of PLANNING_DATA_STRUCTURE_MODERNIZATION.md)
- Public bridge APIs for MADMESHR/ADMESH/ADMESH-Domains
- Improved performance (O(n²) → O(n log n) on edge building)

**0.3.0-0.9.x (Future)**
- Advanced features as needed by downstream projects
- Spatial indexing if required
- Parallel operations (mesh refinement, smoothing)

**1.0.0 (Eventually)**
- API stabilization
- Full backward compatibility guarantee
- Mature ecosystem and integrations

### 7.2 Success Metrics

- ✅ All tests pass consistently
- ✅ Zero regressions reported from downstream projects
- ✅ Performance improves monotonically with each release
- ✅ Code clarity improves (fewer bugs, easier to understand)
- ✅ Documentation completeness reaches 95%+
- ✅ Community engagement (GitHub issues/discussions active)

### 7.3 Technical Debt Management

**Acceptable:**
- Short-term hacks with clear FIXME comments
- Intentional simplifications that can be optimized later
- Legacy code paths that will be refactored in next version

**Unacceptable:**
- Undocumented behavior changes
- Silent performance regressions
- Breaking changes without deprecation period
- Commented-out code (delete it)

---

## 8. DISPUTE RESOLUTION

### 8.1 Disagreement Process

If developers disagree on technical direction:
1. **Document both positions** in GitHub issue
2. **Present evidence**: performance data, scientific rationale, downstream impact
3. **Seek consensus**: iterate until agreement or compromise
4. **Escalate to maintainer** if consensus not reached within 1 week
5. **Maintainer decides** with explicit rationale in issue

### 8.2 Reverting Changes

Changes can be reverted if:
- Tests fail unexpectedly
- Performance regresses >10% on large fixtures
- Breaking change not caught in review
- Scientific error discovered (algorithmic flaw)

**Process:**
1. Identify regression
2. Create issue with reproducible example
3. Revert commit
4. Post-mortem: improve tests or review process
5. Re-implement correctly

---

## 9. AMENDMENTS TO CONSTITUTION

**Change Process:**
1. Propose amendment in GitHub issue labeled `constitution`
2. Justify why current principles are insufficient
3. Draft new principle/policy
4. Solicit feedback (minimum 1 week)
5. Maintainer approves or suggests refinement
6. Merge amended constitution.md to main branch

**Vote Required For:**
- Changes to core principles (Section 2)
- Changes to backward compatibility policy
- Changes to release/versioning scheme

**Maintainer Decision For:**
- Governance structure changes
- Process improvements
- Clarifications of existing principles

---

## 10. DOCUMENT CONTROL

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2026-04-26 | Dominik Mattioli | Initial constitution (Planning Phase) |

---

## Appendix A: Principles in Action

### Example 1: Fixing Bug in Skeletonization

**Scenario:** Tests discover that skeletonization produces overlapping layers.

**Correct Process:**
1. Create GitHub issue with minimal reproducible example
2. Write regression test (fails with current code)
3. Analyze root cause
4. Fix algorithm with explanation
5. Verify test now passes
6. Document in CHANGELOG and commit message
7. Communicate fix to downstream projects if they rely on this behavior

**Incorrect Process:**
- Quietly "fix" without tests
- Change behavior without documentation
- Break tests and suppress them

### Example 2: Performance Optimization

**Scenario:** O(n²) edge building is too slow on large meshes.

**Correct Process:**
1. Profile to confirm bottleneck
2. Design new approach (hash maps, sorted arrays, etc.)
3. Implement with tests verifying same output
4. Benchmark before/after on all fixtures
5. Document what changed and why (in code and commit message)
6. Preserve API (internal refactoring, not user-visible)
7. Report performance improvement in release notes

**Incorrect Process:**
- Change adjacency structure without updating all code using it
- Improve speed but break test semantics
- Introduce randomness or non-determinism
- Hide changes behind opaque "optimization"

### Example 3: Adding Downstream Integration

**Scenario:** MADMESHR needs fast neighbor queries.

**Correct Process:**
1. Ask MADMESHR authors: what API would you like?
2. Document interface in PLANNING document
3. Implement with clear semantics
4. Add integration tests simulating MADMESHR usage
5. Communicate API guarantees (what will never change)
6. Provide migration guide if replacing old approach

**Incorrect Process:**
- Guess what they need without asking
- Make public API directly expose internal structures
- Change behavior based on implicit understanding
- Surprise them with incompatibility in next release

---

**This constitution is a living document. It reflects our current values and will evolve as CHILmesh matures.**
