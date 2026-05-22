# CHILmesh Constitution

**Version:** 1.1 | **Adopted:** 2026-05-15 | **Last Amended:** 2026-05-18 | **Status:** Governance Contract

This is the canonical project constitution. The older `.planning/constitution.md`
is preserved as a redirect stub. The spec-structure guide at
`.specify/speckit-constitution.md` is a separate document about how to lay out
individual feature specs (`spec.md`, `plan.md`, `tasks.md` etc.) and is NOT a
project-principles document.

---

## Core Principles

### I. Library-First Architectural Principle

Every feature starts as a standalone library before integration. Libraries must be:
- **Self-contained:** No circular dependencies; clear public API surface
- **Independently testable:** Unit tests runnable without external mesh data
- **Documented:** Docstrings for all public methods; usage examples in README

**Rationale:** CHILmesh is a mature mesh library (not a one-off tool). Users import specific modules (`from chilmesh import CHILmesh, write_fort14`). Architectural decisions must defend that contract.

### II. Mesh Immutability by Default (Post-v1.0)

Meshes are read-only post-initialization. Topology changes (`split`, `merge`, `insert_vertex`) are explicitly opt-in via future **Phase 5: Mutation API** and NOT the default path.

**Current state:** CHILmesh is mutation-free (read-only). Vertex smoothing works because it modifies only coordinates, not topology.

**Exception:** Boundary vertex coordinate motion (smoothing). Topology invariant: element count, edge count, layer structure remain constant.

**Rationale:** Immutability = deterministic behavior, cache-friendly, parallel-safe. Mutation is post-v1.0; design separately per #94, #93.

### III. Test-First + Regression-Gated Releases

- **TDD mandatory for public API:** Tests written → User/maintainer approved → Implement
- **Regression tests required:** Every bug fix includes test that would have caught it
- **Release gate:** All tests passing before tag; CI runs on every PR

**Test coverage expected:** ≥90% on public API surface (core CHILmesh class). Visualization (matplotlib) is smoke-test only (rendering hard to unit-test).

**Rationale:** CHILmesh ships downstream to MADMESHR, ADMESH. Reliability is non-negotiable.

### IV. Geometric Correctness Over Performance (v1.0)

When correctness conflicts with speed, choose correctness first. Optimize after verification.

**Examples:**
- Floating-point tolerances must be explicit, documented, user-adjustable
- Mixed-element support (triangles + quads) is correctness-mandatory even if slower
- MATLAB parity tests (e.g., skeletonization layer counts) are regression gates

**Rationale:** Physics-based simulations (coastal modeling, CFD) amplify small errors.

### V. Coordinate System Agnosticism

CHILmesh does not assume lat/lon (geographic) OR Cartesian. Coordinates are opaque:
- `bounding_box` returns `{min_x, max_x, min_y, max_y}` (axis labels only)
- Semantic interpretation (lon/lat vs. UTM vs. local Cartesian) is caller's responsibility

**Rationale:** Users bring meshes from diverse sources (ADCIRC, ADMESH, Triangle, GMsh, SMS). Forcing geographic = false choice.

### VI. Format Pluralism: ADCIRC + SMS 2DM Support

CHILmesh reads/writes multiple formats without privileging one:
- **ADCIRC fort.14:** Primary format; full read/write support
- **SMS 2DM:** Secondary; read support mandatory, write support future
- **Others:** Via adapter pattern

**Rationale:** Coastal modeling ecosystem is fragmented. Forcing a single format alienates users.

### VII. Public API Stability

Once a method enters the public API (in `__init__.py`), changes require:
1. **Deprecation cycle:** 1 minor version before removal
2. **Semantic versioning:** MAJOR.MINOR.PATCH (no 0.x forever)
3. **Changelog entry:** Must document breaking changes, deprecations, migration path

**Rationale:** Downstream users (MADMESHR, ADMESH) depend on stable imports.

### VIII. Documentation = Contract

Public API documentation is normative (not descriptive):
- Docstring says what method does; implementation must match
- Examples in docstring are tested (doctests or linked test cases)
- If docstring and implementation disagree → docstring wins; fix implementation

**Rationale:** Documentation is the contract between library and user. Inconsistency erodes trust.

### IX. Mixed-Element Matrices + Correctness

Mixed-element support (triangles + quads) is a **core competency**, not a nice-to-have:
- Tests must cover all element combinations (tri-only, quad-only, mixed)
- Smoothing, quality metrics, skeletonization handle both types
- Schema stores `element_type` as metadata; reader enforces consistency

**Why:** Quad-dominant meshes (Gulf of Mexico, large bays) are common in coastal modeling.

### X. MATLAB Parity Where It Matters

Original QuADMesh+ (MATLAB) is reference for:
- Skeletonization layer definitions (layer assignment algorithm must match within tolerance)
- Element quality metrics (Jacobian, aspect ratio, scaled Jacobian)
- Mixed-element handling (if MATLAB version had it)

**Non-parity is acceptable:** Performance optimizations (vectorization), API design (Pythonic naming), internal data structures.

**Rationale:** CHILmesh is a Python port, not a rewrite. Geometric correctness is the primary reason for the port.

## Governance

### Decision Authority

| Scope | Authority | Process |
|---|---|---|
| Minor (naming, internal refactor, tests, docs) | Any contributor with commit access | Commit with clear message; tests must pass |
| Medium (API additions, internal structure, deprecations) | Maintainer + consensus | GH issue with `design` label; document tradeoffs; minimum 3-day feedback window |
| Major (version strategy, mission scope, large dependencies) | Maintainer + downstream authors (MADMESHR/ADMESH/ADMESH-Domains) | GH issue with `strategic` label; gather downstream requirements; phased plan recorded in `project_plan.md` and CHANGELOG |

### Breaking Changes

Any change requiring downstream code modification. Policy:
1. Avoid until v1.0; opt-in via Phase-5 mutation APIs (#94) after v1.0.
2. Document: rationale in GH issue, migration guide in CHANGELOG, deprecation warning for one minor version.
3. Coordinate: notify MADMESHR/ADMESH/ADMESH-Domains before merging.
4. Major version bump signals breaking changes.

### Amendments

1. **Justification:** Why existing principle fails or is misaligned with project goals.
2. **Impact analysis:** Which open issues/PRs affected?
3. **Ratification:** Maintainer approval (commit + CHANGELOG entry).
4. **Migration plan:** How do existing code/tests adapt?

### Enforcement

- **Code review checklist:** All PRs audited against constitution.
- **Release gate:** Release notes must demonstrate compliance.
- **Conflict resolution:** If principle conflicts with deadline, maintain principle; extend timeline.

### Dispute Resolution

1. Document both positions in a GH issue with evidence (performance data, scientific rationale, downstream impact).
2. Seek consensus; iterate.
3. Escalate to maintainer if unresolved within 1 week.
4. Maintainer decides with explicit rationale recorded on the issue.

**Revert criteria:** unexpected test failure, performance regression >10% on large fixtures, breaking change not caught, scientific error discovered. Process: open regression issue with reproducer → revert commit → post-mortem → re-implement correctly.

---

## Release Process

### Versioning

`MAJOR.MINOR.PATCH` (semantic versioning).

- **MAJOR:** Reserved for v1.0.0 when API stabilizes.
- **MINOR:** New features, internal refactoring, performance improvements.
- **PATCH:** Bug fixes, tests, documentation.

### Release Checklist

- [ ] All tests pass on Python 3.10 / 3.11 / 3.12 × Ubuntu / macOS / Windows.
- [ ] `CHANGELOG.md` updated with user-visible changes.
- [ ] Version bumped in `pyproject.toml`.
- [ ] README examples verified.
- [ ] Performance benchmarks recorded.
- [ ] GitHub release notes written, including migration guide if needed.
- [ ] PyPI package built and checked with `twine check`.

### Deprecation Timeline

1. Announce in GH issue (1 week minimum).
2. Mark with `warnings.warn(..., DeprecationWarning, stacklevel=2)`.
3. Document in CHANGELOG as "Deprecated in 0.x.y".
4. Support old API for minimum 2 releases.
5. Remove in next MINOR version bump.

### Performance Baselines (v0.2.0)

- Annulus adjacency build: <1 ms
- Donut adjacency build: <10 ms
- Structured adjacency build: <20 ms
- Block_O full initialization: ~14 s (targeting further reduction in v0.3.x)
- Skeletonization: <2× adjacency build time
- Fort.14 I/O: <1 s on all fixtures

---

## Feature-Specific Principles

### I/O and Data Loading

- **Zero-regression I/O:** Fort.14 reader/writer must preserve mesh through roundtrip (geometry equality on all fixtures).
- **Atomic format support:** Quad-element support shipped with triangle support; new format readers (SMS 2DM, etc.) route through the existing public API without behavior changes for existing code.
- **Boundary preservation:** Operations affecting mesh coordinates (ADMESH warm-start, optimization) verify boundary vertices via `np.array_equal()` (bit-exact, not approximate) and document this in their contracts.

### Mesh Smoothing

- Public API `smooth_mesh(method='fem')` signature is stable.
- All element types supported: pure triangles, pure quads, mixed (padded triangles + quads).
- Boundary pinning is exact; interior nodes converge toward equilibrium.
- No element collapse: minimum element quality is preserved across all fixtures.
- **DOMsmooth hybrid fallback** (issues #95, #100): when FEM assembly on a mixed mesh degrades quality, fall back to (1) split tri-only/quad-only submeshes, (2) FEM-smooth each, (3) recombine, (4) angle-based polish (Zhou & Shimada). The public dispatcher selects FEM vs DOMsmooth automatically based on mesh composition.

### Skeletonization

- **MATLAB fidelity:** Layer separation invariant — a vertex in layer k does not appear in layers k+2 or beyond. Source line references to the MATLAB reference are preserved in comments.
- **Algorithm stability:** Skeletonization preserves invariants across triangular, quad, and mixed meshes. Tests verify layer disjoint cover and vertex-layer assignment on every fixture.

---

## External Upstream: DomI

`domattioli/DomI` manages foundational skills and policy used by this repo. `.domi-pin` is committed; session start runs `scripts/instructions_on_start.sh`, which invokes `sync-from-domi/check_pin.sh` and HARD STOPs on drift. `/sync-from-domi` unblocks.

CHILmesh-specific rules (branch policy, API stability) take precedence over DomI universal defaults where they conflict. The precedence flows from the documented session-start read order: `.specify/memory/constitution.md` is read before `.claude/CLAUDE.md`, and both override DomI universal defaults.

---

## Living Document

This constitution reflects current project maturity (v1.0 finalization). As CHILmesh grows:
- **v1.0 → v2.0:** Mutation API will add a new principle (transactional topology edits — see #94).
- **Feedback cycles:** Users discover gaps; amendment proposed + ratified.
- **Versioning:** Constitution version bumped with each amendment.

---

## Document Control

| Version | Date | Changes |
|---|---|---|
| 1.0 | 2026-05-15 | Initial speckit-format constitution: 10 core principles + governance/amendment process. |
| 1.1 | 2026-05-18 | Consolidated operationally-useful content from `.planning/constitution.md` (decision authority matrix, breaking-change policy, dispute resolution, release process, deprecation timeline, performance baselines, feature-specific principles, DomI section). Resolves #107. |

---

**Version**: 1.1 | **Ratified**: 2026-05-15 | **Last Amended**: 2026-05-18
