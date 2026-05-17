# CHILmesh Constitution

**Version:** 1.0 | **Adopted:** 2026-05-15 | **Status:** Governance Contract

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

### Amendments

Constitution amendments require:
1. **Justification:** Why existing principle fails or is misaligned with project goals
2. **Impact analysis:** Which open issues/PRs affected?
3. **Ratification:** Maintainer approval (documented in commit + CHANGELOG)
4. **Migration plan:** How do existing code/tests adapt?

### Enforcement

- **Code review checklist:** All PRs audited against constitution
- **Release gate:** Release notes must demonstrate compliance
- **Conflict resolution:** If principle conflicts with deadline, maintain principle; extend timeline

### Living Document

This constitution reflects current project maturity (v1.0 finalization). As CHILmesh grows:
- **v1.0 → v2.0:** Mutation API will add Principle XI (transactional topology edits)
- **Feedback cycles:** Users discover gaps; amendment proposed + ratified
- **Versioning:** Constitution version bumped with each amendment

---

**Version**: 1.0 | **Ratified**: 2026-05-15 | **Last Amended**: 2026-05-15
