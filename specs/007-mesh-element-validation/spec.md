# Feature Specification: Mesh Element Validity Test Suite (Quad-Dominant, Boundary-Triangle, Planar, Non-Self-Intersecting)

**Feature Branch**: `claude/mesh-quad-triangle-spec-IIvKw`
**Created**: 2026-05-21
**Status**: Draft
**Input**: User description: "Test suite to verify that a mesh is comprised completely of quads, except for triangles on the boundary, all of which are planar 2D non-intersecting elements. Degenerate quads are acceptable but not self-intersecting ones."

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Mesh consumer trusts element-type composition (Priority: P1)

A downstream consumer (MADMESHR, ADMESH, ADMESH-Domains) loads a CHILmesh and assumes the element population is **quad-dominant with triangles confined to the boundary layer**. The consumer should not have to defensively guard against interior triangles, pentagons, or non-conforming polygons.

**Why this priority**: This is the single load-bearing assumption every downstream pipeline (FEM smoother, skeletonization rings, layer paths) makes about CHILmesh-produced meshes. A violation cascades silently into every analysis. Without an explicit test gate this assumption lives only in maintainer memory.

**Independent Test**: For each of the four built-in fixtures (`annulus`, `donut`, `block_o`, `structured`), invoke a single suite entry point — `validate_mesh_elements(mesh)` — and assert it returns `ok=True` with zero violations across the element-type rules. The test produces a structured report of any offending element IDs so a maintainer can localize the failure without re-running.

**Acceptance Scenarios**:

1. **Given** any built-in fixture, **When** the suite enumerates all elements, **Then** every element has either 3 or 4 distinct (or padded) vertex slots — no element has 5+ vertices and no element has fewer than 3 distinct topological corners.
2. **Given** any built-in fixture, **When** the suite classifies each element by type, **Then** every triangle satisfies the boundary-triangle rule (see FR-006) and every quad is a topological quadrilateral.
3. **Given** a synthetic mesh fixture that injects an interior triangle (a triangle whose vertices are all interior, i.e., none lie on `boundary_vertices()`), **When** the suite runs, **Then** it FAILS with the violation `INTERIOR_TRIANGLE_FORBIDDEN` and reports the offending element ID(s).
4. **Given** a synthetic fixture containing a pentagon (5-vertex element), **When** the suite runs, **Then** it FAILS with `UNSUPPORTED_ELEMENT_ARITY` and reports the offending element ID(s).

---

### User Story 2 - Mesh consumer trusts geometric validity (Priority: P1)

A downstream consumer assumes every element is a **planar 2D polygon** that is **not self-intersecting** (no bowtie quads) and that **distinct elements do not overlap** in the interior (they touch only at shared edges or shared vertices).

**Why this priority**: Self-intersecting or overlapping elements silently corrupt every quadrature rule, every signed-area computation, every smoother. Floating-point area sign flips into negative without warning. Skeletonization-layer logic produces nonsense. This must be a regression gate.

**Independent Test**: Run `validate_mesh_elements(mesh)` on each built-in fixture and assert no violations of category `SELF_INTERSECTING`, `INTERIOR_OVERLAP`, or `EDGE_CROSSING`. The suite must also detect planted bowtie quads in a synthetic fixture.

**Acceptance Scenarios**:

1. **Given** any built-in fixture, **When** the suite checks every quad for self-intersection (diagonals-cross-on-both-pairs test or equivalent), **Then** zero violations are reported.
2. **Given** any built-in fixture, **When** the suite computes pairwise edge crossings restricted to non-adjacent elements (broadphase via edge AABB / spatial index, narrowphase via segment-segment crossing with shared-endpoint exemption), **Then** zero crossings are reported.
3. **Given** any built-in fixture, **When** the suite verifies that distinct element interiors are disjoint, **Then** for every pair of elements that share at most one vertex (no shared edge), the polygon interiors do not overlap.
4. **Given** a synthetic fixture containing a bowtie quad (e.g., vertices `[0,0], [1,1], [1,0], [0,1]` in that connectivity order), **When** the suite runs, **Then** it FAILS with `SELF_INTERSECTING_QUAD` and reports the element ID and the two crossing edges.
5. **Given** a synthetic fixture containing two quads whose interiors partially overlap (without sharing an edge), **When** the suite runs, **Then** it FAILS with `INTERIOR_OVERLAP` and reports the element-pair.

---

### User Story 3 - Mesh maintainer keeps degenerate quads as a first-class accepted case (Priority: P1)

A maintainer knows CHILmesh accepts **degenerate quads** — quads with collinear vertices, zero or near-zero area, or a padded-triangle-as-quad encoding — as a legitimate input. The validity suite must distinguish *degenerate* (allowed) from *self-intersecting* (forbidden) and never conflate the two.

**Why this priority**: CHILmesh's mixed-element story (triangles encoded as 4-column padded quads with `-1` sentinel, or as quads with one duplicated vertex) is already tested in `tests/test_degeneracy.py`. The new suite must not regress those flows by tightening the gate too aggressively.

**Independent Test**: Run the new suite against the existing `tests/test_degeneracy.py` scenarios (padded-triangle quad, duplicate-vertex quad, mixed-element mesh). Suite must report `ok=True`. Then plant a true bowtie quad in the same scenarios; suite must fail.

**Acceptance Scenarios**:

1. **Given** a quad with one repeated vertex index (degenerate triangle-as-quad), **When** the suite runs, **Then** it PASSES (degenerate ≠ invalid).
2. **Given** a quad with four collinear vertices (zero signed area, no self-crossing), **When** the suite runs, **Then** it PASSES but emits a `DEGENERATE_ZERO_AREA` informational record (not a failure) so downstream code can find these elements without running the suite again.
3. **Given** the padded-triangle-as-quad fixture in `tests/test_degeneracy.py`, **When** the suite runs alongside the existing degeneracy test, **Then** both tests pass.
4. **Given** a quad whose vertex ordering causes its edges to cross (bowtie), **When** the suite runs, **Then** it FAILS even though the signed area may be near zero (a bowtie can have small signed area; the suite must not use area as a proxy for self-intersection).

---

### User Story 4 - Mesh maintainer can run the suite as a regression gate (Priority: P2)

A maintainer adds the new suite to `pytest -v` and to CI. It runs in seconds on `annulus`, `donut`, `structured`, and in under a minute on `block_o`. Output on failure points at the offending element ID and violation category — no manual mesh inspection required.

**Why this priority**: A test suite that runs too slowly gets disabled. A test suite with unclear failure messages gets ignored. Both kill the gate.

**Independent Test**: Time the suite per fixture; assert wall-clock budgets per fixture. Inject a single planted violation; assert the failure message contains the element ID, the violation category, and at least one geometric specific (edge IDs / vertex coords).

**Acceptance Scenarios**:

1. **Given** the four built-in fixtures, **When** the suite runs under `pytest -v`, **Then** total runtime is under 60 seconds on a developer laptop (block_o ≤ 45 s, others ≤ 5 s each).
2. **Given** a single planted bowtie quad in a fixture, **When** the suite runs, **Then** the failure message names the element ID, the violation category (`SELF_INTERSECTING_QUAD`), the two crossing edges, and at minimum the four vertex IDs.
3. **Given** the suite passes locally, **When** the CI job runs the same suite on the same fixtures, **Then** results match (no flakes from non-determinism in pair enumeration / floating-point reductions).

---

### Edge Cases

- **Triangle-as-padded-quad encoding**: `connectivity_list` row `[a, b, c, -1]` (or `[a, b, c, c]`). Classified as triangle for type-composition rule; classified as degenerate quad for geometric rules; both classifications must pass.
- **Coincident vertices**: A quad with two identical vertex *coordinates* but distinct vertex *IDs* (encoded as a sliver). Degenerate, allowed, no self-intersection — must pass.
- **Numerical-precision near-bowtie**: A quad that is geometrically *almost* bowtie but technically not (edges share endpoint exactly, no proper crossing). Treated as non-self-intersecting; pass. Tolerance is explicit (see FR-013).
- **Floating-point coincident edges**: Two quads whose shared edge endpoints match to within tolerance but not exact. Treated as adjacent for the edge-crossing test (broadphase exemption).
- **Boundary triangle with all vertices on boundary**: Allowed.
- **Boundary triangle with exactly one vertex on the boundary** (other two interior): Allowed per user clarification — "must have at least one vert on the boundary."
- **Interior triangle (zero boundary vertices)**: Forbidden — `INTERIOR_TRIANGLE_FORBIDDEN`.
- **Disconnected components**: Suite runs per-component; violations are global element IDs.
- **Open-ocean / non-closed boundaries**: Treated like any other boundary edge for the triangle-membership rule (a vertex on an open-ocean boundary still counts as a boundary vertex).
- **Empty mesh** (`n_elems == 0`): Suite passes trivially; emits informational `EMPTY_MESH`.
- **Single-element mesh**: Suite checks only the element-validity rules (no pairs to compare).

## Requirements *(mandatory)*

### Functional Requirements

#### Suite entry point

- **FR-001**: A single function `validate_mesh_elements(mesh: CHILmesh, *, tol: float | None = None) -> MeshValidityReport` MUST exist as a **test-suite-internal helper** under `tests/_validity/` (NOT a public `chilmesh.*` module). It MUST NOT mutate the mesh. Long-term promotion to `chilmesh.validate` is tracked in the discussion issue referenced in **Clarifications Q6**; the spec deliberately does NOT introduce a new public API surface at this stage. **(Updated per Clarify Q6.)**
- **FR-002**: `MeshValidityReport` MUST expose at minimum: `ok: bool`, `violations: list[Violation]`, `notes: list[InformationalNote]`, `n_elems_checked: int`, `runtime_s: float`. Each `Violation` MUST include `category: str`, `element_ids: tuple[int, ...]`, `edge_ids: tuple[int, ...] | None`, `detail: str`.

#### Element-type composition

- **FR-003**: Every element MUST classify as exactly one of: `TRI` (3 distinct vertex IDs, fourth slot `-1` or repeats one of the first three), `QUAD` (4 distinct vertex IDs). Any other arity raises `UNSUPPORTED_ELEMENT_ARITY`.
- **FR-004**: For each `QUAD`, the four vertex IDs MUST be pairwise distinct. A quad with a duplicate-but-non-padding vertex (i.e., not in the padded-triangle encoding) is `DEGENERATE_QUAD_DUPLICATE_VERTEX` — recorded as an informational note, NOT a violation (matches existing `test_degeneracy.py` behavior).
- **FR-005**: For each `TRI`, the three distinct vertex IDs MUST be pairwise distinct.
- **FR-006**: Every `TRI` MUST have at least one vertex that is a member of `mesh.boundary_vertices()`. A triangle with zero boundary vertices raises `INTERIOR_TRIANGLE_FORBIDDEN`.
- **FR-007**: Cross-check against skeletonization: every `TRI` MUST belong to **layer 0** (the outermost layer) per `mesh.layers["OE"][0]`. A triangle in `layers["OE"][k]` or `layers["IE"][k]` for any `k ≥ 1` raises `INTERIOR_LAYER_TRIANGLE_FORBIDDEN`. (Per user clarification: "by definition of the chilmesh layers triangles may only exist in the boundary layer.") If `mesh.layers` has not yet been computed, the validator MUST trigger `_skeletonize()` once and rely on the existing layer cache on `mesh` (no re-computation per call). **(Clarified per Clarify Q8.)**

#### Planarity (2D)

- **FR-008**: All vertices MUST have exactly 2 coordinate components (`mesh.points.shape[1] == 2`), or, if 3D coordinates are present, all `z`-values MUST be identical within `tol` (i.e., the mesh lies in a single plane parallel to the xy-plane). A mesh with non-coplanar `z`-values raises `NON_PLANAR_MESH` once at the report level (not per-element).

#### Self-intersection (single element)

- **FR-009**: Every `QUAD` MUST NOT be self-intersecting. A quad `v0, v1, v2, v3` is self-intersecting iff edge `(v0,v1)` properly crosses edge `(v2,v3)`, OR edge `(v1,v2)` properly crosses edge `(v3,v0)`. "Properly crosses" means the two open segments intersect at an interior point of both. Violation: `SELF_INTERSECTING_QUAD`.
- **FR-010**: Every `TRI` is trivially non-self-intersecting; no test required (a triangle with 3 distinct, non-collinear vertices cannot self-intersect; the collinear case is degenerate-but-not-self-intersecting, allowed).

#### Element-element non-overlap

- **FR-011**: For every pair of elements `(e_i, e_j)` that share at most one vertex (i.e., do NOT share an edge), the polygon interiors MUST be disjoint. Interior overlap is detected by checking: (a) any vertex of `e_i` lies strictly inside `e_j` (or vice versa), OR (b) any edge of `e_i` properly crosses any edge of `e_j`. Violation: `INTERIOR_OVERLAP` (when (a) triggers) or `EDGE_CROSSING` (when (b) triggers).
- **FR-012**: Pair enumeration MUST use a broadphase (axis-aligned bounding box overlap, accelerated by Phase-5 spatial index if available, else a uniform grid hash) so worst-case is `O(n log n)` for typical inputs. A literal `O(n²)` pair enumeration is acceptable for `n_elems ≤ 5000` but MUST NOT be the default code path on `block_o`-class meshes.

#### Tolerances & numerics

- **FR-013**: All geometric predicates MUST accept an explicit `tol` parameter. **Default behavior is coordinate-system-agnostic**: when the caller passes `tol=None` (the default), the validator computes `bbox_diag = ||(max_x - min_x, max_y - min_y)||_2` and uses `tol_effective = 1e-12 * bbox_diag` (with a floor of `1e-15` to avoid zero on degenerate empty bboxes). Callers may override with an absolute `tol` if they need stricter or coord-specific behavior. The point-in-polygon test MUST use a robust orientation predicate (signed-area sign with `tol_effective`-aware tie-breaking). Shared-endpoint detection MUST treat two endpoints as identical when their L2 distance is `≤ tol_effective`. Tolerances and the bbox-relative formula MUST be documented in the function docstring. **(Clarified per Clarify Q7; aligns with constitution principle V — coordinate-system agnosticism.)**
- **FR-014**: Degenerate quads (collinear vertices, near-zero signed area) MUST NOT be reported as violations. They MAY be reported as `DEGENERATE_ZERO_AREA` informational notes for downstream consumers.

#### Synthetic fixtures (for negative tests)

- **FR-015**: The test suite MUST ship synthetic fixtures, each constructed in-test, that plant exactly one violation per fixture: `bowtie_quad_mesh`, `interior_triangle_mesh`, `pentagon_mesh`, `overlapping_quads_mesh`, `edge_crossing_mesh`. Each fixture MUST be small (≤ 20 elements) so failure messages are easy to read. Because the `CHILmesh` constructor performs its own degeneracy fallback (see `tests/test_degeneracy.py`) which may silently repair some of the planted violations, fixtures that depend on bypassing constructor validation MUST build the `CHILmesh` instance, then **post-mutate** `connectivity_list` (and any cached adjacencies) directly to inject the violation. The mutation helper lives next to the fixtures in `tests/_validity/` and is not part of the public API.

#### Reporting

- **FR-016**: On test failure, the pytest assertion message MUST list at most the first 10 violations per category (to bound output), with each entry containing element IDs, edge IDs (where applicable), and a one-line `detail` string. The full violation list MUST be accessible programmatically via the returned `MeshValidityReport` for callers invoking `validate_mesh_elements` directly.
- **FR-017**: Pytest test shape — **one test per fixture** that asserts `report.ok`. The failure message MUST aggregate all violation categories triggered on that fixture (not just the first) so a maintainer sees the full picture from a single failing test. Parametrization is over the four built-in fixtures plus the five synthetic negative fixtures from FR-015. **(Clarified per Clarify Q9.)**

### Key Entities

- **`MeshValidityReport`**: Aggregated output of one validation pass. Fields: `ok`, `violations`, `notes`, `n_elems_checked`, `runtime_s`. Hashable summary, JSON-serializable.
- **`Violation`**: A single rule failure. Fields: `category` (one of the FR-defined error codes), `element_ids` (tuple), `edge_ids` (tuple or None), `detail` (free-form string with vertex coords / geometric specifics).
- **`InformationalNote`**: A non-failure observation worth recording (e.g., `DEGENERATE_ZERO_AREA`, `EMPTY_MESH`). Same shape as `Violation` but does not flip `ok` to False.

### Violation Categories (canonical list)

| Code | Severity | FR |
|------|----------|------|
| `UNSUPPORTED_ELEMENT_ARITY` | violation | FR-003 |
| `DEGENERATE_QUAD_DUPLICATE_VERTEX` | note | FR-004 |
| `INTERIOR_TRIANGLE_FORBIDDEN` | violation | FR-006 |
| `INTERIOR_LAYER_TRIANGLE_FORBIDDEN` | violation | FR-007 |
| `NON_PLANAR_MESH` | violation | FR-008 |
| `SELF_INTERSECTING_QUAD` | violation | FR-009 |
| `INTERIOR_OVERLAP` | violation | FR-011 |
| `EDGE_CROSSING` | violation | FR-011 |
| `DEGENERATE_ZERO_AREA` | note | FR-014 |
| `EMPTY_MESH` | note | edge case |

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: All four built-in fixtures (`annulus`, `donut`, `block_o`, `structured`) pass `validate_mesh_elements(mesh).ok == True` with zero violations. Informational notes are permitted (and expected on fixtures known to contain degenerate-quad triangle encodings).
- **SC-002**: All five synthetic negative fixtures (FR-015) FAIL the suite with the expected violation category and at least the offending element ID(s) present in the report.
- **SC-003**: The new suite is added as `tests/test_mesh_element_validity.py` (with helpers in `tests/_validity/`) and runs as part of the default `pytest -v` invocation. CI passes. **No new files under `src/chilmesh/`** (per Clarify Q6 — no public API change at this stage).
- **SC-004**: Wall-clock budget: total runtime across all four fixtures is ≤ 60 s on a developer laptop. Block_o alone is ≤ 45 s; each of the other three is ≤ 5 s. Reported in pytest output via `pytest --durations=10`.
- **SC-005**: The existing 288-test suite (and in particular `tests/test_degeneracy.py`, `tests/test_invariants.py`, `tests/test_signed_area.py`, `tests/test_skeletonization_invariant.py`) continues to pass without modification.
- **SC-006**: A planted bowtie quad injected into the `annulus` fixture causes the suite to fail with `SELF_INTERSECTING_QUAD` and the failure message names the offending element ID, the two crossing edges, and the four vertex coordinates. Verified by a parametrized "negative path" test inside `test_mesh_element_validity.py`.
- **SC-007**: On `block_o` (≈5,200 elements), pair-enumeration broadphase reduces narrowphase calls by at least 100× compared to literal `O(n²)`. Measured via an internal counter exposed only in `DEBUG` builds or via a `--validate-stats` opt-in flag; not asserted in regular CI but documented in `plan.md`.
- **SC-008**: A GitHub discussion issue is opened on `domattioli/CHILmesh` referencing this spec and asking whether the validator should be promoted from `tests/_validity/` into a public `chilmesh.validate` module after the QuADMesh dev cycle settles. Tracked at **#142** (labels: `enhancement`, `priority:low`, `status:voting`, `type:decision`). No README / CHANGELOG entry at this stage (no public API change). **(Updated per Clarify Q6.)**

## Clarifications

**Q1** (user, 2026-05-21): "non planar 2d not intersecting elements" → User clarified: **planar 2D**, with all three intersection checks active: (a) no self-intersecting element, (b) no element-element interior overlap, (c) no edge-edge crossings between non-adjacent elements.

**Q2** (user, 2026-05-21): Boundary triangle definition → "By definition of the chilmesh layers triangles may only exist in the boundary layer so they must have at least one vert on the boundary." Encoded as FR-006 (≥1 boundary vertex) + FR-007 (must reside in `layers["OE"][0]`).

**Q3** (resolved by spec): Degenerate vs self-intersecting → Degenerate (zero area, collinear, padded-triangle, duplicate vertex) is **allowed and reported as a note, not a violation**. Self-intersecting (bowtie connectivity) is **forbidden** regardless of signed area. The two are detected by separate predicates; signed area is NOT used as a proxy for self-intersection.

**Q4** (resolved by spec): Pair-enumeration scaling → Broadphase mandated for meshes above ~5,000 elements; `O(n²)` acceptable below that threshold to keep the implementation simple on small fixtures.

**Q5** (resolved by spec): Tolerance policy → Single explicit `tol` parameter on the entry point, with documented per-predicate behavior. No silent floating-point thresholds. *(Superseded by Q7 — default is now bbox-relative.)*

---

### Clarify pass (2026-05-21)

**Q6 — Module placement / public-API status**: This is a **short-term test-suite-only helper** to guide QuADMesh development; long-term promotion to `chilmesh.validate` is undecided. **Resolution:** Live under `tests/_validity/` (NOT `src/chilmesh/`). No public API change at this stage. A discussion issue is open at **#142** (labels: `enhancement`, `priority:low`, `status:voting`, `type:decision`) for the maintainer + downstream consumers to vote on promotion.

**Q7 — Tolerance default vs coord-system-agnostic constitution**: Default `tol` is now **relative to mesh bbox diagonal** (`1e-12 * bbox_diag`, floored at `1e-15`). Callers may override with an absolute value. Updates FR-013 and aligns with constitution principle V.

**Q8 — FR-007 layer cross-check when skeletonization not yet computed**: Validator **auto-triggers `_skeletonize()` once** and relies on the existing `mesh.layers` cache for subsequent calls. No per-call recomputation. Updates FR-007.

**Q9 — Pytest UX**: **One test per fixture** asserting `report.ok`; failure message aggregates all violation categories triggered. Parametrized over the four built-in fixtures plus the five synthetic negative fixtures from FR-015. Adds FR-017.

## Assumptions

- CHILmesh exposes `mesh.boundary_vertices()`, `mesh.boundary_edges()`, `mesh.layers`, `mesh.connectivity_list`, `mesh.points` — all of which exist in the v0.2.0 public API surface.
- Triangle encoding in `connectivity_list` uses the padded form documented in `.claude/CLAUDE.md` (4 columns, `-1` sentinel for the missing fourth vertex on triangles). The suite handles `-1` padding explicitly per FR-003.
- The fort.14 reader and 2DM reader produce element orderings consistent with `Elem2Vert` (CCW for positive signed area). Negative signed area indicates either reversed winding (handled by `signed_area()`) or bowtie self-intersection (caught by FR-009 independently of winding).
- Spatial index from Phase 5 may or may not be available at suite-write time. The implementation uses a uniform-grid hash as the default broadphase and opportunistically uses the Phase-5 index when present.
- Open-ocean boundaries are treated as ordinary boundary edges for the purpose of FR-006 (a vertex on an open-ocean boundary edge counts as a boundary vertex). This matches the existing `mesh.boundary_vertices()` semantics.
- Visualization is out of scope. The suite emits a structured report only; rendering of violation locations on a matplotlib axis is deferred to a follow-up.
- The suite is read-only: no mesh repair, no auto-fix, no mutation. Failures point at the data; humans (or a separate skill) decide remediation.
