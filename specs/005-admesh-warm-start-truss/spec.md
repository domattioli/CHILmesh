# Feature Specification: ADMESH Warm-Start Truss Optimization

**Feature Branch**: `005-admesh-warm-start-truss` (spec lives on `planning-optimize_modernize` per CHILmesh single-branch policy)
**Created**: 2026-05-02
**Status**: Clarified (ready for `/speckit-plan`)
**Input**: User description: "feed the initial triangulation and boundary into the admesh routine so that we can run the admesh optimization/solver/whatever on that data (the truss algorithm) to achieve an optimal mesh for the underlying size function and domain without ruining the original outer and inner boundaries" + follow-up: "this functionality should be extensible. i should be able to convert any given triangulation via admesh"

## Context (informational)

CHILmesh's row 1 fixture (`chilmesh.examples.annulus()`) is deliberately poor-quality. ADMESH's main optimization engine is the distmesh truss/spring solver in `admesh.distmesh.distmesh2d`. Public entry point `admesh.triangulate(domain, h_max)` always generates fresh points from scratch — no path to pass an existing triangulation.

This spec asks for a **generic warm-start adapter**: given ANY triangulation (`points + triangles`) and domain description (SDF + size function), feed inputs into ADMESH's truss loop with boundary preserved bit-exactly. Annulus is V1 demo, but function MUST be domain-agnostic and source-agnostic — CHILmesh fixtures, `.fort.14`, `.2dm`, or raw numpy arrays.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Warm-Start Optimization Preserves Boundary (Priority: P1)

Researcher has poor-quality but topologically valid triangulation. Wants to run ADMESH's truss solver without regenerating boundary. Single function call + size function → higher-quality mesh with byte-identical boundary.

**Why this priority**: Core capability. Boundary preservation distinguishes warm-start from fresh `admesh.triangulate()`.

**Independent Test**: Call warm-start on `chilmesh.examples.annulus()` with constant size function. Verify (a) median quality strictly improves, (b) every boundary point at same `(x, y)` — bit-exact, not approximate.

**Acceptance Scenarios**:

1. **Given** the row 1 annulus mesh with median quality ≈ 0.49 and 180 boundary points, **When** the user calls the warm-start optimizer with the annulus SDF and a constant size function, **Then** the returned mesh has median quality > 0.49 AND boundary point coordinates match the input bit-exactly (np.array_equal).
2. **Given** an input mesh whose boundary is not on the SDF zero set (a malformed input), **When** the warm-start optimizer is called, **Then** it raises a clear `ValueError` naming the offending vertex indices — it MUST NOT silently project them.
3. **Given** a valid input mesh and a size function that grades from fine near boundaries to coarse in the middle, **When** the warm-start optimizer runs, **Then** elements near the boundary in the output are smaller than elements in the middle, demonstrating the size function was respected.

---

### User Story 2 - Restructured 4-Row README Visualization (Priority: P2)

README reader wants warm-start truss optimizer as Row 2 centerpiece, with FEM smoother and right-isoceles smoother applied to Row 2's result. Fresh-ADMESH-from-bbox row dropped.

**Why this priority**: User explicitly directed new layout. Visual story readers see first.

**New row layout** (replacing the existing 4-row pipeline):

| Row | Mesh | Source |
|-----|------|--------|
| 1 | Raw annulus (as-is) | `chilmesh.examples.annulus()` — unchanged from current |
| 2 | **Warm-start truss applied to Row 1** | This spec's `optimize_with_admesh_truss(row1_mesh, annulus_sdf, size_fn)` |
| 3 | FEM smoother applied to **Row 2** | `mesh_row2.smooth_mesh(method='fem', ...)` |
| 4 | Right-isoceles smoother applied to **Row 2** | `admesh.quad_prep.smooth_for_quadrangulation(mesh_row2.points, mesh_row2.connectivity_list, ...)` |

Note: rows 3 and 4 both branch off row 2 (they are alternative downstream smoothers, not sequential — row 4 is NOT row 3 + right-iso).

**Independent Test**: Render 4-row PNG. Verify (a) row 2 boundary bit-exact equal to row 1, (b) rows 3 and 4 downstream of row 2 (not row 1, not fresh ADMESH), (c) V1-V7 assertions updated to new pipeline.

**Acceptance Scenarios**:

1. **Given** the new 4-row generator script, **When** it runs, **Then** rows 2, 3, and 4 all share row 1's boundary geometry bit-exactly (boundary preservation propagates through the chain).
2. **Given** the rendered PNG, **When** displayed at GitHub README viewport width, **Then** the warm-start improvement (row 1 → row 2) is visually obvious and the FEM vs right-iso branch (row 3 vs row 4) shows distinct downstream effects.
3. **Given** the previous 4-row pipeline script (`generate_4row_admesh.py`), **When** restructured, **Then** the script either replaces the old one in place OR is renamed (e.g., `generate_4row_warmstart.py`) and the README is updated to point at the new script — design choice deferred to plan phase.

---

### User Story 3 - Generic Triangulation-In, Optimized-Mesh-Out (Priority: P1)

Downstream tool author (MADMESHR, ADMESH-Domains, notebook) wants to convert ANY triangulation through ADMESH's truss optimizer. Single function call `(points, triangles, sdf, size_fn)` → optimized triangulation regardless of source.

**Why this priority**: User explicitly stated "i should be able to convert any given triangulation via admesh." Generic-input support makes this reusable building block, not one-shot demo.

**Independent Test**: Function MUST accept inputs in **two equivalent forms**; test suite MUST exercise both:

- **Form A (CHILmesh-native)**: `optimize_with_admesh_truss(mesh: CHILmesh, sdf, size_fn, **kwargs) -> CHILmesh`
- **Form B (raw arrays)**: `optimize_with_admesh_truss_arrays(points: ndarray, triangles: ndarray, sdf, size_fn, **kwargs) -> tuple[ndarray, ndarray]`

CHILmesh form MUST be thin wrapper around raw-arrays form. Both must work on all four bundled fixtures, at least one fort.14-loaded mesh, and at least one synthetic `(points, triangles)`.

**Acceptance Scenarios**:

1. **Given** the donut fixture and a donut SDF / size pair, **When** the CHILmesh-form function is called, **Then** it returns a `CHILmesh` instance whose boundary points match the donut's input boundary bit-exactly AND median quality does not regress.
2. **Given** raw `(points, triangles)` numpy arrays for a square domain (constructed in the test, not loaded from a fixture), **When** the raw-array form is called with a square SDF, **Then** it returns optimized `(points, triangles)` arrays where the four corner points and the boundary chain are bit-exact preserved.
3. **Given** a CHILmesh loaded from a `.fort.14` file, **When** the CHILmesh form is called, **Then** the output preserves the loaded boundary and the function can roundtrip via `write_fort14` without error.
4. **Given** an unsupported input (e.g., a quad mesh), **When** the function is called, **Then** it raises a clear `NotImplementedError` referring to "triangle-only for v1".
5. **Given** a triangulation whose boundary is NOT consistent with the supplied SDF (i.e., user passed the wrong domain), **When** the function is called, **Then** it raises `ValueError` listing offending boundary indices, naming the SDF tolerance used, and surfacing the worst-case absolute SDF value seen — i.e., enough information to diagnose without re-running.

---

### Edge Cases

- **Boundary identification**: `mesh.boundary_edges()` returns edge IDs, not vertex indices. Function MUST use deterministic boundary-detection rule and document it.
- **SDF mismatch**: Boundary not on SDF zero set → raise, not silently drag boundary off geometry.
- **Connectivity vs re-triangulation**: ADMESH re-triangulates every iteration. Output connectivity ≠ input. "Warm-start" means input *points* as initial distribution, NOT preserving input *triangles*.
- **Convergence failure**: Hit max iterations → surface RuntimeWarning, return best-so-far mesh.
- **Mixed/quad inputs**: V1 triangles only. Quad/mixed → `NotImplementedError`.
- **Domain-boundary mismatch**: SDF doesn't contain input boundary points → raise.
- **Degenerate triangles**: Zero or negative area → raise before feeding to ADMESH.
- **Tiny meshes**: <10 interior points → still run, issue RuntimeWarning.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The system MUST expose **two layered Python entry points**:
  - **FR-001a (low-level / generic)**: `optimize_with_admesh_truss_arrays(points: ndarray[N,2|3], triangles: ndarray[M,3], sdf: Callable, size_fn: Callable, **kwargs) -> tuple[ndarray, ndarray]`. This is the canonical implementation; it accepts ANY triangulation regardless of source.
  - **FR-001b (high-level / CHILmesh convenience)**: `optimize_with_admesh_truss(mesh: CHILmesh, sdf: Callable, size_fn: Callable, **kwargs) -> CHILmesh`. This is a thin wrapper that extracts `(points, triangles, boundary_indices)` from the CHILmesh, calls FR-001a, and rewraps the result.
  - Both forms MUST share validation logic and produce equivalent results given equivalent inputs.
- **FR-002**: For the CHILmesh form (FR-001b), boundary vertices MUST be identified deterministically from the input CHILmesh's adjacency structures (using `boundary_edges()` → `Edge2Vert` lookup → `np.unique` on flattened vertex indices). For the raw-arrays form (FR-001a), the caller MAY pass an explicit `boundary_indices` argument; if omitted, boundary indices MUST be inferred by detecting edges that appear in exactly one triangle.
- **FR-003**: The function MUST pin all identified boundary vertices via ADMESH's `pfix` mechanism so the truss loop never moves them.
- **FR-004**: The function MUST inject the input mesh's interior points (non-boundary vertices) as the initial point distribution for the truss loop, bypassing ADMESH's `_initial_distribution` rejection sampler.
- **FR-005**: The function MUST validate, before invoking the truss loop, that every input boundary vertex satisfies `|sdf(v)| < 1e-6` (a documented tolerance) — and raise `ValueError` listing offending indices if not.
- **FR-006**: The function MUST validate that the input is triangle-only (no quads, no mixed) and raise `NotImplementedError` otherwise.
- **FR-007**: The function MUST validate that all input triangles have positive signed area and raise `ValueError` listing offending element indices if not.
- **FR-008**: The function MUST return a CHILmesh whose `points` array contains every input boundary point at the same index AND coordinate as the input. Bit-exact equality on the boundary subset is REQUIRED.
- **FR-009**: The function MUST forward optional truss-solver parameters (max iterations, displacement tolerance, retriangulation cadence, force-scaling factor) to ADMESH's underlying `distmesh2d` call with documented defaults that match how row 3 of the existing 4-row visualization invokes ADMESH.
- **FR-010**: The function MUST NOT modify the input CHILmesh in place. The output is always a fresh instance.
- **FR-011**: The function MUST surface a `RuntimeWarning` (not raise) if the truss loop hits the max-iteration cap without converging, and still return the best-so-far mesh. **Additionally, if the post-loop median quality is *strictly lower* than the input median quality (a regression), the function MUST return the *input* mesh unchanged and emit a `RuntimeWarning` stating "warm-start truss did not improve quality (input median=X, output median=Y); returning input." This is the non-degradation guarantee — no caller ever gets a worse mesh from this function.** [Resolved Q4=b: graceful degradation, not raise.]
- **FR-012**: The existing `generate_4row_admesh.py` script MUST be restructured to the new layout (US2 table): Row 1 = raw, Row 2 = warm-start applied to Row 1, Row 3 = FEM applied to Row 2, Row 4 = right-isoceles applied to Row 2. The fresh-ADMESH-from-bbox row of the current pipeline is dropped. The script's filename MAY be renamed to reflect the new pipeline (e.g., `generate_4row_warmstart.py`); the README MUST be updated to point at whatever the new authoritative script is. **The output PNG path remains `tests/output/annulus_quickstart.png` so README image links don't break.** [Resolved Q2=d: restructure existing rows around warm-start.]
- **FR-013**: The demo script MUST enforce fail-loud assertions before saving any output:
  - **V_BND**: Row 2's boundary points are bit-exact equal to Row 1's boundary points (warm-start preservation).
  - **V_BND_PROP**: Rows 3 and 4 boundary points are bit-exact equal to Row 2's boundary points (preservation propagates through downstream smoothers — to the extent each smoother claims to preserve boundary; if a smoother does not preserve boundary by design, this assertion is relaxed to "boundary unchanged within smoother's documented tolerance").
  - **V_QI**: Row 2 median quality is strictly greater than Row 1 median quality (warm-start actually helps in this demo; if the warm-start regresses on the annulus, the script fails — separate from the function-level non-degradation guarantee in FR-011, which silently returns input on regression).
  - **V_CONN**: All four rows define valid triangulations (every triangle has positive area, no duplicate vertex indices).
  - **V_CHAIN**: Row 3's input was Row 2's `(points, triangles)` (verified by reference identity or by computing input hash). Row 4's input was Row 2's `(points, triangles)`. Neither rows 3 nor 4 secretly used Row 1.
- **FR-014**: The function and demo script MUST be importable / runnable without ADMESH installed system-wide; the existing project convention of `sys.path.insert(0, str(Path(__file__).parent / "ADMESH"))` (or equivalent) is acceptable. If ADMESH cannot be imported, the function MUST raise a clear `ImportError` naming the missing dependency.
- **FR-015**: The integration MUST pin the ADMESH dependency to a **specific known-good commit SHA** (not `main`) and document that SHA at the top of the adapter file. The plan phase MUST identify a commit whose `admesh/distmesh.py` and `admesh/routine.py` import each other cleanly (i.e., before the `MeshOutput` regression). The pinned SHA is recorded in a single source-of-truth location (e.g., a constant in the adapter module or a `requirements.txt`-style file) and referenced everywhere else. Upgrading the pin in the future is a single-place change, not a sweep. [Resolved Q1=a: pin, do not upstream-fix or vendor.]
- **FR-016**: The function MUST be **input-source-agnostic**. It MUST work on:
  - Any of the four bundled CHILmesh fixtures (`annulus`, `donut`, `block_o`, `structured`).
  - A CHILmesh loaded via `read_from_fort14()`.
  - A CHILmesh loaded via `read_from_2dm()`.
  - A CHILmesh loaded via `from_admesh_domain()`.
  - Raw `(points, triangles)` numpy arrays produced ad-hoc (e.g., from `scipy.spatial.Delaunay`).
  Any user-facing message that names "CHILmesh" specifically MUST be in the high-level wrapper only; the low-level function MUST refer to "triangulation" generically.
- **FR-017**: The function MUST be **domain-agnostic**. The SDF and size-function arguments are pure callables; the function MUST NOT hard-code annulus parameters, ring radii, or any other domain-specific assumption. The annulus is the v1 demo, not a special case in the library code.
- **FR-018**: A library-level docstring (and a sibling `usage.md` or equivalent) MUST show at least three worked examples covering: (a) bundled annulus, (b) bundled donut, (c) raw-arrays from a non-CHILmesh source. This documents the extensibility contract for downstream users.

### Key Entities

- **WarmStartInput**: `(CHILmesh, SDF, SizeFn)` + optional kwargs. Optimizer input contract.
- **BoundaryVertexSet**: Vertex indices on domain boundary, identified via `boundary_edges()` → `Edge2Vert`. Becomes `pfix` for ADMESH.
- **InteriorVertexSet**: All vertex indices NOT in BoundaryVertexSet. Free nodes truss loop moves.
- **TrussSolverConfig**: Optional kwargs forwarded to `distmesh2d` (max iter, tolerance, retriangulation cadence, force scaling). Documented defaults.
- **WarmStartOutput**: Fresh `CHILmesh` from optimized point set + ADMESH's final retriangulation. Boundary point indices preserved.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: On the bundled annulus fixture (`chilmesh.examples.annulus()`) with a constant size function, the warm-start optimizer produces a mesh whose median element quality is at least **0.60** (compared to the input's ≈ 0.49 — a ≥22% relative improvement). This is the demo benchmark; the function-level requirement is non-degradation (FR-011 graceful-degradation guarantee). [Resolved Q5=b: 0.60 target, clearly better than the ad-hoc distmesh prototype's 0.5617.]
- **SC-002**: On the same fixture, **100 %** of input boundary points appear in the output at bit-exact equal `(x, y)` coordinates (np.array_equal == True). This is the boundary-preservation guarantee.
- **SC-003**: The warm-start optimizer completes on a 580-element mesh (annulus) in under **30 seconds** of wall time on a developer laptop without ADMESH-side parallelization.
- **SC-004**: The demo visualization renders without invoking any of `_initial_distribution`, the rejection-sampling code path, or any ADMESH internals that generate fresh points — verified via assertion or instrumentation, not vibes.
- **SC-005**: All four bundled CHILmesh fixtures (annulus, donut, block_o, structured) can be passed to the warm-start optimizer without crashing. Quality may not strictly improve on every fixture (some are already excellent), but the function MUST NOT regress median quality below input on any fixture.
- **SC-006**: The README (or a sibling docs page) clearly explains when to use warm-start vs `admesh.triangulate()`-from-scratch, in plain language, with one or two sentences of guidance.
- **SC-007 (extensibility)**: The test suite MUST include at least three input-source variations (bundled fixture, fort.14-loaded mesh, raw numpy arrays from a non-CHILmesh source) and verify all three produce valid optimized output with bit-exact boundary preservation. If any of these three input forms fails, the feature is not done.
- **SC-008 (extensibility)**: The test suite MUST include at least two distinct domains: **annulus + donut** (the CHILmesh `donut` fixture, since it's already in the test corpus, has a non-trivial geometry distinct from annulus, and exercises the same SDF→pfix→truss pipeline). [Resolved Q3=d: donut as second domain.] The test for the donut MUST verify (a) bit-exact boundary preservation, (b) no quality regression, (c) the donut SDF is supplied by the test (not hard-coded inside the library) — this is the proof that the library is domain-agnostic.

## Cross-Repo Tracking Policy

Any extensibility hook, API addition, or upstream fix belonging inside ADMESH MUST be filed on `https://github.com/domattioli/ADMESH/issues` — not implemented as CHILmesh workaround. CHILmesh adapter is thin shim only.

Candidate ADMESH issues (file by user — Claude's MCP scope restricted to `domattioli/chilmesh`):

- **ADMESH-A**: Broken import — `MeshOutput` is imported by `admesh/routine.py` but not defined in `admesh/distmesh.py` on `main` HEAD. Fix: define and export `MeshOutput` (or remove the import). This blocks any external caller of `admesh.triangulate()`.
- **ADMESH-B**: Public warm-start entry point — `admesh.triangulate()` currently always runs `_initial_distribution` to generate fresh points. Proposal: add an optional `initial_points: ndarray | None` argument (or a sibling function `triangulate_from_points(...)`) that bypasses rejection sampling and uses the caller's points as the initial truss state.
- **ADMESH-C**: Public `pfix` documentation — `Domain(..., pfix=...)` is used here as the boundary-pinning mechanism. Confirm that `pfix` points are guaranteed bit-exact preserved by the truss loop (no floating-point drift), and document this contract.
- **ADMESH-D**: Distmesh re-triangulation cadence — document the current behavior (Delaunay re-triangulation every iteration) and expose the cadence as a tunable. Some warm-start callers may want to skip re-triangulation for the first N iters to keep the input topology longer.
- **ADMESH-E**: Convergence diagnostics — surface the per-iteration max-displacement and quality-stat history as a return value from `distmesh2d`, so callers can produce convergence plots and decide whether the loop converged sensibly.

CHILmesh adapter MUST work around unfixed upstream issues. Each workaround MUST link to tracking issue and be marked TEMP in code.

## Assumptions

- ADMESH source reachable at `https://github.com/domattioli/ADMESH.git`; team can clone known-good commit if `main` broken.
- Existing `sys.path.insert(0, "ADMESH")` convention continues to work.
- Annulus SDF accepts row 1 boundary points within tolerance (verified during spec→plan transition).
- "Truss algorithm" = ADMESH's distmesh-style spring/repulsion solver (`distmesh2d`), not generic structural-mechanics solver.
- Boundary preservation hard constraint: bit-exact equality via `pfix` mechanism.
- V1 scope: triangles only. Quad/mixed warm-start explicitly deferred.
- CHILmesh single-branch policy applies; `specs/005-admesh-warm-start-truss/` added without creating separate git branch.
- Warm-start row is *additive*, not replacement for any existing pipeline row.
