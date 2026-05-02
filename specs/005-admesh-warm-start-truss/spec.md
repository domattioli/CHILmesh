# Feature Specification: ADMESH Warm-Start Truss Optimization

**Feature Branch**: `005-admesh-warm-start-truss` (spec lives on `planning-optimize_modernize` per CHILmesh single-branch policy)
**Created**: 2026-05-02
**Status**: Draft
**Input**: User description: "feed the initial triangulation and boundary into the admesh routine so that we can run the admesh optimization/solver/whatever on that data (the truss algorithm) to achieve an optimal mesh for the underlying size function and domain without ruining the original outer and inner boundaries" + follow-up: "this functionality should be extensible. i should be able to convert any given triangulation via admesh"

## Context (informational)

CHILmesh's row 1 fixture (`chilmesh.examples.annulus()`) is a deliberately poor-quality triangulation — a fixed annular boundary with random-Delaunay interior. ADMESH's main optimization engine is the distmesh-style truss/spring solver in `admesh.distmesh.distmesh2d` (referenced at `https://github.com/domattioli/ADMESH/blob/05bc68fc81060f7d710b8f4abb2cc382f85df33f/admesh/routine.py#L27`). Today the public entry point `admesh.triangulate(domain, h_max)` always **generates** an initial point distribution from scratch (rejection sampling within the domain), then runs the truss loop. There is no path to pass in an existing triangulation as the starting state.

This spec asks for a **generic warm-start adapter**: given ANY triangulation (`points + triangles`) and a description of the domain it lives in (SDF + size function), feed those inputs into ADMESH's truss loop so the interior nodes are optimized against the size function while the input boundary is preserved bit-exactly. The annulus is the V1 demo, but the function MUST be domain-agnostic and triangulation-source-agnostic — it must work on the four bundled CHILmesh fixtures, on a CHILmesh built from a `.fort.14` file, on a CHILmesh built from a `.2dm` file, and on a raw `(points, triangles)` numpy pair without a CHILmesh wrapper.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Warm-Start Optimization Preserves Boundary (Priority: P1)

As a researcher studying mesh-quality improvement, I have a poor-quality but topologically valid triangulation (e.g., `chilmesh.examples.annulus()`) and I want to run ADMESH's truss solver on it without regenerating the boundary. I should be able to call a single function with the existing mesh and a size function, and receive a higher-quality mesh whose outer and inner boundary points are byte-identical to the input.

**Why this priority**: This is the core capability. Without it the feature delivers nothing. Boundary preservation is the non-negotiable constraint that distinguishes warm-start from a fresh `admesh.triangulate()`.

**Independent Test**: Call the warm-start function on `chilmesh.examples.annulus()` with a constant size function. Verify (a) median element quality strictly improves, and (b) every boundary point in the output appears at the same `(x, y)` coordinates as the corresponding input boundary point — bit-exact equality, not approximate.

**Acceptance Scenarios**:

1. **Given** the row 1 annulus mesh with median quality ≈ 0.49 and 180 boundary points, **When** the user calls the warm-start optimizer with the annulus SDF and a constant size function, **Then** the returned mesh has median quality > 0.49 AND boundary point coordinates match the input bit-exactly (np.array_equal).
2. **Given** an input mesh whose boundary is not on the SDF zero set (a malformed input), **When** the warm-start optimizer is called, **Then** it raises a clear `ValueError` naming the offending vertex indices — it MUST NOT silently project them.
3. **Given** a valid input mesh and a size function that grades from fine near boundaries to coarse in the middle, **When** the warm-start optimizer runs, **Then** elements near the boundary in the output are smaller than elements in the middle, demonstrating the size function was respected.

---

### User Story 2 - Visible Side-by-Side Comparison in README (Priority: P2)

As a reader of the project README, I want to see a fifth row (or a separate figure) demonstrating the warm-start truss optimization on row 1's mesh, alongside the existing four rows (raw / FEM-smoothed / fresh-ADMESH / right-isoceles), so I can visually compare the warm-start approach to the alternatives at a glance.

**Why this priority**: Documents the new capability and makes the comparison concrete for users. Lower priority than the underlying capability, but high value for adoption.

**Independent Test**: Generate the new visualization, open the PNG, and confirm a row labeled "Warm-Start ADMESH Truss" exists with three columns (mesh / layers / quality) and that the boundary in column 1 of that row visually matches the boundary of the row 1 raw mesh.

**Acceptance Scenarios**:

1. **Given** the existing 4-row visualization script, **When** the warm-start variant is added, **Then** the resulting figure cleanly shows the new row with the same column conventions as the existing rows (left = wireframe, center = parula layers, right = cool quality).
2. **Given** the new row's mesh, **When** rendered alongside row 1, **Then** the boundary contours are visually indistinguishable (because they are bit-exact equal under the hood).

---

### User Story 3 - Generic Triangulation-In, Optimized-Mesh-Out (Priority: P1)

As an author of a downstream tool (MADMESHR, ADMESH-Domains adapter, custom workflow, ad-hoc notebook), I want to convert ANY triangulation through ADMESH's truss optimizer without writing the plumbing myself. I should be able to call one function with `(points, triangles, sdf, size_fn)` (or a CHILmesh equivalent) and receive an optimized triangulation back, regardless of where the input came from — bundled fixture, file load, my own numpy arrays, or a third-party mesh library.

**Why this priority**: Elevated to P1 (co-equal with US1) because the user explicitly stated this functionality MUST be extensible: "i should be able to convert any given triangulation via admesh." Without generic-input support, this is just a one-shot demo. With it, this is a reusable building block.

**Independent Test**: The function MUST accept inputs in **two equivalent forms**, and the test suite MUST exercise both:

- **Form A (CHILmesh-native)**: `optimize_with_admesh_truss(mesh: CHILmesh, sdf, size_fn, **kwargs) -> CHILmesh`
- **Form B (raw arrays)**: `optimize_with_admesh_truss_arrays(points: ndarray, triangles: ndarray, sdf, size_fn, **kwargs) -> tuple[ndarray, ndarray]`

The CHILmesh form MUST be a thin convenience wrapper around the raw-arrays form. Both must work on all four bundled fixtures (annulus, donut, block_o, structured) AND on at least one fort.14-loaded mesh from the test corpus AND on at least one synthetic `(points, triangles)` constructed in the test itself.

**Acceptance Scenarios**:

1. **Given** the donut fixture and a donut SDF / size pair, **When** the CHILmesh-form function is called, **Then** it returns a `CHILmesh` instance whose boundary points match the donut's input boundary bit-exactly AND median quality does not regress.
2. **Given** raw `(points, triangles)` numpy arrays for a square domain (constructed in the test, not loaded from a fixture), **When** the raw-array form is called with a square SDF, **Then** it returns optimized `(points, triangles)` arrays where the four corner points and the boundary chain are bit-exact preserved.
3. **Given** a CHILmesh loaded from a `.fort.14` file, **When** the CHILmesh form is called, **Then** the output preserves the loaded boundary and the function can roundtrip via `write_fort14` without error.
4. **Given** an unsupported input (e.g., a quad mesh), **When** the function is called, **Then** it raises a clear `NotImplementedError` referring to "triangle-only for v1".
5. **Given** a triangulation whose boundary is NOT consistent with the supplied SDF (i.e., user passed the wrong domain), **When** the function is called, **Then** it raises `ValueError` listing offending boundary indices, naming the SDF tolerance used, and surfacing the worst-case absolute SDF value seen — i.e., enough information to diagnose without re-running.

---

### Edge Cases

- **Boundary identification**: How does the function distinguish boundary points from interior points? The current `mesh.boundary_edges()` returns edge IDs, not vertex indices. The function MUST use a deterministic boundary-detection rule (e.g., vertices appearing in any boundary edge) and document it.
- **SDF mismatch**: If the input mesh's boundary is not on the SDF zero set (within numerical tolerance), the function MUST raise rather than silently letting truss forces drag boundaries off the input geometry.
- **Connectivity vs free re-triangulation**: ADMESH's distmesh loop **internally re-triangulates** at each iteration (Delaunay over the current point cloud). The output connectivity will therefore NOT match the input connectivity. The spec MUST clarify that "warm-start" means "use the input *points* as the initial distribution," not "preserve input *triangles*."
- **Convergence failure**: If the truss loop fails to converge within max iterations, the function MUST surface the warning and return the best-so-far mesh, not raise.
- **Mixed-element / quad inputs**: V1 supports triangles only. Quad / mixed inputs MUST raise `NotImplementedError` with a clear message.
- **Domain-boundary mismatch**: If the user passes an SDF that doesn't actually contain all the input boundary points (e.g., wrong domain), the function MUST raise rather than trying to "fix" the points.
- **Degenerate triangles in input**: If the input mesh has triangles with zero or negative area, the function MUST raise rather than feeding them to ADMESH.
- **Tiny input meshes**: Inputs with fewer than ~10 interior points may not benefit from truss optimization. The function MUST still run but issue a warning if interior count is below a threshold.

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
- **FR-011**: The function MUST surface a `RuntimeWarning` (not raise) if the truss loop hits the max-iteration cap without converging, and still return the best-so-far mesh.
- **FR-012**: A demo script (extending or sibling to `generate_4row_admesh.py`) MUST render the warm-start result alongside the existing four rows so users can visually compare. Whether this is an additional row in the existing PNG or a separate PNG is a design choice deferred to the plan phase.
- **FR-013**: The demo script MUST enforce fail-loud assertions before saving any output:
  - **V_BND**: Output boundary points are bit-exact equal to input boundary points.
  - **V_QI**: Output median quality is strictly greater than input median quality (this is a *demo* assertion, not a function-level guarantee — see SC-002 wording).
  - **V_CONN**: Output connectivity defines a valid triangulation (every triangle has positive area, no duplicate vertex indices).
- **FR-014**: The function and demo script MUST be importable / runnable without ADMESH installed system-wide; the existing project convention of `sys.path.insert(0, str(Path(__file__).parent / "ADMESH"))` (or equivalent) is acceptable. If ADMESH cannot be imported, the function MUST raise a clear `ImportError` naming the missing dependency.
- **FR-015**: The integration MUST work with the live ADMESH source tree (the cloned `/tmp/ADMESH` from `https://github.com/domattioli/ADMESH.git`); if upstream ADMESH has the broken `MeshOutput` import (current `main` HEAD), the spec author and the implementer MUST coordinate on whether to upstream a fix or pin to a known-good commit.
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

- **WarmStartInput**: The triple `(CHILmesh, SDF, SizeFn)` plus optional kwargs. Defines the optimizer's input contract.
- **BoundaryVertexSet**: The set of vertex indices in the input mesh that lie on the domain boundary, identified via `boundary_edges()` → `Edge2Vert`. Marked as `pfix` for ADMESH.
- **InteriorVertexSet**: All input vertex indices NOT in BoundaryVertexSet. These are the "free" nodes the truss loop will move.
- **TrussSolverConfig**: The set of optional kwargs forwarded to `distmesh2d` (max iter, tolerance, retriangulation cadence, force scaling). Has documented defaults.
- **WarmStartOutput**: A new `CHILmesh` instance built from the optimized point set and ADMESH's final retriangulation. Indexes preserved for boundary points.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: On the bundled annulus fixture (`chilmesh.examples.annulus()`) with a constant size function, the warm-start optimizer produces a mesh whose median element quality is at least **0.55** (compared to the input's ≈ 0.49 — a ≥12% relative improvement). This is the demo benchmark; the function-level requirement is non-degradation (FR-013 V_QI).
- **SC-002**: On the same fixture, **100 %** of input boundary points appear in the output at bit-exact equal `(x, y)` coordinates (np.array_equal == True). This is the boundary-preservation guarantee.
- **SC-003**: The warm-start optimizer completes on a 580-element mesh (annulus) in under **30 seconds** of wall time on a developer laptop without ADMESH-side parallelization.
- **SC-004**: The demo visualization renders without invoking any of `_initial_distribution`, the rejection-sampling code path, or any ADMESH internals that generate fresh points — verified via assertion or instrumentation, not vibes.
- **SC-005**: All four bundled CHILmesh fixtures (annulus, donut, block_o, structured) can be passed to the warm-start optimizer without crashing. Quality may not strictly improve on every fixture (some are already excellent), but the function MUST NOT regress median quality below input on any fixture.
- **SC-006**: The README (or a sibling docs page) clearly explains when to use warm-start vs `admesh.triangulate()`-from-scratch, in plain language, with one or two sentences of guidance.
- **SC-007 (extensibility)**: The test suite MUST include at least three input-source variations (bundled fixture, fort.14-loaded mesh, raw numpy arrays from a non-CHILmesh source) and verify all three produce valid optimized output with bit-exact boundary preservation. If any of these three input forms fails, the feature is not done.
- **SC-008 (extensibility)**: The test suite MUST include at least two distinct domains (annulus + one other, e.g., a square or an L-shape) to verify the function is genuinely domain-agnostic, not just hard-coded to annulus parameters.

## Assumptions

- ADMESH source code is reachable at `https://github.com/domattioli/ADMESH.git` and the implementation team can clone a known-good commit if `main` is broken.
- The CHILmesh project's existing convention of `sys.path.insert(0, "ADMESH")` for in-tree imports continues to work.
- The annulus SDF used by the existing row 3 / 4 pipeline (`admesh.domains.ANNULUS`) accepts the row 1 boundary points within tolerance — i.e., row 1's boundary is in fact on the annulus zero set. (This will be verified during the spec→plan transition.)
- "Truss algorithm" in the user's request refers specifically to ADMESH's distmesh-style spring/repulsion solver (`distmesh2d` and friends), not a generic structural-mechanics solver. This is the authoritative interpretation throughout this spec.
- Boundary preservation is a hard constraint, not soft: bit-exact equality is the contract. ADMESH's `pfix` mechanism already supports this; we are not asking ADMESH to be modified, only to be invoked with the right inputs.
- V1 scope is triangle meshes only. Quad and mixed-element warm-start is explicitly deferred.
- The CHILmesh single-branch policy applies (`planning-optimize_modernize` is the working branch); this spec's directory `specs/005-admesh-warm-start-truss/` is added to that branch without creating a separate `005-*` git branch.
- The downstream goal of the four-row pipeline (right-isoceles smoothing for tri→quad fusion) is unchanged. The warm-start row is *additive*, not a replacement for any existing row.
