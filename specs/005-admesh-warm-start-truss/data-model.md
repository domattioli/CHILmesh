# Phase 1 Data Model: ADMESH Warm-Start Truss

**Date**: 2026-05-02
**Spec**: [spec.md](spec.md) | **Plan**: [plan.md](plan.md) | **Research**: [research.md](research.md)

This is a single-function library feature plus a demo script. **No persistent data model** — no databases, no schemas, no API entities. What follows are the runtime entities the adapter manipulates.

---

## Overview

The adapter is a pure function: triangulation in, optimized triangulation out. The runtime entities are the inputs, the validated intermediate state, the truss loop's working memory, and the output. Plus the four-row figure the demo produces.

---

## E1. WarmStartInput (function arguments)

The contract surface for `optimize_with_admesh_truss_arrays(...)`.

| Field | Type | Required | Source | Purpose |
|-------|------|----------|--------|---------|
| `points` | `ndarray[N, 2 \| 3]` | yes | caller | Vertex coordinates. If shape is (N,3), the z column is preserved through but never used by the truss. |
| `triangles` | `ndarray[M, 3]` | yes | caller | Triangle vertex indices. Quads/mixed → `NotImplementedError`. |
| `sdf` | `Callable[[(K,2) array], (K,) array]` | yes | caller | Signed distance function (negative inside). |
| `size_fn` | `Callable[[(K,2) array], (K,) array] \| None` | yes (None allowed) | caller | Element size function; None means uniform. |
| `boundary_indices` | `ndarray[B] \| None` | no | caller / adapter | If None, inferred via "edges in exactly one triangle". |
| `h0` | `float` | no (inferred) | caller / adapter | Target edge length. If None, computed as mean input edge length. |
| `bbox` | `tuple[float×4] \| None` | no (inferred) | caller / adapter | If None, inferred from input points (with margin). |
| `**kwargs` | various | no | caller | Forwarded to `distmesh2d_warmstart` (`dptol`, `ttol`, `Fscale`, `deltat`, `geps_factor`, `niter`, `seed`). See R3. |

The CHILmesh form `optimize_with_admesh_truss(mesh, sdf, size_fn, **kwargs) -> CHILmesh` is a thin wrapper that extracts `points`, `triangles`, and `boundary_indices` from the CHILmesh and passes them to the array form.

---

## E2. ValidatedInput (post-validation, pre-truss)

After running validation (FR-005, 006, 007) but before entering the truss loop, the adapter holds a normalized internal representation:

| Field | Type | Computed from | Purpose |
|-------|------|--------------|---------|
| `points_xy` | `ndarray[N, 2]` | `points[:, :2]` | 2D coordinates only (truss is 2D). |
| `triangles` | `ndarray[M, 3]` | input | Validated triangle-only. |
| `boundary_indices` | `ndarray[B]` | input or inferred | Validated to have `\|sdf(point)\| < 1e-6` for every entry. |
| `interior_indices` | `ndarray[N-B]` | `setdiff1d(arange(N), boundary_indices)` | The free nodes for the truss. |
| `boundary_xy` | `ndarray[B, 2]` | `points_xy[boundary_indices]` | Will become `pfix` for ADMESH. |
| `interior_xy` | `ndarray[N-B, 2]` | `points_xy[interior_indices]` | The warm-start free points. |
| `h0` | `float` | provided or inferred from input edge lengths | Target edge length. |
| `bbox` | `tuple[float×4]` | provided or inferred from input points | Bounding box. |
| `truss_kwargs` | `dict` | input kwargs | Forwarded verbatim to vendored loop. |

Validation rules (raise `ValueError` / `NotImplementedError` with offending indices listed):

| Rule | Check | Error |
|------|-------|-------|
| V_TRI | `triangles.shape[1] == 3` | `NotImplementedError("triangle-only for v1")` |
| V_AREA | All triangles have positive signed area | `ValueError("degenerate or negative-area triangles at indices [...]")` |
| V_BND_SDF | `\|sdf(boundary_xy)\| < 1e-6` for all boundary points | `ValueError("boundary not on SDF zero set; max abs sdf=X at index Y")` |
| V_INTERIOR_INSIDE | `sdf(interior_xy) < 0` (warning, not raise) | `RuntimeWarning("N interior points have sdf > 0; truss may push them onto boundary")` |
| V_BBOX | All input points within `bbox` | `ValueError("points outside bbox at indices [...]")` |

---

## E3. TrussLoopState (working memory inside `_vendor_admesh_truss.distmesh2d_warmstart`)

The vendored truss loop holds:

| Field | Type | Initial value | Purpose |
|-------|------|---------------|---------|
| `p` | `ndarray[N, 2]` | `np.vstack([boundary_xy, interior_xy])` | Working points. Boundary at indices 0..B-1, interior at B..N-1. |
| `nfix` | `int` | `B` (= len(boundary_xy)) | Count of pinned points. `Ftot[:nfix] = 0.0` every iter. |
| `pold` | `ndarray[N, 2]` | `np.full_like(p, inf)` | Previous-iteration positions (for ttol check). |
| `t` | `ndarray[Mp, 3]` | `np.empty((0,3), dtype=int64)` | Current triangulation; rebuilt via Delaunay when ttol exceeded. |
| `bars` | `ndarray[K, 2]` | `np.empty((0,2), dtype=int64)` | Edge list derived from `t` (truss bars). |

State transitions per iteration (verbatim from upstream `distmesh2d` lines 140-220):

```text
1. If max_movement / h0 > ttol:
     pold = p.copy()
     tri = Delaunay(p)
     t = tri.simplices, filtered by centroid SDF
     bars = unique edges from t
2. Compute truss forces F along each bar.
3. Sum forces per node (Ftot).
4. Zero forces at pfix indices: Ftot[:nfix] = 0.
5. p_new = p + deltat * Ftot
6. Project boundary-violating interior nodes back onto SDF zero set.
7. Check stopping criterion: max(interior movement) / h0 < dptol → break.
8. p = p_new
```

The vendor copy is **byte-identical** to upstream lines 140-end of `distmesh2d`. The only divergence is the preamble (steps 1-3 of upstream replaced with our warm-start setup).

---

## E4. WarmStartOutput (return value)

For the array form (`optimize_with_admesh_truss_arrays`):

| Field | Type | Source | Property |
|-------|------|--------|----------|
| `points_out` | `ndarray[Np, 2]` | truss loop final state | Boundary at indices 0..B-1 (bit-exact equal to `boundary_xy`); interior at B..Np-1 |
| `triangles_out` | `ndarray[Mp, 3]` | truss loop final Delaunay | Re-triangulated; *not* equal to input triangles |

For the CHILmesh form: a fresh `CHILmesh(connectivity=triangles_out, points=column_stack([points_out, zeros(Np)]))`.

**Index preservation guarantee**: If a caller passed `boundary_indices=[i0, i1, ..., iB-1]` (in any order), the output's first B rows correspond to those points in *that* order (i.e., `points_out[k] == input_points[boundary_indices[k]]`).

---

## E5. Demo: FourRowFigure (refreshed for new layout)

The demo script's runtime entities, replacing the spec-004 4-row layout:

| Field | Value |
|-------|-------|
| `shape` | 4 rows × 3 columns = 12 axes |
| `figsize` | `(15, 18)` inches (matches existing) |
| `dpi` | 100 (matches existing) |
| `output_path` | `tests/output/annulus_quickstart.png` (UNCHANGED) |
| `suptitle` | `"CHILmesh × ADMESH: Warm-Start Truss Optimization Pipeline"` |

Per-row mesh provenance (the new layout from Q2=d):

| Row | Source | Smoother applied |
|-----|--------|------------------|
| 1 | `chilmesh.examples.annulus()` | None (raw Delaunay) |
| 2 | Row 1 → `optimize_with_admesh_truss(row1, ANNULUS_SDF, size_fn=None)` | ADMESH truss (warm-start) |
| 3 | Row 2 → `mesh.smooth_mesh(method='fem', acknowledge_change=True)` | CHILmesh FEM (`direct_smoother`) |
| 4 | Row 2 → `smooth_for_quadrangulation(row2.points, row2.triangles, ANNULUS_SDF, ...)` | ADMESH right-isoceles |

Note rows 3 and 4 both branch off row 2 — they are siblings, not sequential.

Per-cell rendering rule (preserved from spec 004):

| Cell | Renderer | Color source |
|------|----------|--------------|
| `(row, 0)` | `ax.triplot(x, y, t, color='black')` | None (no fill) |
| `(row, 1)` | `ax.fill(...)` per element | `parula_cmap(elem_to_layer[i] / max(1, n_layers-1))` |
| `(row, 2)` | `ax.fill(...)` per element | `cool_r_cmap(norm(quality[i]))` (cool_r = red→blue, see spec 004 Q5 history) |

Subplot titles (per spec 004 contract format `"<Row Label>\n<Column Label>"`):

| Row | Row Label |
|-----|-----------|
| 1 | `Raw Delaunay` |
| 2 | `+ ADMESH Truss (warm-start)` |
| 3 | `Row 2 + FEM Smoother` |
| 4 | `Row 2 + Right-Isoceles Smoother` |

| Column | Column Label |
|--------|-------------|
| 1 | `Mesh` |
| 2 | `Layers` |
| 3 | `Quality` |

---

## State Transitions

The whole feature has a single linear runtime flow per call:

```text
LOAD ValidatedInput
    | (raise ValueError / NotImplementedError on validation failure)
    v
INVOKE _vendor_admesh_truss.distmesh2d_warmstart(...)
    | (RuntimeWarning if niter cap hit without convergence)
    v
COMPUTE q_in (median quality of input)
COMPUTE q_out (median quality of output)
    | if q_out < q_in:
    |     emit RuntimeWarning("warm-start truss did not improve quality...")
    |     return input mesh (FR-011, Q4=b)
    v
RETURN WarmStartOutput
```

The demo script's flow:

```text
GENERATE row1 (chilmesh.examples.annulus)
GENERATE row2 (optimize_with_admesh_truss on row1)
GENERATE row3 (FEM smooth on row2.copy)
GENERATE row4 (right-isoceles on row2.copy)
RUN ASSERTIONS:
    V_BND  - row2 boundary == row1 boundary (np.array_equal)
    V_BND_PROP - row3, row4 boundaries == row2 boundary (or smoother's documented tolerance)
    V_QI   - median(row2 quality) > median(row1 quality)
    V_CONN - all four meshes have positive triangle areas, no duplicate vertex indices
    V_CHAIN - row3 input was row2.points/triangles; row4 input was row2.points/triangles
    | (any failure → RuntimeError, PNG NOT written)
    v
RENDER 12 subplots
SAVE PNG to tests/output/annulus_quickstart.png
EXIT 0
```

---

## Conclusion

The data model is **runtime-only and intentionally minimal**. The adapter holds a few hundred bytes of working state at most; nothing is persisted. The demo's only persistent artifact is the regenerated PNG at the same path the README already references.

**Status**: Ready for `contracts/api-contract.md` and `contracts/visualization-output.md`.
