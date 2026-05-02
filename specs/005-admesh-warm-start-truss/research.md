# Phase 0 Research: ADMESH Warm-Start Truss Optimization

**Date**: 2026-05-02
**Spec**: [spec.md](spec.md)
**Plan**: [plan.md](plan.md)

This document captures the technical investigation done before design. Each finding (R1-R5) maps to a constraint or decision in `plan.md`.

---

## R1: ADMESH `main` HEAD has a broken import — `MeshOutput`

### Question

Can we call `admesh.triangulate()` (the public entry point in `admesh/routine.py`) from CHILmesh?

### Investigation

```bash
cd /tmp/ADMESH && git rev-parse HEAD
# 05bc68fc81060f7d710b8f4abb2cc382f85df33f
```

The file `admesh/routine.py` (line 21) imports:

```python
from admesh.distmesh import MeshOutput, SizeFn, distmesh2d, distmesh2d_admesh
```

But `admesh/distmesh.py` defines only `SizeFn`, `distmesh2d`, `fixmesh`, `_initial_distribution`, `_rejection_method`. There is no `MeshOutput` class and no `distmesh2d_admesh` function on `main`.

Searching git history:

```text
f21e340 Restore missing ADMESH-variant distmesh code (MeshOutput, SizeFn, distmesh2d_admesh)
855670e Fix test_issue_11_ring_sorting: correct numpy array shapes and domain access
```

These commits exist on the `daily-issue-fixing` branch but **have not been merged to `main`**. The user's pinned commit `05bc68f` does not contain them.

### Verification

```bash
python -c "import sys; sys.path.insert(0, '/tmp/ADMESH'); from admesh import triangulate"
# ImportError: cannot import name 'MeshOutput' from 'admesh.distmesh'
```

```bash
python -c "import sys; sys.path.insert(0, '/tmp/ADMESH'); from admesh.distmesh import distmesh2d; print(distmesh2d)"
# <function distmesh2d at 0x...>
```

### Decision

`distmesh2d` is **directly importable**. The CHILmesh adapter calls `admesh.distmesh.distmesh2d` directly and never imports from `admesh.routine`. This sidesteps the broken state without needing a fix upstream.

**Cross-repo issue**: This is candidate ADMESH-A in spec.md's cross-repo tracking section. The user files it on `domattioli/ADMESH/issues`.

---

## R2: `pfix` bit-exact preservation — proof from source

### Question

The spec promises bit-exact equality (`np.array_equal`) on the boundary subset. Can ADMESH's `pfix` mechanism actually deliver this?

### Investigation

Source of `distmesh2d` (`admesh/distmesh.py` lines 117-132, 170-172):

```python
# Lines 117-132: Fixed points get the lowest indices
if pfix is not None:
    pfix_arr = np.asarray(pfix, dtype=float).reshape(-1, 2)
    if pfix_arr.size:
        # Drop any free points that coincide with fixed points.
        if len(p):
            dist_to_fix = np.min(...)
            p = p[dist_to_fix > geps]
        p = np.vstack([pfix_arr, p])    # <-- pfix at indices 0..nfix-1
        nfix = len(pfix_arr)
    else:
        nfix = 0
else:
    nfix = 0

# Lines 170-172: Force at fixed indices is zeroed every iteration
if nfix > 0:
    Ftot[:nfix] = 0.0

# (next line) p_new = p + deltat * Ftot
# Therefore p_new[:nfix] = p[:nfix] + 0 = p[:nfix] exactly.
```

### Risk: Re-triangulation reordering

Inside the truss loop:

```python
tri = Delaunay(p)
t_all = tri.simplices
```

scipy's `Delaunay` produces triangles indexing into `p`. If `Delaunay` reorders `p` internally (it doesn't — it returns `tri.simplices` indexing the input array unchanged), our pfix invariant could break. **Verified**: `scipy.spatial.Delaunay(p)` does not modify `p`; `tri.points` is `p` (same object) and indices in `tri.simplices` reference the input array's order. pfix preservation survives Delaunay.

### Decision

`pfix` is bit-exact preserved by inspection of source. The spec's `np.array_equal` requirement is satisfiable. We add a regression test (V_BND in FR-013, plus a unit test in `test_admesh_warmstart.py`) that exercises this end-to-end so any future scipy or ADMESH change that breaks the invariant fails loudly.

**Cross-repo issue**: Document this contract upstream — candidate ADMESH-C.

---

## R3: `distmesh2d` public signature — what tunables to forward

### Question

Per FR-009, the adapter forwards optional truss-solver parameters with documented defaults. What are they?

### Investigation

```python
from admesh.distmesh import distmesh2d
import inspect
print(inspect.signature(distmesh2d))
# (fd, fh, h0, bbox, pfix=None, *, dptol=0.001, ttol=0.1, Fscale=1.2,
#  deltat=0.2, geps_factor=0.001, niter=500, seed=0)
```

| Param | Default | Meaning |
|-------|---------|---------|
| `fd` | required | Signed distance function (callable) |
| `fh` | required* | Mesh size function (callable; None means uniform) |
| `h0` | required | Target edge length (sets initial lattice spacing) |
| `bbox` | required | (xmin, ymin, xmax, ymax) bounding box |
| `pfix` | `None` | Fixed points (M, 2) — never moved |
| `dptol` | `0.001` | Stopping criterion: relative interior-node movement tolerance |
| `ttol` | `0.1` | Re-triangulation trigger (relative node movement) |
| `Fscale` | `1.2` | Truss internal-pressure factor (> 1) |
| `deltat` | `0.2` | Euler time step for force-displacement update |
| `geps_factor` | `0.001` | Boundary tolerance: `geps = geps_factor * h0` |
| `niter` | `500` | Maximum iterations |
| `seed` | `0` | RNG seed for the rejection step (deterministic) |

\* `fh` accepts None internally; we surface it as `None` default in the adapter too, meaning uniform.

### Decision

The CHILmesh adapter exposes all of these as kwargs and forwards them through. Defaults match upstream so callers who pass nothing get the same behavior as `admesh.triangulate()` would.

**Note on `seed`**: Even though we skip `_rejection_method` (which is the only consumer of the seed), we still accept and forward `seed` for forward-compatibility — if ADMESH adds another stochastic step in the future, our adapter still controls it.

---

## R4: Can we inject custom initial points without modifying ADMESH?

### Question

Is there any way (monkey-patch, optional argument, side door) to make `distmesh2d` skip its `_initial_distribution` + `_rejection_method` steps and use our supplied points instead?

### Investigation

`distmesh2d` body lines 110-115:

```python
# 1. Initial lattice + drop points outside the domain.
p = _initial_distribution(bbox, h0)
p = p[fd(p) < geps]

# 2. Probability-based rejection by fh.
p = _rejection_method(p, fh, rng)

# 3. Fixed points first ...
```

The local variable `p` is hardcoded to the lattice generator's output. There is no parameter, no sentinel value, no callable injection point that lets a caller substitute a different starting `p`. Monkey-patching `_initial_distribution` is possible (it's module-scoped) but fragile and pollutes test isolation.

### Tested workaround attempt

If we pass our existing points as `pfix`, they get pinned (never moved). We could pin only the boundary and let `_initial_distribution` generate fresh interior points — but that's fresh-from-bbox, which is exactly the existing row 3 behavior we're replacing. Doesn't satisfy "warm-start using row 1's interior."

### Decision

The vendor approach is necessary. We copy the inner truss loop (about 50 lines) into `_vendor_admesh_truss.py` with a header comment listing the source SHA. The only behavioral difference: skip lines 110-115 (initial distribution + rejection) and instead start with `p = np.vstack([pfix_boundary, interior_initial_from_caller])`.

**Cross-repo issue**: This is candidate ADMESH-B (add a public warm-start entry point). When that lands, we delete the vendor and call ADMESH directly.

---

## R5: Right-isoceles smoother (`smooth_for_quadrangulation`) — does it import cleanly?

### Question

Spec 004 row 4 used `admesh.quad_prep.smooth_for_quadrangulation`. Our new row 4 also uses it (applied to row 2 instead of row 3). Does the broken `routine.py` block this import?

### Investigation

```bash
python -c "import sys; sys.path.insert(0, '/tmp/ADMESH'); from admesh.quad_prep import smooth_for_quadrangulation; print(smooth_for_quadrangulation)"
# <function smooth_for_quadrangulation at 0x...>
```

`admesh.quad_prep` does not import from `admesh.routine`, so it is unaffected by the `MeshOutput` issue. Our new row 4 works exactly as in spec 004.

### Decision

No change to row 4's invocation pattern. The only difference vs. spec 004 is the *input mesh* (row 2 instead of row 3).

---

## R6 (bonus): Does scipy `Delaunay` preserve point ordering?

### Question

The truss loop calls `Delaunay(p)` repeatedly. If scipy reorders `p`, our `pfix` invariant breaks (R2). Have we verified this?

### Investigation

```python
import numpy as np
from scipy.spatial import Delaunay

p = np.random.RandomState(0).rand(50, 2)
tri = Delaunay(p)
print(tri.points is p)            # True — same object
print(np.array_equal(tri.points, p))  # True — bit-exact
print(tri.simplices.max() == len(p) - 1)  # True — indices are valid for original p
```

scipy `Delaunay` does NOT reorder, copy, or modify the input point array. `tri.points is p`. The simplices reference the input ordering directly.

### Decision

No mitigation needed; the property is contractual upstream. We document this assumption in the adapter docstring so any future scipy major-version bump triggers re-verification.

---

## Summary Table

| Finding | Status | Decision | Cross-Repo Issue |
|---------|--------|----------|------------------|
| R1: routine.py broken | Confirmed | Bypass `routine`; call `distmesh2d` directly | ADMESH-A |
| R2: pfix preservation | Confirmed bit-exact | Use as-is; add regression test | ADMESH-C |
| R3: distmesh2d signature | Documented | Forward all 12 kwargs | (none) |
| R4: no warm-start hook | Confirmed | Vendor the inner truss loop | ADMESH-B |
| R5: quad_prep clean | Confirmed | No changes to row 4 invocation | (none) |
| R6: Delaunay preserves order | Confirmed | Document assumption | (none) |

**Status**: Phase 0 complete. All assumptions for Phase 1 design are validated.
