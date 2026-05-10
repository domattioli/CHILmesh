# Phase 0 Research: ADMESH Warm-Start Truss Optimization

**Date**: 2026-05-02
**Spec**: [spec.md](spec.md) | **Plan**: [plan.md](plan.md)

Technical investigation before design. Each finding (R1-R5) maps to constraint/decision in `plan.md`.

---

## R1: ADMESH `main` HEAD has a broken import — `MeshOutput`

**Question**: Can we call `admesh.triangulate()` (public entry point in `admesh/routine.py`) from CHILmesh?

### Investigation

```bash
cd /tmp/ADMESH && git rev-parse HEAD
# 05bc68fc81060f7d710b8f4abb2cc382f85df33f
```

The file `admesh/routine.py` (line 21) imports:

```python
from admesh.distmesh import MeshOutput, SizeFn, distmesh2d, distmesh2d_admesh
```

`admesh/distmesh.py` defines only `SizeFn`, `distmesh2d`, `fixmesh`, `_initial_distribution`, `_rejection_method`. No `MeshOutput` or `distmesh2d_admesh` on `main`.

Git history:

```text
f21e340 Restore missing ADMESH-variant distmesh code (MeshOutput, SizeFn, distmesh2d_admesh)
855670e Fix test_issue_11_ring_sorting: correct numpy array shapes and domain access
```

These commits on `daily-issue-fixing` — not merged to `main`. Pinned commit `05bc68f` doesn't contain them.

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

`distmesh2d` directly importable. Adapter calls `admesh.distmesh.distmesh2d` directly; never imports `admesh.routine`. **Cross-repo issue**: ADMESH-A.

---

## R2: `pfix` bit-exact preservation — proof from source

**Question**: Can ADMESH's `pfix` mechanism deliver `np.array_equal` bit-exact preservation on boundary subset?

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

`scipy.spatial.Delaunay(p)` does not modify `p`; `tri.points is p` (same object). pfix preservation survives Delaunay — verified.

### Decision

`pfix` bit-exact preserved by source inspection. Regression test (V_BND in FR-013 + unit test) exercises end-to-end. **Cross-repo issue**: ADMESH-C.

---

## R3: `distmesh2d` public signature — tunables to forward

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

\* `fh` accepts None = uniform.

### Decision

Adapter exposes all kwargs, forwards through. Defaults match upstream. `seed` forwarded for forward-compatibility even though `_rejection_method` is skipped.

---

## R4: Can we inject custom initial points without modifying ADMESH?

**Question**: Any way (monkey-patch, optional arg, side door) to skip `_initial_distribution` + `_rejection_method` and use our supplied points?

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

Local var `p` hardcoded to lattice generator output. No injection point. Monkey-patching fragile. Passing existing points as `pfix` pins them (doesn't satisfy warm-start using row 1's interior).

### Decision

Vendor approach necessary. Copy inner truss loop (~50 lines) into `_vendor_admesh_truss.py` with source SHA header. Only change: skip lines 110-115, start with `p = np.vstack([pfix_boundary, interior_initial_from_caller])`. **Cross-repo issue**: ADMESH-B. Delete vendor when that lands.

---

## R5: Right-isoceles smoother (`smooth_for_quadrangulation`) — does it import cleanly?

**Question**: Does broken `routine.py` block `admesh.quad_prep.smooth_for_quadrangulation` import?

### Investigation

```bash
python -c "import sys; sys.path.insert(0, '/tmp/ADMESH'); from admesh.quad_prep import smooth_for_quadrangulation; print(smooth_for_quadrangulation)"
# <function smooth_for_quadrangulation at 0x...>
```

`admesh.quad_prep` doesn't import `admesh.routine` — unaffected by `MeshOutput` issue.

### Decision

No change to row 4 invocation. Only difference vs spec 004: input mesh is row 2 instead of row 3.

---

## R6 (bonus): Does scipy `Delaunay` preserve point ordering?

**Question**: If scipy reorders `p`, `pfix` invariant breaks. Verified?

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

scipy `Delaunay` does NOT reorder/copy/modify input array. `tri.points is p`. Simplices reference input ordering directly.

### Decision

No mitigation needed. Document assumption in adapter docstring — triggers re-verification on future scipy major-version bump.

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

**Status**: Phase 0 complete. All Phase 1 design assumptions validated.
