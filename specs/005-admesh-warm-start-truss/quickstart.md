# Quickstart: Using the ADMESH Warm-Start Truss Adapter

**Date**: 2026-05-02
**Spec**: [spec.md](spec.md) | **Plan**: [plan.md](plan.md) | **API contract**: [contracts/api-contract.md](contracts/api-contract.md)

This guide shows three worked examples covering the extensibility contract (FR-018):

1. Bundled annulus (the simplest case — uses the CHILmesh form)
2. Bundled donut (proves domain-agnostic — different SDF, same API)
3. Raw `(points, triangles)` from a non-CHILmesh source (proves source-agnostic — uses the array form)

---

## Installation prerequisites

```bash
pip install -e /path/to/CHILmesh        # this repo
pip install -e /path/to/ADMESH          # sibling repo, pinned to SHA 05bc68f
```

Or in development:

```bash
git clone https://github.com/domattioli/CHILmesh.git
cd CHILmesh
pip install -e .
git clone https://github.com/domattioli/ADMESH.git ./ADMESH
cd ADMESH && git checkout 05bc68f && cd ..
```

The CHILmesh adapter resolves ADMESH via `sys.path.insert(0, "ADMESH")` at import time if not pip-installed.

---

## Example 1: Bundled annulus (CHILmesh form)

The minimal one-liner. This is also exactly what Row 2 of the README's 4-row visualization does.

```python
import chilmesh
import matplotlib.pyplot as plt

# Load the bundled annulus fixture (poor-quality random Delaunay)
mesh = chilmesh.examples.annulus()
print(f"Input median quality: {mesh.elem_quality()[0].mean():.4f}")
# Input median quality: 0.4906

# Define the annulus SDF (negative inside the ring)
def annulus_sdf(p):
    import numpy as np
    r = np.linalg.norm(p, axis=1)
    return np.maximum(r - 1.0, 0.5 - r)  # outer R=1, inner r=0.5

# Run warm-start truss optimization (boundary preserved bit-exactly)
optimized = chilmesh.optimize_with_admesh_truss(
    mesh,
    sdf=annulus_sdf,
    size_fn=None,    # uniform size (= h0)
    seed=0,          # deterministic for repeatable demos
)

print(f"Output median quality: {optimized.elem_quality()[0].mean():.4f}")
# Output median quality: ~0.60 (per SC-001 Q5=b benchmark)

# Verify boundary preservation
import numpy as np
boundary_edges = mesh.boundary_edges()
boundary_indices = np.unique(mesh.adjacencies['Edge2Vert'][boundary_edges].flatten())
assert np.array_equal(
    optimized.points[:len(boundary_indices), :2],
    mesh.points[boundary_indices, :2]
), "Boundary preservation failed!"

# Plot the result
optimized.plot()
plt.show()
```

**What happened under the hood**:
1. `optimize_with_admesh_truss` extracted `(points, triangles, boundary_indices)` from the CHILmesh.
2. Validated triangle-only, positive areas, boundary-on-SDF, all points within bbox.
3. Called the vendored `distmesh2d_warmstart` with our boundary as `pfix` and our interior as the warm-start free points (skipping `_initial_distribution`).
4. Truss loop ran ~50-200 iterations until convergence (default `dptol=1e-3`).
5. Wrapped the output `(points, triangles)` in a fresh CHILmesh.

---

## Example 2: Bundled donut (different domain, same API)

Same function, completely different geometry. This demonstrates FR-017 (domain-agnostic).

```python
import chilmesh
import numpy as np

# Load the bundled donut fixture
mesh = chilmesh.examples.donut()
print(f"Donut input quality: {mesh.elem_quality()[0].mean():.4f}")

# Define the donut SDF (donut = annulus with extra detail near rings)
# The exact SDF depends on the donut's geometry; here's a representative example:
def donut_sdf(p):
    """Donut shape: outer ring R=1.0, inner hole r=0.4."""
    r = np.linalg.norm(p, axis=1)
    return np.maximum(r - 1.0, 0.4 - r)

# Define a graded size function: finer near boundaries
def donut_size_fn(p):
    """Smaller elements near the rings, larger in the middle."""
    r = np.linalg.norm(p, axis=1)
    d_outer = 1.0 - r          # distance to outer ring
    d_inner = r - 0.4          # distance to inner ring
    d_boundary = np.minimum(d_outer, d_inner).clip(min=0)
    # min element size 0.05, max element size 0.15
    return 0.05 + 0.10 * (d_boundary / 0.3).clip(max=1.0)

# Run warm-start with a graded size function this time
optimized = chilmesh.optimize_with_admesh_truss(
    mesh,
    sdf=donut_sdf,
    size_fn=donut_size_fn,
    h0=0.1,           # explicit base spacing
    seed=0,
)

print(f"Donut output quality: {optimized.elem_quality()[0].mean():.4f}")
# Quality should not regress (FR-011 graceful-degradation guarantee).

# Verify the size function was respected: elements near the boundary should be
# smaller than elements in the middle.
elem_centroids = optimized.points[optimized.connectivity_list].mean(axis=1)[:, :2]
elem_radii = np.linalg.norm(elem_centroids, axis=1)
elem_areas = optimized.elem_quality()[0]  # quality, but we want sizes — use a different call
# (See `chilmesh.signed_area()` for actual area computation; pseudocode here.)
```

**Key point**: The exact same function works on a completely different domain. The only thing that changed is the `sdf` and `size_fn` callables — the library code knows nothing about annuli vs donuts.

---

## Example 3: Raw arrays from a non-CHILmesh source (array form)

This is the most general case — prove you can warm-start *any* triangulation, not just ones produced by CHILmesh.

```python
from chilmesh import optimize_with_admesh_truss_arrays
import numpy as np
from scipy.spatial import Delaunay

# Generate a poor-quality square triangulation from scratch
np.random.seed(42)
n_interior = 80
side = 2.0  # unit square from [-1, 1]^2

# Boundary points: 40 evenly-spaced points around the square boundary
n_per_side = 10
edge = np.linspace(-1, 1, n_per_side, endpoint=False)
boundary = np.concatenate([
    np.column_stack([edge, np.full(n_per_side, -1)]),  # bottom
    np.column_stack([np.full(n_per_side, 1), edge]),   # right
    np.column_stack([edge[::-1], np.full(n_per_side, 1)]),  # top
    np.column_stack([np.full(n_per_side, -1), edge[::-1]])  # left
])
n_boundary = len(boundary)

# Interior points: random
interior = np.random.uniform(-0.95, 0.95, size=(n_interior, 2))

# Combine and triangulate
all_points = np.vstack([boundary, interior])
tri = Delaunay(all_points)
triangles = tri.simplices

# Boundary indices = first n_boundary in our concatenated array
boundary_indices = np.arange(n_boundary)

# Define the square SDF
def square_sdf(p):
    return np.max(np.abs(p) - 1.0, axis=1)

# Optimize with warm-start truss
points_out, triangles_out = optimize_with_admesh_truss_arrays(
    points=all_points,
    triangles=triangles,
    sdf=square_sdf,
    size_fn=None,                       # uniform
    boundary_indices=boundary_indices,  # explicit (not inferred)
    h0=0.2,
    bbox=(-1.05, -1.05, 1.05, 1.05),
    seed=0,
)

# Verify the four corners and boundary chain are bit-exact preserved
assert np.array_equal(points_out[:n_boundary], boundary), "Boundary not preserved!"
print(f"Optimized: {len(points_out)} points, {len(triangles_out)} triangles")
print(f"Boundary preserved: {n_boundary} points bit-exact")
```

**What this proves**:
- The function accepts plain numpy arrays — no CHILmesh dependency on the *input* side.
- A user-supplied SDF for a brand-new domain (the unit square) works without any library changes.
- An explicit `boundary_indices` argument bypasses the auto-detection (handy when you already know which points are boundary).

---

## Common patterns

### Pattern: Inferring boundary from the triangulation

If you don't know which indices are boundary, omit `boundary_indices` and the adapter infers them via "edges in exactly one triangle":

```python
points_out, triangles_out = optimize_with_admesh_truss_arrays(
    points=my_points,
    triangles=my_triangles,
    sdf=my_sdf,
    # boundary_indices omitted — auto-detected
)
```

This works for any closed-domain triangulation but assumes the input is well-formed (no duplicate triangles, no self-intersecting edges).

### Pattern: Tuning truss convergence

For a stubborn input where the default 500 iterations isn't enough:

```python
optimized = chilmesh.optimize_with_admesh_truss(
    mesh, sdf, size_fn,
    niter=2000,        # more iterations
    dptol=1e-4,        # tighter convergence
    Fscale=1.3,        # stronger truss pressure
    deltat=0.15,       # smaller time step (more stable)
)
```

### Pattern: Suppressing the non-degradation guard

If you specifically want to see what the truss *would* produce, even if it's worse:

```python
optimized = chilmesh.optimize_with_admesh_truss(
    mesh, sdf, size_fn,
    enforce_non_degradation=False,
)
# May return a worse mesh, but no surprise input-substitution.
```

### Pattern: Catching specific errors

```python
try:
    optimized = chilmesh.optimize_with_admesh_truss(mesh, sdf, size_fn)
except NotImplementedError:
    # Mixed-element or quad input
    pass
except ValueError as e:
    # Boundary off SDF, degenerate triangles, or out-of-bbox points
    if "boundary not on SDF zero set" in str(e):
        # User passed wrong domain
        ...
    elif "non-positive signed area" in str(e):
        # Mesh is malformed — fix CCW orientation first
        ...
```

### Pattern: Asserting bit-exact boundary preservation

```python
import numpy as np

# Identify boundary in input
boundary_edges = mesh.boundary_edges()
input_boundary_indices = np.unique(
    mesh.adjacencies['Edge2Vert'][boundary_edges].flatten()
)
input_boundary_xy = mesh.points[input_boundary_indices, :2]

optimized = chilmesh.optimize_with_admesh_truss(mesh, sdf, size_fn)

# Boundary appears at indices 0..B-1 of output (in input order)
output_boundary_xy = optimized.points[:len(input_boundary_indices), :2]

assert np.array_equal(input_boundary_xy, output_boundary_xy), \
    f"Boundary preservation violated; max delta = " \
    f"{np.max(np.abs(input_boundary_xy - output_boundary_xy)):.2e}"
```

---

## Verification checklist

After running the adapter, verify:

- [ ] Output is a valid `CHILmesh` instance (CHILmesh form) or a `(points, triangles)` tuple (array form).
- [ ] No `RuntimeError` was raised; warnings are acceptable but tell you something happened.
- [ ] `np.array_equal(output_boundary, input_boundary)` is True for the boundary subset.
- [ ] `output.elem_quality()[0].mean() >= input.elem_quality()[0].mean()` (with `enforce_non_degradation=True`, default).
- [ ] The output's connectivity is a valid triangulation (positive areas everywhere).

If any of these fail, file an issue at `https://github.com/domattioli/chilmesh/issues` with a minimal reproducer.
