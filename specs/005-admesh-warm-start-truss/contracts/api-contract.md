# Contract: Adapter API

**Date**: 2026-05-02
**Phase**: 1 — Design & Contracts
**Type**: Public Python API

---

## Purpose

Exact function signatures, validation contract, error catalog, and behavioral guarantees for both warm-start adapter entry points.

---

## Module Location

```python
# Importable from chilmesh top-level (re-exported in __init__)
from chilmesh import optimize_with_admesh_truss              # high-level CHILmesh form
from chilmesh import optimize_with_admesh_truss_arrays       # low-level array form

# Or from the explicit module
from chilmesh.admesh_warmstart import (
    optimize_with_admesh_truss,
    optimize_with_admesh_truss_arrays,
)
```

---

## Low-Level: `optimize_with_admesh_truss_arrays`

### Signature

```python
def optimize_with_admesh_truss_arrays(
    points: np.ndarray,                  # shape (N, 2) or (N, 3)
    triangles: np.ndarray,               # shape (M, 3)
    sdf: Callable[[np.ndarray], np.ndarray],
    size_fn: Optional[Callable[[np.ndarray], np.ndarray]] = None,
    *,
    boundary_indices: Optional[np.ndarray] = None,    # shape (B,)
    h0: Optional[float] = None,
    bbox: Optional[tuple[float, float, float, float]] = None,
    # forwarded to vendored truss loop (R3)
    dptol: float = 1e-3,
    ttol: float = 0.1,
    Fscale: float = 1.2,
    deltat: float = 0.2,
    geps_factor: float = 1e-3,
    niter: int = 500,
    seed: int = 0,
    # adapter-specific
    sdf_tolerance: float = 1e-6,
    enforce_non_degradation: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Optimize an existing triangulation using ADMESH's distmesh truss algorithm,
    preserving boundary points bit-exactly.

    Parameters
    ----------
    points : (N, 2) or (N, 3) ndarray
        Vertex coordinates. Z column (if present) is preserved as zero in the
        output but never used by the truss.
    triangles : (M, 3) ndarray of int
        Triangle vertex indices. Quad / mixed-element inputs raise NotImplementedError.
    sdf : callable
        Signed distance function: `sdf(points: (K,2)) -> (K,)`. Negative inside.
    size_fn : callable or None
        Element size function: `size_fn(points: (K,2)) -> (K,)`. None means uniform (h0).
    boundary_indices : (B,) ndarray of int or None
        Indices of points to pin (preserved bit-exactly). If None, inferred from
        triangulation (edges in exactly one triangle).
    h0 : float or None
        Target edge length. If None, computed as the mean input edge length.
    bbox : (xmin, ymin, xmax, ymax) or None
        Bounding box. If None, inferred from input points with 5% margin.
    dptol, ttol, Fscale, deltat, geps_factor, niter, seed : truss-loop tunables
        Forwarded verbatim to the vendored distmesh2d_warmstart. See R3 in
        research.md for descriptions.
    sdf_tolerance : float
        Tolerance for boundary-on-SDF check. Default 1e-6.
    enforce_non_degradation : bool
        If True (default), output median quality < input median quality causes
        the input to be returned unchanged plus a RuntimeWarning. If False,
        always return the truss-loop output regardless of quality.

    Returns
    -------
    points_out : (Nout, 2) ndarray
        Optimized vertex coordinates. Boundary points (the ones at indices
        `boundary_indices` in the input) appear at indices 0..B-1 of the output
        in the same order, byte-equal to the input boundary subset.
    triangles_out : (Mout, 3) ndarray of int
        Final triangulation from the truss loop's last Delaunay rebuild.

    Raises
    ------
    NotImplementedError
        If `triangles.shape[1] != 3`.
    ValueError
        - Input triangles have non-positive signed area.
        - Boundary points fail the SDF-zero-set check (|sdf(p)| < sdf_tolerance).
        - Input points fall outside the supplied/inferred bbox.
    ImportError
        If ADMESH cannot be imported. The error message names the missing
        module and references the pinned commit SHA from FR-015.
    RuntimeWarning (emitted, not raised)
        - Truss loop hit niter cap without convergence (returns best-so-far).
        - Output quality regressed and enforce_non_degradation=True (returns input).
    """
```

### Validation order (deterministic, first failure raises)

1. `NotImplementedError` if `triangles.shape[1] != 3`.
2. `ValueError` if any triangle has signed area ≤ 0.
3. Compute `boundary_indices` if None.
4. `ValueError` if any boundary point has `|sdf(p)| ≥ sdf_tolerance`. Error message: `"boundary not on SDF zero set: {n} of {B} boundary points fail; max |sdf| = {value:.2e} at index {idx}; tolerance = {sdf_tolerance:.2e}"`.
5. Compute `bbox` if None; verify all input points lie within it; `ValueError` if not.
6. Compute `h0` if None.
7. (warn-only) `RuntimeWarning` if any interior point has `sdf > 0`.
8. (warn-only) `RuntimeWarning` if `len(interior_indices) < 10` (per spec edge case).

### Behavioral guarantees

| Guarantee | Verifiable how |
|-----------|---------------|
| Input boundary is bit-exact preserved | `np.array_equal(points_out[:B], input_points[boundary_indices])` |
| Input is never mutated | Caller can pass read-only arrays (`flags.writeable = False`) and call succeeds |
| RNG is deterministic for fixed `seed` | Two calls with identical inputs and identical `seed` produce identical outputs |
| Non-degradation when enforced | `median_q(out) >= median_q(in)` always holds |
| ADMESH version is pinned | Function imports from `_vendor_admesh_truss`; vendored loop's header docstring lists the SHA |

---

## High-Level: `optimize_with_admesh_truss`

### Signature

```python
def optimize_with_admesh_truss(
    mesh: CHILmesh,
    sdf: Callable[[np.ndarray], np.ndarray],
    size_fn: Optional[Callable[[np.ndarray], np.ndarray]] = None,
    **kwargs,
) -> CHILmesh:
    """
    CHILmesh-native wrapper around `optimize_with_admesh_truss_arrays`.

    Extracts (points, triangles, boundary_indices) from `mesh`, calls the
    array form, and rewraps the result in a fresh CHILmesh.

    All kwargs are forwarded to the array form. See its docstring for the
    full parameter list and behavior.

    Parameters
    ----------
    mesh : CHILmesh
        Input mesh. Must be triangle-only. Must have a non-empty boundary
        (i.e., be an open / non-closed manifold).
    sdf, size_fn, **kwargs
        See `optimize_with_admesh_truss_arrays`.

    Returns
    -------
    CHILmesh
        Fresh CHILmesh instance. Adjacencies and layers are recomputed from
        scratch. The input mesh is not modified.

    Raises
    ------
    All exceptions from `optimize_with_admesh_truss_arrays`, plus:
    NotImplementedError
        If `mesh` is mixed-element (has both triangles and quads). The error
        message says "high-level wrapper requires triangle-only CHILmesh; use
        the array form with manual filtering for mixed meshes".
    """
```

### Boundary identification (CHILmesh form)

```python
boundary_edge_ids = mesh.boundary_edges()
edge2vert = mesh.adjacencies['Edge2Vert']
boundary_indices = np.unique(edge2vert[boundary_edge_ids].flatten())
```

This is the same pattern used in spec 004's `generate_4row_admesh.py` and is the canonical CHILmesh way to identify boundary vertices.

---

## Error Catalog

| Exception | Trigger | Message format |
|-----------|---------|---------------|
| `NotImplementedError` | non-triangle input | `"warm-start truss requires triangle-only mesh; got {shape[1]}-vertex elements"` |
| `ValueError` (V_AREA) | degenerate triangles | `"input has {n} triangles with non-positive signed area: indices {list[:10]}{...}"` |
| `ValueError` (V_BND_SDF) | boundary not on SDF | `"boundary not on SDF zero set: {n} of {B} boundary points fail; max \|sdf\| = {v:.2e} at index {idx}; tolerance = {sdf_tolerance:.2e}"` |
| `ValueError` (V_BBOX) | points outside bbox | `"input has {n} points outside bbox {bbox}: indices {list[:10]}{...}"` |
| `ImportError` | ADMESH missing | `"ADMESH (admesh.distmesh) is required for warm-start truss; pinned to SHA 05bc68f. Install via 'pip install -e /path/to/admesh' or 'pip install admesh2D'."` |
| `RuntimeWarning` | non-convergence | `"distmesh truss did not converge in {niter} iterations; returning best-so-far mesh"` |
| `RuntimeWarning` | quality regression | `"warm-start truss did not improve quality (input median={q_in:.4f}, output median={q_out:.4f}); returning input unchanged (set enforce_non_degradation=False to override)"` |
| `RuntimeWarning` | small interior set | `"warm-start truss called with only {n} interior points (threshold: 10); optimization may have minimal effect"` |
| `RuntimeWarning` | interior outside SDF | `"{n} interior points have sdf > 0 (outside domain); truss may push them onto boundary or stall"` |

---

## Determinism Contract

```python
out1 = optimize_with_admesh_truss_arrays(p, t, sdf, size_fn, seed=42)
out2 = optimize_with_admesh_truss_arrays(p, t, sdf, size_fn, seed=42)
assert np.array_equal(out1[0], out2[0])
assert np.array_equal(out1[1], out2[1])
```

Determinism holds when: callables are pure, numpy/scipy versions unchanged, pinned ADMESH commit unchanged. Default `seed=0` ensures reproducible demo PNG.

---

## Versioning & Stability

New feature, not modifying existing API. Contract is **public and stable**. Breaking change requires: new spec/plan/tasks pair, `DeprecationWarning` cycle for at least one minor version.

`_vendor_admesh_truss` is **private**. Callers MUST NOT import from it. Deleted without deprecation when ADMESH-B lands.
