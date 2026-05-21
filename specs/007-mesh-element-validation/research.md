# Research: Mesh Element Validity Test Suite

## Geometric predicates

### Segment proper-crossing

Two open segments `(p, q)` and `(r, s)` cross properly iff `p, q` lie on strict opposite sides of line `rs` AND `r, s` lie on strict opposite sides of line `pq`. "Strict" means signed-area magnitude `> tol_effective`; equal-to-tol counts as collinear/touching → NO proper crossing (so shared endpoints don't trigger).

```python
def _orient(a, b, c, tol):
    cross = (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])
    if cross >  tol: return  1
    if cross < -tol: return -1
    return 0

def segment_proper_cross(p, q, r, s, tol):
    o1 = _orient(p, q, r, tol)
    o2 = _orient(p, q, s, tol)
    o3 = _orient(r, s, p, tol)
    o4 = _orient(r, s, q, tol)
    return o1 != 0 and o2 != 0 and o3 != 0 and o4 != 0 and o1 != o2 and o3 != o4
```

Vectorized variant operates on `(N, 2, 2)` segment arrays for bulk pair checks.

### Bowtie quad

Quad `v0, v1, v2, v3` is self-intersecting iff:
- edge `(v0, v1)` proper-crosses edge `(v2, v3)` OR
- edge `(v1, v2)` proper-crosses edge `(v3, v0)`

(These are the two "non-adjacent" edge pairs of the quad. The two adjacent pairs share a vertex and can't proper-cross by definition.)

### Point in polygon (winding number)

For a point `p` and polygon `[v0, ..., v(n-1)]`, compute winding number by counting signed crossings of the horizontal ray from `p`. Strict interior (winding ≠ 0 AND point not on any edge within `tol`). On-edge → "outside" (so shared-vertex / shared-edge cases don't trigger interior-overlap false positives).

## Broadphase: uniform grid hash

```text
cell_size ≈ 2 × mean(element bbox diag) = 2 × √(mesh_area / n_elems)
n_cells_x = ceil((bbox_max_x - bbox_min_x) / cell_size)
```

Each element registers in every cell its AABB overlaps. Pair generation: for each cell, all C(k, 2) pairs where k = elements in cell; dedup via `frozenset({i, j})`.

Expected: on `block_o` (n_elems ≈ 5200, roughly square domain), cell_size produces ~1200 cells, ~4-5 elements/cell, → ~25k unique candidate pairs vs `O(n²) ≈ 2.7e7`. **>1000× reduction** (SC-007 ≥100× passes with margin).

Fallback for `n_elems ≤ 5000`: `O(n²)` direct (per FR-012). Annulus / donut / structured fit this; only block_o needs the broadphase.

## Element-type classification

```python
def classify_element(row: np.ndarray) -> str:
    a, b, c, d = row
    if d == -1 or d == a or d == b or d == c:
        # padded triangle or duplicate-fourth-vertex degenerate triangle
        return "TRI"
    if len({a, b, c, d}) < 4:
        return "DEGENERATE_QUAD"   # duplicate-but-non-padding → still QUAD shape, noted
    return "QUAD"
```

`DEGENERATE_QUAD` classifier rolls into `QUAD` for the per-element validity tests but emits a `DEGENERATE_QUAD_DUPLICATE_VERTEX` informational note (FR-004).

## Edge sharing map

```python
edge_to_elems: dict[frozenset[int], list[int]] = {}
for elem_id, row in enumerate(connectivity_list):
    verts = [v for v in row if v != -1]
    n = len(verts)
    for i in range(n):
        key = frozenset({verts[i], verts[(i+1) % n]})
        edge_to_elems.setdefault(key, []).append(elem_id)
```

Two elements share an edge iff they appear together in any `edge_to_elems` value. Edge-sharing pairs are exempt from the `EDGE_CROSSING` check.

## Constructor bypass for synthetic fixtures

`CHILmesh.__init__` calls `_degeneracy_fallback()` (per `tests/test_degeneracy.py`) which silently repairs some malformed connectivity. To plant bowtie / overlap / pentagon fixtures, build a valid mesh first, then:

```python
def corrupt_to(mesh: CHILmesh, row_id: int, new_row: list[int]) -> CHILmesh:
    mesh.connectivity_list[row_id] = new_row
    # Invalidate caches that depend on connectivity:
    if hasattr(mesh, "_layers_cache"):
        mesh._layers_cache = None
    if hasattr(mesh, "adjacencies"):
        mesh.adjacencies = None
    return mesh
```

(Exact cache-key names verified during Phase 3 against `src/chilmesh/CHILmesh.py`.)

## Existing-test cross-check

`tests/test_degeneracy.py` asserts that padded-triangle and duplicate-vertex quads load successfully. The new validator MUST agree (Q3, FR-014): degenerate ≠ violation. Verified by parametrizing the new suite over the same fixtures used by `test_degeneracy.py` and asserting `report.ok == True`.

`tests/test_invariants.py:test_layers_disjoint_cover` already enforces that every element appears in some layer. The new FR-007 layer-membership check piggybacks on this invariant.

## Tolerance scaling reference

| Coord system | Typical bbox_diag | `tol_effective = 1e-12 * bbox_diag` |
|--------------|-------------------|------------------------------------|
| meter-scale Cartesian | 1e4 m | 1e-8 m |
| lat/lon (degrees)     | 1e1 ° | 1e-11 ° |
| UTM (meters)          | 1e5 m | 1e-7 m |
| km-scale (e.g. global) | 1e7 m | 1e-5 m |

Floor at `1e-15` guards against degenerate empty bbox.

## Open items deferred to plan-execution

- Exact `mesh._layers_cache` / `mesh.adjacencies` attribute names verified at implementation time by reading `src/chilmesh/CHILmesh.py`.
- Whether `_invalidate_cache()` is a public helper or needs to be inlined.
