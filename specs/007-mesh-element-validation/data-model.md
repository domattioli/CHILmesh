# Data Model: Mesh Element Validity

No persistent storage. Three frozen dataclasses in `tests/_validity/types.py`:

## `Violation`

| Field | Type | Notes |
|-------|------|-------|
| `category` | `str` | One of the FR-defined error codes (see spec.md violation table) |
| `element_ids` | `tuple[int, ...]` | Offending element ID(s). Singleton for per-element rules; pair for pair-wise rules |
| `edge_ids` | `tuple[int, ...] \| None` | Offending edge ID(s) when applicable (FR-009, FR-011 EDGE_CROSSING) |
| `detail` | `str` | One-line geometric specifics (vertex coords, edge IDs, cross point) |

Hashable (frozen). JSON-serializable via `dataclasses.asdict()`.

## `InformationalNote`

| Field | Type | Notes |
|-------|------|-------|
| `category` | `str` | `DEGENERATE_QUAD_DUPLICATE_VERTEX`, `DEGENERATE_ZERO_AREA`, `EMPTY_MESH`, `LAYERS_AUTO_TRIGGERED` |
| `element_ids` | `tuple[int, ...]` | Affected elements (empty for mesh-wide notes) |
| `detail` | `str` | Free-form description |

Same shape as `Violation` minus `edge_ids`. Notes do NOT flip `ok` to False.

## `MeshValidityReport`

| Field | Type | Notes |
|-------|------|-------|
| `ok` | `bool` | True iff `len(violations) == 0` |
| `violations` | `tuple[Violation, ...]` | All violations found (no per-category cap; FR-016 caps only the pytest error message) |
| `notes` | `tuple[InformationalNote, ...]` | Informational records |
| `n_elems_checked` | `int` | Equals `mesh.n_elems` after broadphase deduping |
| `runtime_s` | `float` | Wall-clock from entry to return |

## Helper signatures (predicates.py)

```python
def bbox_diag(points: np.ndarray) -> float: ...
def classify_element(row: np.ndarray) -> str: ...  # returns 'TRI' | 'QUAD' | 'PENTAGON_OR_ABOVE'
def segment_proper_cross(p, q, r, s, tol: float) -> bool: ...
def point_strictly_in_polygon(p: np.ndarray, polygon: np.ndarray, tol: float) -> bool: ...
def is_self_intersecting_quad(quad_pts: np.ndarray, tol: float) -> bool: ...
```

## Broadphase (broadphase.py)

```python
class UniformGridIndex:
    def __init__(self, elem_bboxes: np.ndarray, cell_size: float, bbox_min: np.ndarray): ...
    def candidate_pairs(self) -> Iterator[tuple[int, int]]: ...
```

`elem_bboxes` shape `(n_elems, 4)` columns `[min_x, min_y, max_x, max_y]`. Yields each unique pair exactly once.

## Fixtures (fixtures.py)

| Fixture | Construction | Planted violation |
|---------|--------------|-------------------|
| `bowtie_quad_mesh()` | `chilmesh.examples.structured()` minimal 2x2 → `corrupt_to(mesh, 0, [v0,v1,v3,v2])` | `SELF_INTERSECTING_QUAD` |
| `interior_triangle_mesh()` | Build 3x3 quad grid → `corrupt_to(mesh, 4, [v_center, v_a, v_b, -1])` | `INTERIOR_TRIANGLE_FORBIDDEN` |
| `pentagon_mesh()` | NOT representable in 4-column connectivity → simulate by spying a faked row in a fake CHILmesh subclass that exposes 5-column connectivity (or skip with `xfail`) | `UNSUPPORTED_ELEMENT_ARITY` |
| `overlapping_quads_mesh()` | Two valid quads, then translate one to overlap the other but share no vertices | `INTERIOR_OVERLAP` |
| `edge_crossing_mesh()` | Two valid quads in different topological cells but with edges that cross geometrically (no shared vertex) | `EDGE_CROSSING` |

For `pentagon_mesh`, the spec FR-003 requires the validator to catch arity violations even if CHILmesh itself can't construct one. We expose a `_FakeMesh` test double that satisfies the duck-typed interface (`connectivity_list`, `points`, `n_elems`, `boundary_vertices`, `layers`) and ships a 5-column row. This is the only place a test double is used; everywhere else uses a real `CHILmesh`.
