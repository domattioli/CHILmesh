# Data Model: ADMESH-Domains I/O Parity

**Phase 1 output for `/speckit-plan`**
**Branch**: `planning-optimize_modernize` | **Date**: 2026-04-26

---

## Core Entities

### `CHILmesh` (extended)

The single in-memory representation of a loaded mesh. All changes in this feature are additive extensions.

| Attribute | Type | Description | Change |
|---|---|---|---|
| `connectivity_list` | `ndarray (n_elems, 3\|4)` | Element connectivity, 0-based. Triangles in 4-col array use pad: `[v0,v1,v2,v0]` | Extended to support 4-col |
| `points` | `ndarray (n_verts, 3)` | Node coordinates `[x, y, z]` | Unchanged |
| `n_verts` | `int` | Node count | Unchanged |
| `n_elems` | `int` | Element count | Unchanged |
| `n_layers` | `int` | Skeletonization depth; `0` when `compute_layers=False` | Guard added to `get_layer()` |
| `type` | `str \| None` | Element type: `"Triangular"`, `"Quadrilateral"`, or `"Mixed-Element"` | Now set in `_initialize_mesh()` |
| `layers` | `dict` | OE, IE, OV, IV, bEdgeIDs lists | Unchanged (empty when `compute_layers=False`) |
| `grid_name` | `str \| None` | Mesh title from file header | Unchanged |

**New `__init__` signature**:
```python
def __init__(
    self,
    connectivity: np.ndarray | None = None,
    points: np.ndarray | None = None,
    grid_name: str | None = None,
    compute_layers: bool = True,       # NEW
) -> None
```

---

### Metadata Dict

Returned by `admesh_metadata()`. Format-agnostic; compatible with ADMESH-Domains schema field names (note: ADMESH-Domains maps `min_x→min_lon`, `min_y→min_lat` in its own adapter).

```python
{
    "node_count":     int,           # == mesh.n_verts
    "element_count":  int,           # == mesh.n_elems
    "element_type":   str,           # "Triangular" | "Quadrilateral" | "Mixed-Element"
    "bounding_box": {
        "min_x": float,
        "max_x": float,
        "min_y": float,
        "max_y": float,
    }
}
```

---

### ADMESH-Domains Record (duck-typed, external)

Not owned by CHILmesh. `from_admesh_domain()` reads these attributes via `getattr` — no import of `admesh-domains` required.

| Attribute | Type | Usage |
|---|---|---|
| `filename` | `str` | Path to mesh file on disk |
| `type` | `str \| None` | Format hint: `"ADCIRC"`, `"SMS_2DM"`, `"ADCIRC_GRD"` |
| `kind` | `str \| None` | Optional: `"boundary"` suppresses layer warning |

---

## State Transitions

### `CHILmesh` initialisation lifecycle

```
__init__(compute_layers=True)
  └── _initialize_mesh(compute_layers=True)
        ├── n_verts, n_elems ← set from arrays
        ├── _ensure_ccw_orientation()
        ├── _build_adjacencies()       ← skipped when compute_layers=False
        ├── _skeletonize()             ← skipped when compute_layers=False
        └── type ← set from _elem_type() result    [NEW]

get_layer(idx)
  ├── n_layers == 0 → RuntimeError("Layers not computed...")  [NEW guard]
  └── idx out of range → ValueError (existing)
```

### `from_admesh_domain()` routing

```
from_admesh_domain(record, compute_layers=True)
  ├── Path(record.filename).exists() → False → FileNotFoundError
  ├── record.type == "SMS_2DM" → read_from_2dm(filename)
  ├── record.type in {None, "ADCIRC", "ADCIRC_GRD"} → read_from_fort14(filename)
  └── record.type is other → warnings.warn(...) → read_from_fort14(filename)
  [then pass compute_layers kwarg to CHILmesh.__init__]
```

---

## Validation Rules

| Rule | Enforcement |
|---|---|
| `num_nodes` must be 3 or 4 in `.fort.14` | `read_from_fort14()` raises `ValueError` with element ID |
| `num_nodes` must be 3 or 4 in `.2dm` | `read_from_2dm()` raises `ValueError` with element ID |
| `get_layer()` requires `compute_layers=True` at init time | `RuntimeError` raised in `get_layer()` when `n_layers == 0` |
| File must exist before any reader is called | `FileNotFoundError` raised in `from_admesh_domain()` |
| `element_type` vocabulary is constrained | `admesh_metadata()` always returns one of the three canonical strings |

---

## Connectivity Convention (padded triangles)

For a mixed-element mesh stored in a 4-column `connectivity_list`:

- **Quad**: `[v0, v1, v2, v3]` — all four distinct vertices
- **Triangle**: `[v0, v1, v2, v0]` — fourth slot repeats first vertex

This convention is already established in `_ensure_ccw_orientation()` (line 140) and `_elem_type()`. All new code must honour it.

`_elem_type()` identifies a row as a triangle when `row[3] == row[0]` (or the array has only 3 columns).
