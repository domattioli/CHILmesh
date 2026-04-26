# Public API Contracts: ADMESH-Domains I/O Parity

**Phase 1 output for `/speckit-plan`**
**Branch**: `planning-optimize_modernize` | **Date**: 2026-04-26

These are the contracts between CHILmesh and its callers (researchers, ADMESH-Domains, MADMESHR). Each entry defines the function signature, pre/post-conditions, and error behaviour. Implementation must not deviate from these contracts without updating this document.

---

## `CHILmesh.__init__` (extended)

```python
CHILmesh(
    connectivity: np.ndarray | None = None,
    points: np.ndarray | None = None,
    grid_name: str | None = None,
    compute_layers: bool = True,
) -> CHILmesh
```

**Pre-conditions**:
- `connectivity` shape: `(n_elems, 3)` or `(n_elems, 4)`, dtype int, 0-based indices
- `points` shape: `(n_verts, 2)` or `(n_verts, 3)`, dtype float
- `compute_layers=False` is valid for any mesh size/type

**Post-conditions**:
- `mesh.n_verts == len(points)`
- `mesh.n_elems == len(connectivity)`
- `mesh.type` is one of `"Triangular"`, `"Quadrilateral"`, `"Mixed-Element"`, or `None` if no connectivity
- When `compute_layers=False`: `mesh.n_layers == 0`, `mesh.adjacencies == {}`
- When `compute_layers=True`: `mesh.n_layers >= 1`, adjacencies fully populated

**Errors**: none beyond pre-condition violations (handled by callers)

---

## `CHILmesh.read_from_fort14` (extended)

```python
@staticmethod
CHILmesh.read_from_fort14(
    full_file_name: Path,
    compute_layers: bool = True,
) -> CHILmesh
```

**Pre-conditions**:
- `full_file_name` points to a readable ADCIRC `.fort.14` file
- File elements have `num_nodes` of 3 or 4

**Post-conditions**:
- Returns initialised `CHILmesh` with all nodes and elements from the file
- Round-trip fidelity: `n_verts` and `n_elems` match file header counts
- `compute_layers` forwarded to `CHILmesh.__init__`

**Errors**:
- `FileNotFoundError` — file does not exist
- `ValueError(f"Unsupported element type: element {id} has {n} nodes (only 3 or 4 supported).")` — element with 5+ nodes

---

## `CHILmesh.read_from_2dm` (new)

```python
@staticmethod
CHILmesh.read_from_2dm(
    full_file_name: Path,
    compute_layers: bool = True,
) -> CHILmesh
```

**Pre-conditions**:
- `full_file_name` points to a readable SMS `.2dm` file
- Node lines: `ND id x y z`
- Element lines: `E3T id n1 n2 n3 mat` or `E4Q id n1 n2 n3 n4 mat`

**Post-conditions**:
- Returns initialised `CHILmesh`; all geometry preserved
- Mixed meshes (E3T + E4Q) returned as 4-column padded-triangle connectivity
- `compute_layers` forwarded to `CHILmesh.__init__`

**Errors**:
- `FileNotFoundError` — file does not exist
- `ValueError(f"Unsupported element type: element {id} has {n} nodes.")` — element keyword other than E3T/E4Q

---

## `CHILmesh.from_admesh_domain` (new)

```python
@classmethod
CHILmesh.from_admesh_domain(
    record,                        # duck-typed ADMESH-Domains Mesh record
    compute_layers: bool = True,
) -> CHILmesh
```

**Pre-conditions**:
- `record.filename` (str) — path to mesh file on disk
- `record.type` (str, optional) — one of `"ADCIRC"`, `"SMS_2DM"`, `"ADCIRC_GRD"`, or absent/None
- File at `record.filename` exists

**Post-conditions**:
- Returns correctly initialised `CHILmesh`
- `"SMS_2DM"` routes to `read_from_2dm()`; all other types route to `read_from_fort14()`
- Unknown `record.type` emits `warnings.warn(...)` then falls back to `read_from_fort14()`

**Errors**:
- `FileNotFoundError("File not found: {path}. If using ADMESH-Domains, call mesh_record.load() first.")` — file absent

---

## `CHILmesh.admesh_metadata` (new)

```python
def admesh_metadata(self) -> dict
```

**Pre-conditions**: mesh must be initialised (`n_verts > 0`)

**Post-conditions**:
```python
{
    "node_count":    int,            # equals self.n_verts
    "element_count": int,            # equals self.n_elems
    "element_type":  str,            # "Triangular" | "Quadrilateral" | "Mixed-Element"
    "bounding_box": {
        "min_x": float,              # min of points[:, 0]
        "max_x": float,              # max of points[:, 0]
        "min_y": float,              # min of points[:, 1]
        "max_y": float,              # max of points[:, 1]
    }
}
```

**Errors**: none (safe to call with `compute_layers=False`)

---

## `CHILmesh.get_layer` (guard added)

```python
def get_layer(self, layer_idx: int) -> dict
```

**Pre-conditions**: `compute_layers=True` was used at initialisation

**Post-conditions**: unchanged from existing contract

**Errors** (extended):
- `RuntimeError("Layers not computed. Re-initialise with compute_layers=True.")` — when `self.n_layers == 0`
- `ValueError(f"Layer index {layer_idx} out of range [0, {self.n_layers-1}]")` — existing, unchanged

---

## `write_fort14` (module-level, fixed)

```python
write_fort14(
    filename: Path,
    points: np.ndarray,
    elements: np.ndarray,   # (n_elems, 3) or (n_elems, 4) with padded triangles
    grid_name: str,
) -> bool
```

**Pre-conditions**:
- `elements` may have 3 or 4 columns
- 4-column triangles follow padded convention (`row[3] == row[0]`)

**Post-conditions**:
- Written file is valid ADCIRC `.fort.14`
- Triangle rows write `"i 3 n1 n2 n3"`, quad rows write `"i 4 n1 n2 n3 n4"`
- Round-trip: reload of written file produces same `n_verts` and `n_elems`

**Errors**: returns `False` on I/O exception (existing behaviour, unchanged)

---

## Synthetic quad test fixture

**File**: `src/chilmesh/data/quad_2x2.fort.14`

Format:
```
2x2 quad test mesh
4 9
1 0.00000000 0.00000000 0.00000000
2 1.00000000 0.00000000 0.00000000
3 2.00000000 0.00000000 0.00000000
4 0.00000000 1.00000000 0.00000000
5 1.00000000 1.00000000 0.00000000
6 2.00000000 1.00000000 0.00000000
7 0.00000000 2.00000000 0.00000000
8 1.00000000 2.00000000 0.00000000
9 2.00000000 2.00000000 0.00000000
1 4 1 2 5 4
2 4 2 3 6 5
3 4 4 5 8 7
4 4 5 6 9 8
```

9 nodes, 4 quad elements. `bounding_box`: `{min_x:0, max_x:2, min_y:0, max_y:2}`. `element_type`: `"Quadrilateral"`.
