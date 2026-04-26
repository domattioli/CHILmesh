# Research: ADMESH-Domains I/O Parity

**Phase 0 output for `/speckit-plan`**
**Branch**: `planning-optimize_modernize` | **Date**: 2026-04-26

---

## 1. Reader fix — quad/mixed elements in `.fort.14`

**Question**: What is the minimal change to `read_from_fort14()` that handles 3- and 4-node elements without breaking existing triangular-only meshes?

**Finding** (`src/chilmesh/CHILmesh.py:688–698`):
The reader pre-allocates `elements = np.zeros((n_elements, 3), ...)` then raises on `num_nodes != 3`. For mixed-element support:
- Pre-scan the element block to determine max `num_nodes` (1 pass), or
- Use a Python list and convert after reading (simpler, one pass).

**Decision**: Two-pass approach — scan for max node count first, then allocate. Cleaner than a list because `connectivity_list` is always a dense NumPy array downstream.

**Padded-triangle convention** (already used in CCW code at line 140): triangles in a 4-column array store `[v0, v1, v2, v0]` (repeat first vertex). This matches the project's existing convention for mixed-element meshes.

**Alternatives considered**:
- Ragged Python lists — rejected; all downstream code (adjacency, skeletonize) assumes dense ndarray.
- Variable-width dtype — rejected; unsupported by NumPy dense arrays.

---

## 2. Writer fix — quad/mixed elements in `write_fort14()`

**Question**: What does a valid quad line look like in ADCIRC `.fort.14` format?

**Finding**: The ADCIRC grid format uses `element_id num_nodes v1 v2 v3 [v4]`. For a triangle: `1 3 n1 n2 n3`. For a quad: `1 4 n1 n2 n3 n4`. The `num_nodes` field drives the count.

**Current bug** (`src/chilmesh/CHILmesh.py:989–991`):
```python
for i, tri in enumerate(elements, start=1):
    n1, n2, n3 = tri + 1   # crashes for quads
    f.write(f"{i} 3 {n1} {n2} {n3}\n")  # hardcoded 3
```

**Decision**: Detect element type per row using `_elem_type()` — or more efficiently, check whether the 4th column equals the 1st column (pad check). Write 3 or 4 node indices accordingly. This is atomic with the reader fix (same PR / commit).

---

## 3. `compute_layers=False` kwarg — fast initialisation

**Question**: Where in `_initialize_mesh()` does the slow path live, and what is the minimal change?

**Finding** (`src/chilmesh/CHILmesh.py:101–118`): `_initialize_mesh()` always calls `_skeletonize()` at line 118, which drives the ~30s Block_O cost. `_build_adjacencies()` at line 115 is also O(n²) but needed for `admesh_metadata()` — actually no: `admesh_metadata()` only needs `n_verts`, `n_elems`, `type`, and `points` extrema, none of which require adjacency computation.

**Decision**: Add `compute_layers: bool = True` to `__init__` and `_initialize_mesh`. When `False`, skip both `_build_adjacencies()` and `_skeletonize()`. Gate `get_layer()` with a clear error when `n_layers == 0`. This achieves <2s for all fixture sizes.

**Guard message** (from clarify phase F-4):
```
"Layers not computed. Re-initialise with compute_layers=True."
```

---

## 4. `admesh_metadata()` method

**Question**: Which fields and key names does ADMESH-Domains need?

**Finding** (`.planning/ADMESH-DOMAINS-CLARIFY.md`, Finding C-3 and Q-2):
- ADMESH-Domains schema uses `min_lat/max_lat/min_lon/max_lon` but CHILmesh fixtures are Cartesian.
- **Decision** (Q-2 Option A): CHILmesh returns `min_x/max_x/min_y/max_y`. The adapter layer in ADMESH-Domains maps x→lon, y→lat.
- `element_type` vocabulary: `"Triangular"` | `"Quadrilateral"` | `"Mixed-Element"` (matches `_elem_type()` output logic).

**Method signature**:
```python
def admesh_metadata(self) -> dict:
    # Returns: node_count, element_count, element_type, bounding_box
```

`element_type` is derived by calling `_elem_type()` on all elements and checking counts — if all tri → "Triangular", all quad → "Quadrilateral", mixed → "Mixed-Element".

**`type` property** (set at mesh load time): The `CHILmesh.type` attribute is already declared in `__init__` (line 81) as `None`. Set it to the `element_type` string during `_initialize_mesh()`.

---

## 5. `from_admesh_domain()` class method — routing

**Question**: How should `from_admesh_domain()` route between readers?

**Finding** (`.planning/ADMESH-DOMAINS-CLARIFY.md`, Q-3 and C-4):
- Duck-typed record: `getattr(record, "filename", None)` and `getattr(record, "type", None)`.
- Route `"SMS_2DM"` → `read_from_2dm()`, everything else → `read_from_fort14()` (covers `"ADCIRC"` and `"ADCIRC_GRD"`).
- `compute_layers` kwarg forwarded through.
- File-not-found: check `Path(filename).exists()` before opening; raise `FileNotFoundError` with guidance.
- Unknown type: fall back to ADCIRC reader with `warnings.warn()`.

**Duck-typed attributes used**:
- `record.filename` (str) — path on disk
- `record.type` (str, optional) — format hint
- `record.kind` (str, optional) — `"boundary"` suppresses layer warning

---

## 6. SMS `.2dm` reader

**Question**: What is the SMS `.2dm` file format?

**Finding**: SMS `.2dm` is a text format with keyword-prefixed lines:
- `MESH2D` — header (first line, optional)
- `ND id x y z` — node definition
- `E3T id n1 n2 n3 mat` — triangular element (mat = material ID, ignored)
- `E4Q id n1 n2 n3 n4 mat` — quad element
- Lines starting with `#` are comments

**Decision**: Parse node and element lines in a single pass using a keyword dispatch dict. Nodes are 1-indexed in the file; subtract 1 for 0-based storage. Material ID is discarded. Mixed meshes (E3T + E4Q in same file) use padded-triangle convention.

**Alternatives considered**: Using a library (meshio) — rejected; adds a heavy dependency for a format with a trivial parser.

---

## 7. Error handling strategy

**Question**: What error types and messages are most actionable for each failure mode?

**Decision** (from clarify phase F-3, F-4, spec FR-007):

| Scenario | Exception | Message |
|---|---|---|
| File not found | `FileNotFoundError` | `"File not found: {path}. If using ADMESH-Domains, call mesh_record.load() first."` |
| Unsupported element type (5+ nodes) | `ValueError` | `"Unsupported element type: element {id} has {n} nodes (only 3 or 4 supported)."` |
| `get_layer()` with no layers | `RuntimeError` | `"Layers not computed. Re-initialise with compute_layers=True."` |
| Unknown record type in `from_admesh_domain` | `warnings.warn` (not error) | `"Unrecognised mesh type '{type}', falling back to ADCIRC reader."` |

---

## 8. Test fixture — synthetic quad mesh

**Question**: What is the minimal quad `.fort.14` fixture that tests all code paths?

**Decision** (from clarify phase Q-1, Option C):
- Synthetic: 2×2 regular quad grid → 9 nodes, 4 quad elements. Committed to `src/chilmesh/data/quad_2x2.fort.14`.
- Deterministic, no download required, covers reader + writer + roundtrip + metadata tests.
- Real ADMESH-Domains mesh: deferred to a future `@pytest.mark.integration` test (separate issue).

**Format** (9 nodes, 4 quads, all 1-based in file):
```
2x2 quad test mesh
4 9
1  0.0  0.0  0.0
2  1.0  0.0  0.0
...
1 4 1 2 5 4
2 4 2 3 6 5
...
```

---

## Summary of decisions

| # | Decision | Rationale |
|---|---|---|
| 1 | Two-pass reader (scan then allocate) | Single dense ndarray, clean |
| 2 | Per-row node count in writer | Handles tri/quad/mixed uniformly |
| 3 | `compute_layers=False` skips adjacency + skeletonize | <2s target achievable |
| 4 | `bounding_box` uses `min_x/max_x/min_y/max_y` | Coordinate-agnostic |
| 5 | Duck-typed `from_admesh_domain` | No import of admesh-domains |
| 6 | Simple keyword-dispatch `.2dm` parser | Avoids heavy dependency |
| 7 | Structured error table | Each error is actionable |
| 8 | Synthetic 2×2 quad fixture | Fast unit tests, no download |
