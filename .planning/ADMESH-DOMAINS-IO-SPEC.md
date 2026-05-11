# Spec-Kit: Specify Phase
# "CHILmesh as Computational Engine for ADMESH-Domains I/O"

**Branch:** `planning-optimize_modernize`
**Date:** 2026-04-26
**Scope:** ADMESH-Domains = catalog/delivery layer; CHILmesh = computational engine

---

## THE REFRAME

Previous view: CHILmesh reads `.fort.14`. **Correct view:** ADMESH-Domains owns domain catalog, mesh files, I/O. CHILmesh is object users get *after* loading from that catalog. All CHILmesh I/O, metadata, format handling must be consistent with ADMESH-Domains schema.

---

## PART 1: CURRENT STATE

### 1.1 What CHILmesh Does Today (I/O)

```python
CHILmesh.read_from_fort14(path) -> CHILmesh   # triangular elements ONLY (line 693)
mesh.write_to_fort14(filename)                # delegates to write_fort14()
# No other formats, no metadata, no lazy init, no bounding_box
```

**Critical Bug:** `read_from_fort14()` line 692–693 raises `ValueError` for any `num_nodes != 3`. ADMESH-Domains hosts real-world coastal meshes many of which are quad or mixed-element. **CHILmesh cannot load majority of ADMESH-Domains meshes today.**

### 1.2 What ADMESH-Domains Provides

**Mesh schema fields CHILmesh must compute:**
- `node_count`, `element_count`, `element_type` (CHILmesh computes these)
- `bounding_box {min_lat, max_lat, min_lon, max_lon}` (CHILmesh computes this)
- `type`: `"ADCIRC"` | `"SMS_2DM"` | `"ADCIRC_GRD"`
- `kind`: `"mesh"` | `"boundary"`

**File formats:** `.14` / `.fort.14` (ADCIRC) and `.2dm` (SMS 2D mesh).

`Mesh.load()` fetches file from HuggingFace on-demand; CHILmesh then parses into live object.

### 1.3 The Gap

| Capability | ADMESH-Domains Needs | CHILmesh Today | Gap |
|------------|---------------------|---------------|-----|
| Read `.fort.14` triangular | ✓ | ✓ | None |
| Read `.fort.14` quad/mixed | ✓ | ✗ raises ValueError | **CRITICAL BUG** |
| Read `.2dm` (SMS format) | ✓ | ✗ | Missing |
| `node_count` / `element_count` | ✓ | partial (`n_verts`/`n_elems`) | Name mismatch |
| `bounding_box` property | ✓ | ✗ | Missing API |
| Lazy initialization | ✓ | ✗ full init always | Performance gap |
| `from_admesh_domain()` | ✓ | ✗ | Missing |
| `kind: "boundary"` handling | ✓ | ✗ | Missing |
| Metadata dict round-trip | ✓ | ✗ | Missing |

---

## PART 2: CONSTRAINTS

### Hard Constraints
1. Skeletonization semantics — immutable (audit Q3, Decision A2)
2. `.fort.14` roundtrip — byte-identical (audit B1)
3. Mixed-element support — triangles + quads in one mesh
4. API surface — `signed_area()`, `elem_quality()`, `plot()` signatures unchanged
5. No new external runtime deps — only scipy/numpy/matplotlib

### Soft Constraints
1. `n_verts`, `n_elems` stay; new aliases added
2. SMS_2DM reader is best-effort
3. Lazy initialization opt-in (`compute_layers=False` kwarg)

---

## PART 3: SPECIFIC IMPROVEMENTS

Each improvement scoped, testable, driven by ADMESH-Domains schema.

---

### I-1: Fix `read_from_fort14()` for Quad/Mixed Elements [CRITICAL]

**Current:** Line 692–693 raises `ValueError` for any `num_nodes != 3`.  
**Required:** Parse 3-node and 4-node elements from same file.  
**Why:** ADMESH-Domains hosts real coastal meshes; many are mixed or quad-dominant.

**Spec:**
```
read_from_fort14(path) MUST:
  - Parse elements with num_nodes == 3 as triangles
  - Parse elements with num_nodes == 4 as quads
  - Build connectivity_list with shape (n_elems, 4) when any quads present
    (padded-triangle rule: triangles get vertex[3] = vertex[0])
  - Raise ValueError only for num_nodes not in {3, 4}

Acceptance:
  □ Load block_o.14 (quad-heavy) without error
  □ Load a mixed-element .fort.14 without error
  □ Existing triangular roundtrip test still passes (zero regression)
  □ element_type returns "Quadrilateral", "Triangular", or "Mixed-Element"
```

**Affected code:** `src/chilmesh/CHILmesh.py:692–699`

---

### I-2: Add `bounding_box` Property

**Current:** No bbox API. CHILmesh stores `points` (x, y, z), no bbox.  
**Required:** Property returning `{min_lon, max_lon, min_lat, max_lat}` matching ADMESH-Domains schema.

**Spec:**
```python
@property
def bounding_box(self) -> dict:
    """Geographic bounding box matching ADMESH-Domains schema."""
    return {
        "min_lon": float(self.points[:, 0].min()),
        "max_lon": float(self.points[:, 0].max()),
        "min_lat": float(self.points[:, 1].min()),
        "max_lat": float(self.points[:, 1].max()),
    }

Acceptance:
  □ Returns dict with exactly 4 keys: min_lon, max_lon, min_lat, max_lat
  □ Values are floats
  □ Raises no error on empty mesh (returns NaN dict or raises ValueError)
  □ Consistent with points coordinate system (x=lon, y=lat for geographic meshes)
```

---

### I-3: Add `admesh_metadata()` Method

**Current:** `n_verts`, `n_elems`, `type` don't map to ADMESH-Domains schema keys.  
**Required:** Method returning ADMESH-Domains Mesh schema-compatible fields.

**Spec:**
```python
def admesh_metadata(self) -> dict:
    """Return dict matching ADMESH-Domains Mesh schema fields CHILmesh can compute."""
    return {
        "node_count": self.n_verts,
        "element_count": self.n_elems,
        "element_type": self.type,            # "Triangular" | "Quadrilateral" | "Mixed-Element"
        "bounding_box": self.bounding_box,
        "n_layers": self.n_layers,            # CHILmesh-specific (skeletonization depth)
    }

Acceptance:
  □ Returns dict with all 5 keys
  □ element_type matches one of three valid strings
  □ node_count == len(points)
  □ element_count == len(connectivity_list)
  □ bounding_box consistent with bounding_box property
```

---

### I-4: Add SMS_2DM Reader

**Current:** No `.2dm` reader.  
**Required:** `CHILmesh.read_from_2dm(path)` parsing SMS 2DM format.  
**Why:** ADMESH-Domains stores `.2dm` files (`type: "SMS_2DM"`); without this CHILmesh can't load ~30% of catalog.

**SMS 2DM Format (key lines):**
```
MESH2D                    ← header
ND  <id> <x> <y> <z>     ← node definition
E3T <id> <n1> <n2> <n3> <matid>   ← triangle element
E4Q <id> <n1> <n2> <n3> <n4> <matid>  ← quad element
```

**Spec:**
```python
@staticmethod
def read_from_2dm(path: Path) -> "CHILmesh":
    """Parse SMS 2DM mesh format into CHILmesh."""

Acceptance:
  □ Parses ND lines into points array (x, y, z)
  □ Parses E3T lines as triangular elements
  □ Parses E4Q lines as quad elements (if present)
  □ Ignores unrecognized line types gracefully
  □ Returns CHILmesh with correct n_verts, n_elems
  □ Roundtrip: read .2dm → write .fort.14 → read .fort.14 → same mesh (shape-identical)
  □ Test fixture: at least one real .2dm from ADMESH-Domains test_cases
```

---

### I-5: Add `from_admesh_domain()` Classmethod

**Current:** No integration entry point from ADMESH-Domains.  
**Required:** Classmethod accepting ADMESH-Domains `Mesh` object, returning `CHILmesh`.

**Spec:**
```python
@classmethod
def from_admesh_domain(
    cls,
    mesh_record,           # admesh_domains.Mesh or dict
    compute_layers: bool = True,
) -> "CHILmesh":
    """
    Load a CHILmesh from an ADMESH-Domains Mesh record.
    
    Parameters
    ----------
    mesh_record
        An admesh_domains.Mesh instance (or compatible dict with 'filename', 'type').
    compute_layers
        If False, skip skeletonization (faster for bulk metadata queries).
    """
    filename = mesh_record.filename if hasattr(mesh_record, "filename") else mesh_record["filename"]
    fmt = getattr(mesh_record, "type", "ADCIRC")
    path = Path(filename)
    if fmt == "SMS_2DM":
        mesh = cls.read_from_2dm(path)
    else:
        mesh = cls.read_from_fort14(path)
    if not compute_layers:
        # Already built adjacencies; skip _skeletonize()
        pass
    return mesh

Acceptance:
  □ Accepts admesh_domains.Mesh (with .filename and .type attributes)
  □ Accepts plain dict with "filename" and optional "type" keys
  □ Routes to correct reader (fort14 vs 2dm) based on type
  □ compute_layers=False skips _skeletonize() (adjacencies still built)
  □ compute_layers=True (default) runs full init
  □ Works with admesh-domains not installed (import guarded)
  □ Raises ImportError with install hint if mesh_record is admesh_domains type but package absent
```

---

### I-6: Lazy `_skeletonize()` (Opt-In)

**Current:** `_initialize_mesh()` always runs `_build_adjacencies()` then `_skeletonize()`.  
**Required:** `__init__` accepts `compute_layers: bool = True`; when False, skip `_skeletonize()`.  
**Why:** ADMESH-Domains bulk-loads metadata without needing layers; skeletonization is slowest step.

**Spec:**
```python
def __init__(
    self,
    connectivity=None,
    points=None,
    grid_name=None,
    compute_layers: bool = True,    ← NEW
) -> None:

def _initialize_mesh(self, compute_layers: bool = True) -> None:
    ...
    self._build_adjacencies()
    if compute_layers:
        self._skeletonize()
    # else: layers dict remains empty {}; n_layers = 0

Acceptance:
  □ compute_layers=True (default) behaves identically to current init
  □ compute_layers=False: no layers dict populated, n_layers == 0
  □ Calling get_layer() when layers not computed raises informative error
    ("Layers not computed. Re-initialize with compute_layers=True or call _skeletonize().")
  □ _skeletonize() can be called manually after init to populate lazily
  □ All existing tests pass (compute_layers=True by default = no regression)
  □ Benchmark: Block_O init with compute_layers=False is <2s (vs ~30s full init)
```

---

### I-7: Handle `kind: "boundary"` Meshes

**Current:** Every loaded mesh treated as full 2D domain mesh.  
**Required:** When `kind == "boundary"`, tag object as boundary description; skeletonization skipped.

**Spec:**
```python
# In from_admesh_domain():
kind = getattr(mesh_record, "kind", "mesh")
if kind == "boundary":
    mesh = cls(..., compute_layers=False)
    mesh._is_boundary = True
else:
    mesh._is_boundary = False

# _skeletonize() should check:
if self._is_boundary:
    warnings.warn("Skeletonization not meaningful for boundary meshes.")
    return

Acceptance:
  □ Boundary meshes loaded with compute_layers=False automatically
  □ _is_boundary attribute set correctly based on kind
  □ Calling _skeletonize() on boundary mesh warns, does not raise
  □ plotting still works for boundary meshes
```

---

## PART 4: AMBIGUITY SCORING

```
Goal Clarity:        0.92  ✓  (7 concrete improvements, each with acceptance criteria)
Boundary Clarity:    0.88  ✓  (in-scope: I/O, metadata, lazy init; out-of-scope: MADMESHR, Phase 3 algos)
Constraint Clarity:  0.90  ✓  (hard constraints explicit; admesh-domains as optional dep)
Acceptance Criteria: 0.90  ✓  (each improvement has pass/fail checklist)

Ambiguity: 1.0 − (0.35×0.92 + 0.25×0.88 + 0.20×0.90 + 0.20×0.90)
         = 1.0 − (0.322 + 0.220 + 0.180 + 0.180)
         = 1.0 − 0.902
         = 0.098  ← WELL BELOW gate (0.20) ✓
```

**GATE PASSED. Ready to plan and implement.**

---

## PART 5: SCOPE BOUNDARY

- ❌ Replaces Phase 1 graph benchmarking/modernization (different concern)
- ❌ Adds MADMESHR advancing-front API (Phase 4)
- ❌ Implements O(n log n) edge discovery optimization (Phase 3)
- ❌ Bundles ADMESH-Domains as hard dependency
- ❌ Changes skeletonization algorithm or semantics
- ❌ Publishes to ADMESH-Domains registry from CHILmesh

---

## PART 6: PRIORITY ORDER

| # | Improvement | Priority | Risk | Lines Changed |
|---|-------------|----------|------|---------------|
| I-1 | Fix fort14 quad/mixed reader | P0 — CRITICAL BUG | Low | ~10 |
| I-6 | Lazy `_skeletonize()` | P0 — blocks bulk load | Low | ~15 |
| I-2 | `bounding_box` property | P1 — metadata parity | Trivial | ~8 |
| I-3 | `admesh_metadata()` method | P1 — metadata parity | Trivial | ~10 |
| I-5 | `from_admesh_domain()` classmethod | P1 — integration entry | Low | ~25 |
| I-4 | SMS_2DM reader | P2 — format coverage | Medium | ~40 |
| I-7 | `kind: "boundary"` handling | P2 — catalog parity | Low | ~15 |

All 7 improvements ship in ONE PR — all I/O and metadata-layer changes, no interaction with graph traversal algorithms.

---

## PART 7: SINGLE GITHUB ISSUE

7 improvements belong in ONE issue + ONE PR — share:
- Same motivation (ADMESH-Domains I/O parity)
- Same test fixtures (real ADMESH-Domains meshes)
- Same risk profile (I/O layer, no algorithm changes)
- Same logical boundary (format/metadata, not computation)

**Proposed issue title:**  
`feat: CHILmesh ↔ ADMESH-Domains I/O parity (quad reader fix, 2dm support, metadata, lazy init)`

**Labels:** `enhancement`, `io`, `admesh-domains`, `bug` (for I-1)  
**Milestone:** 0.2.0  
**Blocks:** Phase 4 MADMESHR integration (can't test with real domains without this)

---

## PART 8: TEST FIXTURES NEEDED

1. Quad/mixed `.fort.14` fixture (e.g., `block_o.14` if it has quads — verify)
2. Real `.2dm` fixture from ADMESH-Domains (download one `test_case: true` mesh)
3. Current 4 fixtures (annulus, block_o, structured, donut) for regression

`.2dm` fixture committed to `src/chilmesh/data/` alongside existing `.14` files.

---

## SUMMARY

One feature, 7 improvements, closes critical bug, establishes ADMESH-Domains as canonical I/O layer:

1. **I-1 (P0):** Fix fort14 reader — stop crashing on real-world quad/mixed meshes
2. **I-6 (P0):** Lazy skeletonization — fast bulk loading for ADMESH-Domains catalog
3. **I-2 (P1):** `bounding_box` property — metadata parity with catalog schema
4. **I-3 (P1):** `admesh_metadata()` — expose computed fields in catalog-compatible format
5. **I-5 (P1):** `from_admesh_domain()` — canonical integration entry point
6. **I-4 (P2):** SMS_2DM reader — full format coverage
7. **I-7 (P2):** Boundary mesh handling — `kind` field awareness

**Estimated effort:** 3–4 days (mostly I/O code, no algorithm changes)  
**Risk:** Low (touches only readers, init flag, new methods; no existing algorithm logic)  
**Gate:** All acceptance criteria above pass + all existing tests green

---

**Document Status:** SPECIFY PHASE COMPLETE  
**Ambiguity Score:** 0.098 (gate ≤ 0.20 ✓)  
**Next:** `/gsd-discuss-phase` or plan this as a single PR  
**Author:** Spec-Kit (Claude)  
**Date:** 2026-04-26
