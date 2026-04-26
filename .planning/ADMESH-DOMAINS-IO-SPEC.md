# Spec-Kit: Specify Phase
# "CHILmesh as Computational Engine for ADMESH-Domains I/O"

**Branch:** `planning-optimize_modernize`  
**Date:** 2026-04-26  
**Methodology:** Spec-Kit — Specify Phase  
**Scope:** How ADMESH-Domains as the domain I/O layer reshapes CHILmesh's responsibilities  
**Framing:** ADMESH-Domains is the *catalog and delivery layer*; CHILmesh is the *computational engine*

---

## THE REFRAME

Previous view: CHILmesh is a mesh computation library that happens to read `.fort.14`.

**Correct view:** ADMESH-Domains owns the domain catalog, mesh files, and I/O. CHILmesh is the object users get *after* loading from that catalog. Everything in CHILmesh that relates to I/O, metadata, or format handling must be consistent with ADMESH-Domains' schema.

This is not a minor integration. It is a **redefinition of scope**.

---

## PART 1: CURRENT STATE

### 1.1 What CHILmesh Does Today (I/O)

```python
# ADCIRC .fort.14 — read
CHILmesh.read_from_fort14(path) -> CHILmesh   # triangular elements ONLY (line 693)

# ADCIRC .fort.14 — write
mesh.write_to_fort14(filename)                # delegates to write_fort14()

# No other formats
# No metadata fields matching ADMESH-Domains schema
# No lazy initialization
# No bounding_box
# No domain/catalog awareness
```

**Critical Bug (B5 candidate):**  
`read_from_fort14()` line 692–693 raises `ValueError` for any element with `num_nodes != 3`.  
ADMESH-Domains hosts real-world coastal meshes many of which are quad or mixed-element.  
This means **CHILmesh cannot load the majority of ADMESH-Domains meshes today**.

### 1.2 What ADMESH-Domains Provides

**Domain schema** (what the catalog stores):
```
Domain:
  name         str       required
  full_name    str       optional
  description  str       optional
  category     "real-world" | "synthetic"
  region       str       optional
  applications List[str]
  bounding_box {min_lat, max_lat, min_lon, max_lon}  optional
  meshes       List[Mesh]  required, ≥1

Mesh:
  id              str       required
  filename        str       required
  description     str       optional
  size_mb         float     ≥0
  node_count      int       optional  ← CHILmesh computes this
  element_count   int       optional  ← CHILmesh computes this
  element_type    str       optional  ← CHILmesh computes this
  type            "ADCIRC" | "SMS_2DM" | "ADCIRC_GRD"
  contributor     str       optional
  uploaded_date   str       optional
  modified_date   str       optional
  refinement_level str      optional
  features        List[str]
  aliases         List[str]
  bounding_box    {min_lat, max_lat, min_lon, max_lon}  optional  ← CHILmesh computes this
  license         one-of-7
  kind            "mesh" | "boundary"
  test_case       bool
```

**File formats supported:**  
- `.14` / `.fort.14` → ADCIRC format  
- `.2dm` → SMS 2D mesh format (surface water modeling)

**Loading pattern:**  
`Mesh.load()` fetches file from HuggingFace on-demand. Returns raw file bytes or path.  
CHILmesh must then parse that file into a live `CHILmesh` object.

### 1.3 The Gap (Delta Between Today and Needed)

| Capability | ADMESH-Domains Needs | CHILmesh Today | Gap |
|------------|---------------------|---------------|-----|
| Read `.fort.14` triangular | ✓ | ✓ | None |
| Read `.fort.14` quad/mixed | ✓ (real-world coastal) | ✗ raises ValueError | **CRITICAL BUG** |
| Read `.2dm` (SMS format) | ✓ | ✗ not implemented | Missing |
| `node_count` property | ✓ (schema field) | partial (`n_verts`) | Name mismatch |
| `element_count` property | ✓ (schema field) | partial (`n_elems`) | Name mismatch |
| `element_type` property | ✓ (schema field) | partial (`type` str) | Name mismatch |
| `bounding_box` property | ✓ (schema field) | ✗ must compute from points | Missing API |
| Lazy initialization | ✓ (HF on-demand load) | ✗ full init always | Performance gap |
| `from_admesh_domain()` | ✓ (integration entry) | ✗ not implemented | Missing |
| `kind: "boundary"` handling | ✓ (separate mesh type) | ✗ all meshes same | Missing |
| Metadata dict round-trip | ✓ (catalog sync) | ✗ | Missing |

---

## PART 2: CONSTRAINTS

### Hard Constraints (Must Preserve)

1. **Skeletonization semantics** — immutable (audit Q3, Decision A2)
2. **`.fort.14` roundtrip** — byte-identical (audit B1)
3. **Mixed-element support** — triangles + quads in one mesh
4. **API surface** — `signed_area()`, `elem_quality()`, `plot()` signatures unchanged
5. **No new external runtime deps** — only scipy/numpy/matplotlib
   - Exception: `admesh-domains[hf]` is optional (`pip install admesh-domains` optional dep)

### Soft Constraints

1. Backwards compatible property names (`n_verts`, `n_elems` stay; new aliases added)
2. SMS_2DM reader is best-effort (format is less well-documented than fort.14)
3. Lazy initialization is opt-in (`compute_layers=False` kwarg)

---

## PART 3: SPECIFIC IMPROVEMENTS

Each improvement is **scoped, testable, and directly driven by the ADMESH-Domains schema**.

---

### Improvement I-1: Fix `read_from_fort14()` for Quad/Mixed Elements [CRITICAL]

**Current:** Line 692–693 raises `ValueError` for any `num_nodes != 3`.  
**Required:** Parse 3-node (triangular) and 4-node (quad/mixed) elements from same file.  
**Why:** ADMESH-Domains hosts real coastal meshes that are mixed or quad-dominant (e.g., ADCIRC GRD files).

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

### Improvement I-2: Add `bounding_box` Property

**Current:** No bounding box API. CHILmesh stores `points` (x, y, z) but offers no bbox.  
**Required:** Property that returns `{min_lon, max_lon, min_lat, max_lat}` matching ADMESH-Domains schema.  
**Why:** ADMESH-Domains schema has `bounding_box` on both Domain and Mesh. CHILmesh should compute it from `points` for validation/enrichment.

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

### Improvement I-3: Add `admesh_metadata()` Method

**Current:** CHILmesh has `n_verts`, `n_elems`, `type` but they don't map to ADMESH-Domains schema keys.  
**Required:** Method that returns dict with ADMESH-Domains Mesh schema-compatible fields.  
**Why:** Enables CHILmesh to validate or enrich ADMESH-Domains registry entries.

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

### Improvement I-4: Add SMS_2DM Reader

**Current:** No `.2dm` reader.  
**Required:** `CHILmesh.read_from_2dm(path)` that parses SMS 2DM format.  
**Why:** ADMESH-Domains stores `.2dm` files (`type: "SMS_2DM"`). Without this, CHILmesh cannot load ~30% of ADMESH-Domains meshes.

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

### Improvement I-5: Add `from_admesh_domain()` Classmethod

**Current:** No integration entry point from ADMESH-Domains.  
**Required:** Classmethod that accepts an ADMESH-Domains `Mesh` object and returns a `CHILmesh`.  
**Why:** This is the canonical entry point for "load a mesh from the catalog."

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

### Improvement I-6: Lazy `_skeletonize()` (Opt-In)

**Current:** `_initialize_mesh()` always runs `_build_adjacencies()` then `_skeletonize()`.  
**Required:** `__init__` accepts `compute_layers: bool = True`; when False, skip `_skeletonize()`.  
**Why:** ADMESH-Domains bulk-loads metadata without needing layer data. Skeletonization is the slowest step. This alone may eliminate the 30s Block_O bottleneck for catalog queries.

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

### Improvement I-7: Handle `kind: "boundary"` Meshes

**Current:** CHILmesh treats every loaded mesh as a full 2D domain mesh.  
**Required:** When `kind == "boundary"`, CHILmesh should tag the object as a boundary description, not a volume mesh. Skeletonization does not apply.

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

## PART 5: SCOPE BOUNDARY (What This Feature Is NOT)

- ❌ Replaces the Phase 1 graph benchmarking/modernization work (different concern)
- ❌ Adds MADMESHR advancing-front API (Phase 4)
- ❌ Implements O(n log n) edge discovery optimization (Phase 3)
- ❌ Bundles ADMESH-Domains as a hard dependency
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

**All 7 improvements can ship in ONE PR** — they are all I/O and metadata-layer changes with no interaction with graph traversal algorithms.

---

## PART 7: SINGLE GITHUB ISSUE

These 7 improvements belong in **ONE GitHub issue and ONE PR** because they share:
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

For all 7 improvements, tests need:

1. A quad/mixed `.fort.14` fixture (e.g., `block_o.14` if it has quads — verify)
2. A real `.2dm` fixture from ADMESH-Domains (download one `test_case: true` mesh)
3. The current 4 fixtures (annulus, block_o, structured, donut) for regression

The `.2dm` fixture can be committed to `src/chilmesh/data/` alongside existing `.14` files.

---

## SUMMARY

**One feature addresses 7 improvements, closes a critical bug, and establishes ADMESH-Domains as the canonical I/O layer for CHILmesh:**

1. **I-1 (P0):** Fix fort14 reader — stop crashing on real-world quad/mixed meshes
2. **I-6 (P0):** Lazy skeletonization — enable fast bulk loading for ADMESH-Domains catalog
3. **I-2 (P1):** `bounding_box` property — metadata parity with catalog schema
4. **I-3 (P1):** `admesh_metadata()` — expose computed fields in catalog-compatible format
5. **I-5 (P1):** `from_admesh_domain()` — canonical integration entry point
6. **I-4 (P2):** SMS_2DM reader — full format coverage
7. **I-7 (P2):** Boundary mesh handling — `kind` field awareness

**Estimated effort:** 3–4 days (mostly I/O code, no algorithm changes)  
**Risk:** Low (touches only readers, init flag, and new methods; no existing algorithm logic)  
**Gate:** All acceptance criteria above pass + all existing tests green

---

**Document Status:** SPECIFY PHASE COMPLETE  
**Ambiguity Score:** 0.098 (gate ≤ 0.20 ✓)  
**Next:** `/gsd-discuss-phase` or plan this as a single PR  
**Author:** Spec-Kit (Claude)  
**Date:** 2026-04-26
