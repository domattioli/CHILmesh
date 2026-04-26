# ADMESH-Domains: Suggested Issues (Cross-Repo)

**From:** CHILmesh I/O spec analysis  
**Create these at:** https://github.com/domattioli/ADMESH-Domains/issues/new  
**Date:** 2026-04-26

---

## Issue 1 — FEATURE: `Mesh.to_chilmesh()` integration method

**Title:** `feat: Mesh.to_chilmesh() — canonical CHILmesh integration entry point`  
**Labels:** `enhancement`

**Body:**
ADMESH-Domains is the catalog/delivery layer; CHILmesh is the computational engine. Right now there is no bridge. `Mesh.load()` returns raw bytes/path and users must manually call `CHILmesh.read_from_fort14()` or `read_from_2dm()`.

**Proposed:**
```python
# Option A: separate method
mesh_obj = domain.meshes[0].to_chilmesh()

# Option B: kwarg on existing load()
mesh_obj = domain.meshes[0].load(as_chilmesh=True)
```

**Implementation notes:**
- `chilmesh` as soft optional dep: `pip install admesh-domains[chilmesh]`
- Route ADCIRC/ADCIRC_GRD → `CHILmesh.read_from_fort14()`
- Route SMS_2DM → `CHILmesh.read_from_2dm()`
- Pass `compute_layers=False` for `kind="boundary"` meshes automatically
- CHILmesh-side entry point: `CHILmesh.from_admesh_domain()` (CHILmesh issue #43)

**Acceptance:**
- [ ] Returns a live `CHILmesh` instance
- [ ] Works for ADCIRC and SMS_2DM types
- [ ] `ImportError` with install hint if `chilmesh` not installed
- [ ] Documented in README and HF dataset card

---

## Issue 2 — FEATURE: Auto-compute metadata from mesh file on contribution

**Title:** `feat: CI auto-computes node_count, element_count, bounding_box, element_type from mesh file`  
**Labels:** `enhancement`, `ci`

**Body:**
`node_count`, `element_count`, and `bounding_box` are optional schema fields filled in manually. This creates:
1. **Stale data** — contributor enters wrong count; file changes; schema drifts
2. **Missing data** — many meshes have `null`, breaking queries like "find meshes with < 100K nodes"

**Proposed:** On PR submission, run CHILmesh on each new/modified mesh file to auto-compute:
```python
m = CHILmesh.from_admesh_domain(mesh_record, compute_layers=False)
computed = m.admesh_metadata()
# → {node_count, element_count, element_type, bounding_box}
```
Compare against schema values and warn/block on mismatch.

**Acceptance:**
- [ ] CI detects mismatches between schema `node_count` and actual file
- [ ] Contribution workflow auto-fills the 4 computed fields
- [ ] Backfill migration for existing entries
- [ ] One-time migration script provided

---

## Issue 3 — BUG: `element_type` is unconstrained free-form string

**Title:** `bug: element_type schema field needs standardised vocabulary (enum constraint)`  
**Labels:** `bug`, `schema`

**Body:**
`element_type` is `Optional[str]` with no enum. Contributors write anything: "triangles", "tri", "TRI", "Triangular", "mixed", etc. This breaks programmatic queries and CHILmesh integration.

**CHILmesh uses a specific vocabulary:**
- `"Triangular"` — all triangular elements
- `"Quadrilateral"` — all quad elements
- `"Mixed-Element"` — both in same mesh

**Proposed fix:**
```python
from typing import Literal, Optional
element_type: Optional[Literal["Triangular", "Quadrilateral", "Mixed-Element"]] = None
```

**Acceptance:**
- [ ] Enum constraint added to schema
- [ ] Validation rejects any other string value
- [ ] Existing entries migrated to enum vocabulary
- [ ] Documented in schema reference

---

## Issue 4 — FEATURE: Add `n_layers` (skeletonization depth) to Mesh schema

**Title:** `feat: add chilmesh_n_layers as optional Mesh schema field`  
**Labels:** `enhancement`, `schema`

**Body:**
CHILmesh computes `n_layers` — the number of concentric boundary-peeling layers (skeletonization depth). This is a meaningful topological complexity measure:
- Thin estuary: `n_layers = 3`
- Wide open-ocean domain: `n_layers = 40`
- Deeply nested channel at pinch point: `n_layers = 2`

It enables MADMESHR training curriculum (filter by complexity tier) and quality screening (very low n_layers may indicate degenerate geometry).

**Proposed:**
```python
chilmesh_n_layers: Optional[int] = None
```

Prefix `chilmesh_` signals that it's computed by CHILmesh, not a universal mesh format field.

**Acceptance:**
- [ ] Field added as `Optional[int]` with note on provenance
- [ ] Contribution CI populates it via CHILmesh when installed
- [ ] `find_meshes()` supports filtering by `n_layers` range
- [ ] `kind="boundary"` always yields `n_layers = 0`

---

## Issue 5 — BUG: `type: "ADCIRC_GRD"` enum value is undocumented

**Title:** `bug: ADCIRC_GRD type enum is undocumented — clarify format and CHILmesh routing`  
**Labels:** `bug`, `documentation`

**Body:**
The Mesh `type` field accepts `"ADCIRC"`, `"SMS_2DM"`, and `"ADCIRC_GRD"`. The third is not documented anywhere.

Questions that block CHILmesh integration:
1. What file extension(s) does `"ADCIRC_GRD"` correspond to?
2. Is it structurally identical to `"ADCIRC"` (`.fort.14`) or a different layout?
3. Should `CHILmesh.from_admesh_domain()` route it to `read_from_fort14()`?

Without this, CHILmesh must either refuse to load `ADCIRC_GRD` meshes or make an undocumented assumption.

**Acceptance:**
- [ ] All three `type` values documented with file extensions and format notes
- [ ] If `ADCIRC_GRD` ≡ `ADCIRC`, mark as alias and document
- [ ] CHILmesh routing clarified (CHILmesh issue #43 updated accordingly)

---

## Issue 6 — FEATURE: Integration test — load all `test_case=true` meshes into CHILmesh

**Title:** `feat: integration test suite — validate test_case meshes load correctly into CHILmesh`  
**Labels:** `enhancement`, `testing`

**Body:**
`test_case: true` meshes are the canonical fixtures. There is no test that actually loads them into CHILmesh and validates the result. An integration suite would catch:
- Schema metadata mismatches (wrong `node_count`, `element_type`)
- CHILmesh reader failures on real-world files
- Bounding box drift after file updates

**Proposed test:**
```python
# tests/test_chilmesh_integration.py
@pytest.mark.integration
@pytest.mark.parametrize("mesh_record", find_meshes(test_case=True))
def test_load_into_chilmesh(mesh_record):
    m = chilmesh.CHILmesh.from_admesh_domain(mesh_record)
    assert m.n_verts == mesh_record.node_count
    assert m.n_elems == mesh_record.element_count
    assert m.type    == mesh_record.element_type
    bb = m.bounding_box
    assert bb["min_lon"] <= bb["max_lon"]
```

Marked `@pytest.mark.integration` — runs nightly or on new-mesh PRs, not in fast CI.

**Depends on:** CHILmesh #40 (quad reader), #42 (admesh_metadata), #43 (from_admesh_domain)

**Acceptance:**
- [ ] Integration test file added to `tests/`
- [ ] All `test_case=True` meshes pass
- [ ] Schema metadata validated against CHILmesh computed values
- [ ] CI job configured

---

*These issues are the ADMESH-Domains side of the CHILmesh↔ADMESH-Domains integration story.*  
*CHILmesh issues #40–44 are the CHILmesh-side counterparts.*
