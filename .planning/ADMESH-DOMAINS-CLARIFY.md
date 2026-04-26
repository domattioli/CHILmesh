# Spec-Kit: Clarify Phase
# "CHILmesh ↔ ADMESH-Domains I/O Parity"

**Branch:** `planning-optimize_modernize`
**Date:** 2026-04-26
**Perspective:** Boundary Keeper → Failure Analyst → Seed Closer
**Prerequisite:** ADMESH-DOMAINS-IO-SPEC.md (Specify phase, ambiguity 0.098)

---

## CLARIFY ROUND 1 — Boundary Keeper

*"What did the scouting reveal that changes the spec boundary?"*

### Finding C-1: `write_fort14()` ALSO crashes on quads (second critical bug)

**Evidence:** `src/chilmesh/CHILmesh.py:990`
```python
for i, tri in enumerate(elements, start=1):
    n1, n2, n3 = tri + 1   # ← crashes if tri has 4 values
    f.write(f"{i} 3 {n1} {n2} {n3}\n")
```

**Impact:** Even after fixing the reader (I-1), loading a quad mesh and
then writing it back will raise `ValueError: too many values to unpack`.
The `write_to_fort14()` docstring explicitly says `elements: (n_elems, 3)`.

**Clarification:** I-1 (fix reader) **must be paired with fixing the writer**.
These are one atomic fix, not two separate issues.

**Scope change:** Issue #40 expanded to cover both reader AND writer.

---

### Finding C-2: No quad/mixed test fixture exists

**Evidence:** All 4 bundled fixtures have `3` in every element line:
```
Block_O     5214 elems, all "i 3 n1 n2 n3"
annulus      580 elems, all "i 3 n1 n2 n3"
donut        276 elems, all "i 3 n1 n2 n3"
structured   660 elems, all "i 3 n1 n2 n3"
```

**Impact:** There is no way to regression-test quad/mixed-element code today.
Issue #40 cannot be verified without a quad fixture.

**Clarification required (OPEN QUESTION Q-1):**
> Should we create a synthetic quad `.fort.14` fixture committed to
> `src/chilmesh/data/`, or acquire a real one from ADMESH-Domains?

**Recommendation:** Create a minimal synthetic 2×2 quad grid as a test
fixture (4 nodes, 1 quad element). Simple, deterministic, no download needed.
A real-world quad from ADMESH-Domains can be added later as an integration test.

---

### Finding C-3: Coordinates are Cartesian, NOT geographic

**Evidence:** Fixture coordinates:
```
structured: x in [-2.0, 1.0],  y in [-1.0, 12.0]   ← clearly Cartesian
Block_O:    x in [5.0, 10.0],  y in [9.9, 44.0]     ← Cartesian (km scale)
```

**Impact:** The Specify phase spec for `bounding_box` assumed `x=lon, y=lat`.
That is wrong for the bundled fixtures and likely wrong for many ADMESH-Domains
meshes that use projected coordinates (UTM, state plane, etc.).

**Clarification:** `bounding_box` should return `{min_x, max_x, min_y, max_y}`,
not `{min_lon, max_lon, min_lat, max_lat}`. Coordinate interpretation
(geographic vs. Cartesian) is the caller's responsibility.

**ADMESH-Domains schema note:** The schema's `bounding_box` uses lat/lon field
names — this is a mismatch to raise as a cross-repo question (see Q-2 below).

---

### Finding C-4: ADCIRC_GRD is likely identical to ADCIRC format

**Evidence:** All bundled ADCIRC-style files use identical line format
regardless of extension (`.14` vs `.fort.14`). The ADMESH-Domains schema
lists `ADCIRC_GRD` as a distinct type but provides no documentation.
The `.grd` extension in the ADCIRC ecosystem is used as an informal alias
for `.14` (the "grid" file).

**Clarification:** Treat `ADCIRC_GRD` as routed to `read_from_fort14()` in
`from_admesh_domain()`, same as `ADCIRC`. Add a code comment and flag in the
ADMESH-Domains issue (already written to `ADMESH-DOMAINS-ISSUES.md`).

**No spec change needed** — routing logic in I-5 already defaults to
`read_from_fort14()` for any non-SMS_2DM type.

---

## CLARIFY ROUND 2 — Failure Analyst

*"What is the worst thing that breaks if the spec is wrong?"*

### Failure Mode F-1: quad write crashes silently on old code paths

**Scenario:** User loads a quad `.fort.14` from ADMESH-Domains (after I-1
fixes the reader), runs a computation, then calls `write_to_fort14()`.
The old `write_fort14()` (line 990) crashes with:
```
ValueError: too many values to unpack (expected 3)
```

**How bad:** Crash with data loss (computed mesh not saved). Silent regression
if users have patched workarounds.

**Prevention:** Fix writer in same PR as reader (atomic). Add regression test:
`test_fort14_quad_roundtrip` — load quad fixture, write, reload, compare.

---

### Failure Mode F-2: `bounding_box` breaks for meshes with z ≠ 0

**Scenario:** Some coastal meshes (e.g., bathymetric) have real z-values
(depth). `bounding_box` only reads `points[:, 0]` and `points[:, 1]`.
This is correct — z is not part of the horizontal bounding box.

**How bad:** No failure — z is correctly excluded. Document this explicitly
in the docstring to prevent future "why isn't z included?" questions.

---

### Failure Mode F-3: `from_admesh_domain()` called before file is downloaded

**Scenario:** ADMESH-Domains uses HuggingFace on-demand download. If a user
calls `from_admesh_domain(mesh_record)` before calling `mesh_record.load()`,
the file may not exist on disk.

**How bad:** `FileNotFoundError` — acceptable, but error message should guide
the user: *"File not found. Call mesh_record.load() first to download it."*

**Prevention:** Wrap `Path(filename).exists()` check in `from_admesh_domain()`
and raise with explicit guidance.

---

### Failure Mode F-4: `compute_layers=False` then user calls `get_layer()`

**Scenario:** User loads with `compute_layers=False` for speed, forgets,
then calls `mesh.get_layer(0)` — currently returns `{}` or raises `KeyError`.

**Current code:** `get_layer()` checks `layer_idx < 0 or layer_idx >= self.n_layers`.
With `n_layers = 0`, any `layer_idx` triggers `ValueError("out of range [0, -1]")`.
The error message `-1` is confusing.

**Prevention:** Add explicit guard in `get_layer()`:
```python
if self.n_layers == 0:
    raise RuntimeError(
        "Layers not computed. Call _skeletonize() or "
        "re-initialize with compute_layers=True."
    )
```

---

## CLARIFY ROUND 3 — Seed Closer

*"What unresolved territory must be decided before planning?"*

### Open Question Q-1: Quad test fixture — synthetic or real-world?

**Options:**
- A. Synthetic 2×2 quad grid (4 nodes, 1 element) — minimal, committed to repo
- B. Smallest quad mesh from ADMESH-Domains (requires download in CI)
- C. Both — synthetic for fast unit tests, ADMESH-Domains mesh for integration

**Recommendation:** Option C. Synthetic fixture for #40 unit test. Integration
test (separate issue) uses real ADMESH-Domains mesh.

**Decision needed before:** Writing tests for #40.

---

### Open Question Q-2: `bounding_box` key names — x/y vs lon/lat?

**Context:** CHILmesh fixtures use Cartesian coordinates.
ADMESH-Domains schema uses `{min_lat, max_lat, min_lon, max_lon}`.

**Options:**
- A. CHILmesh uses `{min_x, max_x, min_y, max_y}` — correct for Cartesian
- B. CHILmesh uses `{min_lon, max_lon, min_lat, max_lat}` — matches ADMESH-Domains schema
- C. CHILmesh uses both via a `geographic: bool` kwarg

**Recommendation:** Option A — CHILmesh stays coordinate-agnostic (`x`/`y`).
ADMESH-Domains schema can map `x→lon`, `y→lat` in `Mesh.to_chilmesh()` or
via a thin adapter. Document the convention mismatch in `admesh_metadata()`.

**Decision needed before:** Implementing I-2 (bounding_box property).

---

### Open Question Q-3: Does `admesh-domains` become an optional dep?

**Context:** `from_admesh_domain()` accepts a `Mesh` record duck-typed
via `getattr`. It does NOT require importing `admesh_domains`.

**Options:**
- A. No import ever — duck-typing only, `admesh-domains` stays external
- B. Optional import — `try: import admesh_domains` for type hints only
- C. Hard dependency — `pip install chilmesh` pulls in `admesh-domains`

**Recommendation:** Option A. No import. Pure duck-typing on `getattr(record, "filename", None)`.
The `from_admesh_domain()` docstring notes the expected attributes. No `extras_require` needed.

**Decision needed before:** Implementing I-5.

---

## UPDATED AMBIGUITY SCORES (Post-Clarify)

```
Goal Clarity:        0.95  ✓  (writer bug found and scoped, no new unknowns)
Boundary Clarity:    0.92  ✓  (quad fixture decision deferred but bounded)
Constraint Clarity:  0.93  ✓  (coordinate convention resolved: x/y not lon/lat)
Acceptance Criteria: 0.93  ✓  (F-4 guard added, roundtrip test added)

Ambiguity: 1.0 − (0.35×0.95 + 0.25×0.92 + 0.20×0.93 + 0.20×0.93)
         = 1.0 − (0.3325 + 0.23 + 0.186 + 0.186)
         = 1.0 − 0.9345
         = 0.0655  ← gate ≤ 0.20 ✓
```

**Gate passed. Spec + Clarify complete.**
