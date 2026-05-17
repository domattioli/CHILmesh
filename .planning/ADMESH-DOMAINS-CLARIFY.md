# Spec-Kit: Clarify Phase
# "CHILmesh ↔ ADMESH-Domains I/O Parity"

**Branch:** `planning-optimize_modernize`
**Date:** 2026-04-26
**Prerequisite:** ADMESH-DOMAINS-IO-SPEC.md (Specify phase, ambiguity 0.098)

---

## CLARIFY ROUND 1 — Boundary Keeper

### Finding C-1: `write_fort14()` ALSO crashes on quads (second critical bug)

**Evidence:** `src/chilmesh/CHILmesh.py:990`
```python
for i, tri in enumerate(elements, start=1):
    n1, n2, n3 = tri + 1   # ← crashes if tri has 4 values
    f.write(f"{i} 3 {n1} {n2} {n3}\n")
```

**Impact:** After fixing reader (I-1), writing quad mesh back raises `ValueError: too many values to unpack`.

**Clarification:** I-1 (fix reader) **must be paired with fixing writer** — atomic fix. Issue #40 expanded to cover both.

---

### Finding C-2: No quad/mixed test fixture exists

All 4 bundled fixtures have `3` in every element line. No way to regression-test quad/mixed code. Issue #40 unverifiable without quad fixture.

**Recommendation:** Create minimal synthetic 2×2 quad grid committed to `src/chilmesh/data/`. Real-world quad from ADMESH-Domains added later as integration test.

---

### Finding C-3: Coordinates are Cartesian, NOT geographic

Fixture coordinates (structured: x in [-2.0,1.0], Block_O: x in [5.0,10.0]) are clearly Cartesian. Specify phase assumed `x=lon, y=lat` — wrong.

**Clarification:** `bounding_box` returns `{min_x, max_x, min_y, max_y}`. Coordinate interpretation is caller's responsibility.

**ADMESH-Domains schema note:** Schema uses lat/lon — mismatch; raise as cross-repo question (Q-2).

---

### Finding C-4: ADCIRC_GRD is likely identical to ADCIRC format

All bundled ADCIRC-style files use identical line format regardless of extension. `.grd` is informal alias for `.14`.

**Clarification:** Route `ADCIRC_GRD` to `read_from_fort14()` in `from_admesh_domain()`. No spec change needed — routing in I-5 already defaults to `read_from_fort14()` for non-SMS_2DM types.

---

## CLARIFY ROUND 2 — Failure Analyst

### F-1: Quad write crashes on old code paths

After I-1 fixes reader, `write_to_fort14()` on quad mesh → `ValueError: too many values to unpack (expected 3)`.

**Prevention:** Fix writer in same PR as reader. Add `test_fort14_quad_roundtrip`.

### F-2: `bounding_box` for meshes with z ≠ 0

No failure — z excluded correctly (only `points[:,0]` and `points[:,1]` used). Document in docstring.

### F-3: `from_admesh_domain()` called before file downloaded

ADMESH-Domains uses HuggingFace on-demand download. `FileNotFoundError` if file absent.

**Prevention:** Check `Path(filename).exists()` in `from_admesh_domain()`; raise with guidance: *"File not found. Call mesh_record.load() first."*

### F-4: `compute_layers=False` then user calls `get_layer()`

With `n_layers=0`, `get_layer(0)` triggers `ValueError("out of range [0, -1]")` — confusing.

**Prevention:**
```python
if self.n_layers == 0:
    raise RuntimeError(
        "Layers not computed. Call _skeletonize() or "
        "re-initialize with compute_layers=True."
    )
```

---

## CLARIFY ROUND 3 — Seed Closer

### Q-1: Quad test fixture — synthetic or real-world?

**Decision:** Option C — synthetic (2×2 quad grid) for #40 unit test + ADMESH-Domains mesh for integration test (separate issue).

### Q-2: `bounding_box` key names — x/y vs lon/lat?

**Decision:** Option A — CHILmesh uses `{min_x, max_x, min_y, max_y}` (coordinate-agnostic). ADMESH-Domains maps `x→lon`, `y→lat` in its adapter.

### Q-3: Does `admesh-domains` become optional dep?

**Decision:** Option A — no import ever. Pure duck-typing via `getattr(record, "filename", None)`. No `extras_require`.

---

## UPDATED AMBIGUITY SCORES (Post-Clarify)

```
Goal Clarity:        0.95  ✓
Boundary Clarity:    0.92  ✓
Constraint Clarity:  0.93  ✓
Acceptance Criteria: 0.93  ✓

Ambiguity: 1.0 − (0.35×0.95 + 0.25×0.92 + 0.20×0.93 + 0.20×0.93) = 0.0655  ← gate ≤ 0.20 ✓
```

**Gate passed. Spec + Clarify complete.**
