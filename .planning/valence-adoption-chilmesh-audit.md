# CHILmesh audit — Valence CHILmesh-adoption blockers (Valence #212)

**Date:** 2026-07-04
**Scope:** Determine, per the four Valence issues blocked on upstream CHILmesh, what is
already implemented on `main`/`development` vs still missing; recommend the release(s) that
unblock each. Read-only against Valence (`spec.md`, `tests/parity/`, `_chilmesh.py`).

> **Unblock reality check (applies to every row below):** an unblock is only *real* once
> Valence bumps its `chilmesh` pin **and** its corpus-wide byte-parity gate (`tests/parity/`,
> currently 41/41) still passes. That gate is not run here (Valence is read-only this session).

## Branch topology (established this session)

- `development` and `main` have **identical `src/` trees**; the 120-commit `main..development`
  gap is docs/specs/skills churn, not source. So every source finding below holds for both.
- `main` is `version = 1.2.2` and already contains `fort15_io.py` + `summary_io.py`, but there
  is **no `v1.2.2` tag** (tags present: `v1.0.0a1`, `v1.1.0`, `v1.2.1`; PyPI latest = `1.2.1`).
- fort.14 read/write is **byte-for-byte unchanged since `v1.2.1`** — the only commits touching
  `CHILmesh.py` since the tag are FEM-solver / `_skeletonize→_layerize` rename / sdf-smoothing /
  perf work. So the divergences Valence measured against `1.2.1` still hold on `development`.

## Task 1 — per-requirement status + file:line evidence

### Valence #214 · P2 fort.14 read — **NOT-STARTED** (4 real gaps; 1 non-gap)

| chilmesh sub-requirement | Status | Evidence |
|---|---|---|
| preserve fort.14 id columns (node + element) | ❌ NOT-STARTED | `CHILmesh.py:2632` nodes read `p[1:4]` (ignore `p[0]` node-id); `:2640` elements read `p[2:2+n]` (ignore `p[0]` elem-id) — positional-only, id columns discarded |
| honor the header element count | ✅ **already honored** (not a chilmesh gap) | `:2627` reads `n_elems` from the header; `:2637` loops `range(n_elems)`. The `structuredMesh1-4` "element-count" divergence is a **Valence-baseline artifact**, not a chilmesh defect: Valence `plan.md:537` — header `660 374` → chilmesh correctly yields 660 elements, Valence's frozen baseline yields 1 |
| preserve file winding | ❌ NOT-STARTED | constructor runs `_ensure_ccw_orientation` (`:362` → `:404`), which reorders CW elements to CCW (`Test_Case_4.2`: all 4622 elements winding-permuted) |
| keep mixed tri/quad element arity | ❌ NOT-STARTED | `:2647-2652` pads triangles to 4 columns `[v0,v1,v2,v0]`; `connectivity_list` then exposes padded quads (`Mixed_Test` → padded-`0` vertex → Valence `NormalizationError`) |
| parse NBOU barrier/weir columns | ❌ NOT-STARTED | `:2687-2688` reads only column 0 (node id) per boundary node → drops barrier/weir `back_nodes`/`heights`/`coeffs`; open-segment `ibtype` also dropped |

### Valence #216 · P5 fort.14 write — **NOT-STARTED** → **implemented this PR**

| chilmesh sub-requirement | Status (audit) | Evidence |
|---|---|---|
| emit trailing NOPE/NBOU even when 0/0 | ❌ was NOT-STARTED | `CHILmesh.py:2750` guarded `if …boundary_segments:` → block omitted when empty; test `test_fort14_boundary_types.py:87` pinned the omit contract |
| return a success signal | ❌ was NOT-STARTED (module fn) | module `write_fort14` `:2720 -> None`. (method `write_to_fort14 :1504` already returned `True`, but Valence delegates to the **module** fn — Valence `plan.md:573`) |

### Valence #217 · P6 quality kernels — **PARTIAL** → **implemented this PR**

| chilmesh sub-requirement | Status (audit) | Evidence |
|---|---|---|
| expose a true equiangle-skewness metric | ⚠️ was PARTIAL (inverted) | `quality.py:76-79` (tri) + `:140-143` (quad) computed `1 − Qeas` under the name `skew` (1 = ideal). The underlying `(θmax−60)/120`, `(60−θmin)/60` terms are exact, but no *raw* Qeas metric (0 = ideal) was exposed — Valence would need an inverting shim (`plan.md:550`) |
| align/document signed_area numerics | ❌ was NOT-STARTED (undocumented) | `CHILmesh.py:470-474` float64 shoelace, signed, includes ½; near-degenerate → tiny non-zero, not exact `0.0` → Valence's `== 0.0` HARD connectivity gate is fragile (`plan.md:550,582`). No documented contract |

### Valence #218 · P4 fort.15 read — **DONE-on-main** (code); release-gated

| chilmesh sub-requirement | Status | Evidence |
|---|---|---|
| `read_fort15` / `fort15_io` shipped | ✅ DONE-on-main | `fort15_io.py:52 read_fort15`, `:108 write_fort15`, `:33 Fort15`; exported `__init__.py:21`; test `tests/test_fort15_roundtrip.py`. Present on **both** `main` + `development`. **Unreleased** (no `v1.2.2` tag; PyPI = 1.2.1). Tracked as CHILmesh #242 (R1) |

## Task 2 — does a single release ship the done work?

- **#218:** yes — a single **`v1.2.2` tag off the current `main` commit** ships it. **No PR needed**:
  the code, `version = 1.2.2`, and CHANGELOG `## [1.2.2] — 2026-06-15` are already on `main`.
  R1 / CHILmesh #242 is a **tag + PyPI publish**, nothing more.
- **#216 / #217:** were **not** on `main`; this PR adds them to `development`. They ship only
  after the rolling `development → main` PR (#243) merges, then a subsequent tag.

## Task 3 — implemented here vs proposed

### Implemented on `development` this session (low-risk, tested)

- **#216 writer** — `write_fort14` (module fn, `CHILmesh.py:2720`) now **always** emits the
  NOPE/NBOU block (`0/0/0/0` when empty) and **returns `True`**. Test contract updated:
  `test_write_without_boundaries_omits_section` → `test_write_without_boundaries_emits_zero_block`,
  plus `test_write_fort14_module_fn_returns_true`.
  ⚠️ This is a deliberate **public-behavior change** to a shipped writer (operator-requested).
  Round-trips are unaffected (a `0/0/0/0` trailer reads back as empty `boundary_segments`).
- **#217 quality** — new **`equiangle_skewness`** metric (aliases `"equiangle skewness"`, `"eas"`)
  in `quality.element_quality` returning raw Qeas ∈ [0,1] (0 = ideal), the exact complement of
  `skew`. Additive — existing `skew`/`skewness` unchanged. `signed_area` **numeric contract
  documented** (signed, ½-factor, float64 near-degenerate ≠ exact `0.0`). New tests:
  `test_equiangle_skewness.py`, two `signed_area` contract tests.
- Verified: targeted tests green; full fast suite **1327 passed, 77 skipped, 0 failed**
  (`-k "not block_o"`).

### Proposed, NOT implemented — #214 (high-risk / large)

The fort.14 reader gaps cannot be fixed by tweaking `read_from_fort14`, because the
CCW-normalizing `CHILmesh` constructor is load-bearing for the mesh algorithms (skeletonize,
signed_area orientation, adjacency). Proposal — add a **byte-preserving raw reader**, mirroring
the existing `fort13_io` / `fort15_io` / `summary_io` precedent:

- **New module `src/chilmesh/fort14_io.py`**, `read_fort14_raw(path) -> Fort14Raw` dataclass that
  faithfully preserves: node-id + element-id columns; **file winding** (no CCW reorder);
  **true per-element arity** (no tri→quad padding — keep a ragged list or an explicit
  `elem_arity` array); the **full NBOU barrier/weir columns** (`back_nodes`/`heights`/`coeffs`)
  and open-segment `ibtype`.
- A separate `.to_mesh()` builds the CCW-normalized `CHILmesh` when algorithms need it — the
  faithful record and the algorithm-ready mesh are distinct objects. The existing
  `read_from_fort14` / constructor are untouched.
- **Tests:** byte-parity across a fixture corpus including CW-wound, mixed tri/quad, and
  barrier/weir meshes, plus non-contiguous id columns.
- **Gate:** Valence swaps `parse_fort14` to delegate, then re-runs its 41/41 byte-parity. Only
  then is #214 truly unblocked (Valence `plan.md:541,566-567` explicitly defers the swap until
  these upstream fixes land).

## Task 4 — recommended chilmesh releases

| Release | Contents | Unblocks | Notes |
|---|---|---|---|
| **v1.2.2** (tag now) | `fort15_io` + `summary_io` (already on `main`) | **Valence #218** (P4) | = R1 / CHILmesh #242. Tag the **current `main` commit** (it is version 1.2.2). Publish to PyPI. No code work. |
| **v1.2.3** (after #243 merges) | this PR: #216 writer + #217 `equiangle_skewness` + `signed_area` docs | **Valence #216, #217** (P5, P6) | Additive/low-risk; satisfies the operator's stated #216/#217 requirements. |
| **v1.3.0** (later) | byte-preserving `fort14_io` reader (#214 proposal) + optionally the U1–U4 upstream gap-fillers (CHILmesh #238–241) | **Valence #214** (P2) + Valence P7–P9 | Larger; needs Valence to re-verify 41/41 after the pin bump. Valence spec earmarks 1.3.0 for U1–U4 — either co-release, or renumber that bundle to 1.4.0. |

**Sequencing footgun:** tags point at commits, not at "current HEAD version". Tag `v1.2.2` at the
existing `main` commit **before or independent of** merging this PR — do not wait for `main` to
"be" 1.2.2 after #216/#217 land (by then `main` will read 1.2.3). Then tag `v1.2.3` at the
post-merge `main`.

## Alternatives flagged

- **#216/#217 versioning:** shipped here as **1.2.3** (patch; additive, non-breaking). If strict
  "new public feature ⇒ minor" is preferred, cut as **1.3.0** instead — operator's call. The
  `release-integrity` api-semver gate is warn-first and token-gated, so a patch bump for these
  additive changes will not hard-block.
- **#216 minimalism:** Valence `plan.md` Task 9 only needs `write_fort14` to write parseable
  bytes (it reads them back); the always-emit `0/0/0/0` block + `True` return are the operator's
  stricter, canonical-ADCIRC requirement, implemented here.
- **#217 without this release:** Valence *could* swap P6 on `v1.2.1` today using the inverted
  `skew` + `abs(signed_area)` with a verdict-preserving shim (`plan.md:582-583`). The
  `equiangle_skewness` metric shipped here removes the need for that shim.
