# Visual Quality Pre-Merge Checklist

**Feature**: README Hero GIF Refinement (spec-002)  
**Input**: spec.md (FR-001…FR-009, SC-001…SC-007, user stories), plan.md (Architecture A–D, Test Plan), tasks.md (T001–T027, 8 phases)  
**When to use**: After implementation of all tasks is complete, before commit + push. Walk through each section; every item must be green before merge.

This checklist is organized around the six operator-identified defects (FR-001…FR-006) plus cross-cutting concerns (FR-007…FR-009, frame budget, determinism, process). Each item is an atomic verification step with a concrete command or artifact to inspect.

---

## 1. Solver Change Safety (FR-001 / FR-002)

The solver kwargs are additive, default-off, and cause zero regression in the default path.

- [ ] **CHK-101** Inspect `src/chilmesh/_vendor_admesh_truss.py` line 75: three new keyword-only parameters exist immediately after `history_every: int = 1,` with correct names and defaults:
  - `snapshot_retriangulate: bool = False`
  - `snapshot_strict_interior: bool = False`
  - `snapshot_hole_tau: float = 0.02`

- [ ] **CHK-102** Verify docstring: each parameter documented in the `distmesh2d_warmstart` docstring as "capture-side-only" and "default-off" (Constitution VIII).

- [ ] **CHK-103** Inspect capture block (lines 285–286): the branching logic routes to `Delaunay(p)` + filtering only when `snapshot_retriangulate=True`; else unchanged `t.copy()` default path.

- [ ] **CHK-104** Run default-path byte-identity test (from T002–T004):
  ```bash
  MPLBACKEND=Agg pytest tests/test_admesh_warmstart.py::TestHistoryCapture::test_history_none_is_default_noop -v
  MPLBACKEND=Agg pytest tests/test_admesh_warmstart.py::TestHistoryCapture -v -k "default_path"
  ```
  Assert all pass, no pre-existing tests modified.

- [ ] **CHK-105** Run opt-in behavior test (from T002–T004):
  ```bash
  MPLBACKEND=Agg pytest tests/test_admesh_warmstart.py::TestSnapshotRetriangulate -v
  ```
  Assert every captured `(p_i, t_i)` validates: centroid inside domain, zero inverted triangles, boundary rows pinned.

- [ ] **CHK-106** Run return-value invariance test (from T002–T004):
  ```bash
  MPLBACKEND=Agg pytest tests/test_admesh_warmstart.py::TestSnapshotRetriangulate::test_return_value_invariance -v
  ```
  Assert `(p_out, t_out)` with `snapshot_retriangulate=True` equals the hook-off run.

- [ ] **CHK-107** Verify no public API touch: `src/chilmesh/__init__.py` imports unchanged; `_vendor_admesh_truss` remains underscore-prefixed (not re-exported).

---

## 2. Truss Visual Quality (SC-001)

No truss-stage frame contains massively degenerate elements — hole-spanning, oversized, or inverted triangles.

- [ ] **CHK-201** Regenerate the GIF with the solver kwargs enabled:
  ```bash
  MPLBACKEND=Agg python scripts/generate_hero_animation.py
  ```
  Runtime ~2 minutes; GIF written to `docs/gallery/readme_pipeline_annulus.gif`.

- [ ] **CHK-202** Run the degenerate-element audit over all captured snapshots:
  ```bash
  python scripts/audit_hero_gif.py --check degenerate
  ```
  Assert output shows: `TOTAL hole=0, oversize=0, inverted=0` and `MAX-per-snapshot = 0`.

- [ ] **CHK-203** Verify corrected metrics in `audit_hero_gif.py::audit_degenerate`:
  - `hole-span = edge-midpoint fd > +0.02` (NOT `> -geps`)
  - `oversize = triangle area > 8× per-snapshot-median` (NOT 4×)
  - `inverted = signed-area sign flip vs. snapshot majority`

- [ ] **CHK-204** Frame-by-frame visual spot-check (manual): open the regenerated GIF in an image viewer; step through the truss-stage frames (roughly frames 78–145 out of ~265 total); confirm no element visually spans the annulus hole or overlaps neighbors at a grotesque scale.

- [ ] **CHK-205** Verify node positions remain real solver iterates: in `scripts/generate_hero_animation.py::_stage_data()`, the returned `hist` contains only real captured `(p, t)` pairs from `distmesh2d_warmstart`; no interpolated positions added.

---

## 3. Truss Pacing Uniformity (SC-002 / FR-003)

Truss-stage playback speed is approximately constant — no visible acceleration toward convergence. Thirds displacement ratio ≤ 3.

- [ ] **CHK-301** Inspect `scripts/generate_hero_animation.py::_stage_data()`: the snapshot selection uses equal-motion quantile resampling (plan.md § Architecture B), not fixed indices:
  - Compute per-iteration mean node displacement `d_iter`
  - Cumulative sum `Dcum`
  - Select `N_SNAP=32` at equal quantiles via `np.searchsorted`
  - Zero-motion fallback and min-count floor present

- [ ] **CHK-302** Run the pacing audit:
  ```bash
  python scripts/audit_hero_gif.py --check pacing
  ```
  Assert output shows: `thirds_displacement_ratio = X.XXX` where `X.XXX ≤ 3.0` (target ~1.26 per research.md Probe 2).

- [ ] **CHK-303** Verify truss-stage frame count and hold: inspect `scripts/generate_hero_animation.py::main()`:
  - Truss hold set to `TRUSS_HOLD = 2` frames per snapshot
  - Total truss frames = 3 preroll seed + `TRUSS_HOLD × len(snapshots)` (roughly 50–80 frames / 5–8 s)

- [ ] **CHK-304** Visual pacing check (manual): step through the rendered truss frames; perceived mesh relaxation speed should feel steady from start to finish, not accelerating toward the end.

---

## 4. Stage Choreography (SC-003 / FR-004 / FR-005)

Every stage boundary (Seed→Truss, Truss→FEM, FEM→Peel) and pre-loop-end hold shows a pause ≥ 1.5 s (≥ 15 frames @ 10 fps). Peel reveals layers boundary-inward with ≥ 0.8 s between reveals; final hold ≥ 2 s; no re-peel.

- [ ] **CHK-401** Inspect `scripts/generate_hero_animation.py::main()` static hold durations (F-13 standardizes the three inter-stage holds to 18 f / 1.8 s):
  - Seed→Truss hold: 18 frames (1.8 s)
  - Truss→FEM hold: 18 frames (1.8 s)
  - FEM→Peel hold (k=0): 18 frames (1.8 s) — from T014's peel schedule
  - Peel final hold (k=n_layers): ≥ 20 frames (≥ 2.0 s)

- [ ] **CHK-402** Verify peel layer-reveal schedule in `main()` from T014:
  - k=0 (start, pre-reveal): 18 frames / 1.8 s (F-13)
  - k=1, k=2, k=3 (layer reveals): 8 frames / 0.8 s each
  - k=n_layers (final hold): 20 frames / 2.0 s
  - Total peel stage: ~62 frames

- [ ] **CHK-403** Confirm layer reveal order: `_stage_data()` includes rim-distance monotonicity assert verifying boundary-inward ordering (layer 0 is outermost, strictly increasing mean-rim-distance per layer).

- [ ] **CHK-404** Verify no re-peel: the frame schedule ends on `k=n_layers` hold and does not loop back to k=0.

- [ ] **CHK-405** Manual stage-boundary audit: render the GIF and count held (pixel-identical) frames at each of the 4 boundaries in the `play` schedule; each must show ≥ 15 frames held.

---

## 5. Histogram Integrity (SC-004 / FR-006)

During the peel stage, the quality histogram converts to stacked-bar composition layer-by-layer. Per-bin totals remain constant; start state (k=0) is byte-identical to the FEM-hold appearance (convert-in-place).

- [ ] **CHK-501** Inspect `scripts/generate_hero_animation.py` new functions from T013:
  - `_layer_color(li, n_layers)` returns single-layer viridis color
  - `_peel_facecolors(q, elem_layer, n_layers, k)` starts from quality-colored `cool_r(q)`, overwrites `elem_layer < k` to layer colors
  - `_draw_peel_hist(ax, q, elem_layer, n_layers, k, D)` renders stacked bars: bottom-up layers 0..k-1 in viridis, remainder per bin in `cool_r(Normalize(0,1)(bin_midpoints))` pinned EXACTLY to `_draw_hist`'s `bar_colors` (F-04); per-element binning via `np.clip(np.floor(q*HBINS).astype(int), 0, HBINS-1)` (F-10)

- [ ] **CHK-502** Verify convert-in-place semantics: at k=0 (entry to peel stage), `_draw_peel_hist` renders the histogram **pixel-for-pixel identical** to the final FEM-stage frame's `_draw_hist(q_fem)` — the T015 audit renders both to the Agg `fig.canvas` buffer and asserts `np.array_equal` (F-04). Quality-colored bars, no layer colors yet.

- [ ] **CHK-503** Run the peel-invariance audit:
  ```bash
  python scripts/audit_hero_gif.py --check peel-invariance
  ```
  Assert output shows: `per-bin-totals-invariant = True` (per-bin sum constant across k=0..n_layers, matching the plain quality histogram — bins computed **independently** via `np.histogram`, not reused from the render path, F-10) **and** `k0-pixel-equal = True` (the k=0 vs FEM-hold canvas-buffer `np.array_equal` proof, F-04).

- [ ] **CHK-504** Manual per-frame stacked-bar inspection: step through the peel-stage frames in the rendered GIF; visually confirm each histogram bar's total height stays constant as layer colors appear from bottom-up in each bar.

- [ ] **CHK-505** Verify HBINS constant: `HBINS = 40` defined in the generator (matches research.md Probe 3's setup).

---

## 6. Metric Annotations (SC-005 / FR-007)

In every histogram-bearing frame (all 4 stages), Median is marked by a green dotted line + green text, Mean by a red dotted line + red text, Min stays neutral.

- [ ] **CHK-601** Inspect `scripts/generate_hero_animation.py` color constants from T019 (F-08 — exactly ONE new constant):
  - Median: reuse existing `GOOD = "#2ecc71"` (do NOT add a `GREEN` duplicate)
  - Mean: `RED = "#ff5555"` (the one new constant)
  - Min: reuse existing `TEXT = "#e8e8ee"` (do NOT add a `NEUTRAL` duplicate)

- [ ] **CHK-602** Inspect `_draw_metrics(ax, q, alpha=1.0)` function from T019:
  - Green dotted vertical line at median: `axvline(med, color=GOOD, linestyle=":", ...)`
  - Red dotted vertical line at mean: `axvline(mean, color=RED, linestyle=":", ...)`
  - Metric texts rendered at axes-fraction x=0.02 (Median), x=0.40 (Mean), x=0.76 (Min), y=1.06
  - Median text `GOOD`-green, Mean text `RED`, Min text `TEXT`-neutral
  - Every line/text scales by the `alpha` arg (F-09) so `f_wire` fades the metrics in with its bars
  - zorder=20 to keep lines/texts individually distinguishable

- [ ] **CHK-603** Verify `_draw_metrics` is called from `_draw_hist()` (replacing the old single-line `set_title` + `axvline`), which reaches the **seed/truss/FEM** stages. **F-01 — the peel-reveal frames call `_draw_peel_hist()`, NOT `_draw_hist()`**, so verify the peel frame function *also* calls `_draw_metrics(ax_hist, q_fem)` explicitly, immediately after `_draw_peel_hist(...)`. Confirm the metrics are NOT silently absent from peel frames.

- [ ] **CHK-604** Metric-text color inspection: render still-frame PNGs from each of the 4 stages (one from seed/wire, one from truss, one from FEM, and **specifically a peel `k=2` frame**); confirm green dotted line + green text for Median, red line + red text for Mean, neutral Min text in every histogram frame. For the peel `k=2` PNG, assert programmatically that the green + red vertical lines and the three colored metric texts are present (F-01 — metrics must render in peel k-frames, not only seed/truss/FEM).

- [ ] **CHK-605** Edge case: verify in the peel stage where median ≈ mean (roughly 2.4 bins apart), both lines and texts remain individually distinguishable (distinct colors + zorder prevent overlap).

---

## 7. Output Gates (SC-006 / SC-007 / FR-008 / FR-009)

The regenerated GIF is ≤ 4.2 MB, 864×396 px, renders on GitHub after push, and is deterministic run-to-run.

- [ ] **CHK-701** Check file size and dimensions:
  ```bash
  ls -lh docs/gallery/readme_pipeline_annulus.gif
  python -c "from PIL import Image; print(Image.open('docs/gallery/readme_pipeline_annulus.gif').size)"
  ```
  Assert size ≤ 4.2 MB and dimensions = (864, 396).

- [ ] **CHK-702** Verify post-process re-encode ran: inspect `scripts/generate_hero_animation.py::main()` — after `anim.save(...)` call, `_optimize_gif(str(OUT_PATH), colors=256)` is invoked (from T021). Function uses adaptive palette + optimize=True.

- [ ] **CHK-703** Verify `_optimize_gif` implementation from T021:
  - Builds a single global ADAPTIVE palette from a sample of frames
  - Re-quantizes every frame against it with dither=Image.NONE
  - Re-saves with optimize=True, preserving per-frame durations and loop count

- [ ] **CHK-704** Run determinism audit:
  ```bash
  python scripts/audit_hero_gif.py --check determinism
  ```
  Assert output shows: `determinism_check = PASS` (double run produces same frame count and schedule sequence).

- [ ] **CHK-705** Manual GitHub README rendering: push the changes to `development` (T024), then navigate to the README on GitHub and visually confirm the GIF animates smoothly without corruption (SC-006 manual verification).

- [ ] **CHK-706** Verify RNG seed unchanged: `scripts/generate_hero_animation.py` uses fixed seed 11; seed is unchanged from the baseline run.

---

## 8. Full Regression Suite (FR-001…FR-009)

All functional and success criteria validated together. No pre-existing tests modified.

- [ ] **CHK-801** Run full combined audit:
  ```bash
  python scripts/audit_hero_gif.py --check all
  ```
  Assert output shows all four sub-checks PASS: `degenerate`, `pacing`, `peel-invariance`, `determinism`.

- [ ] **CHK-802** Run full test suite:
  ```bash
  MPLBACKEND=Agg pytest tests/test_admesh_warmstart.py tests/test_generate_hero_animation.py -v
  ```
  Assert all tests pass (both new and pre-existing, no modifications to existing tests).

- [ ] **CHK-803** Verify generator smoke test exists (T005): `tests/test_generate_hero_animation.py` runs `main()` (or a truncated version) headless under `MPLBACKEND=Agg` and asserts output GIF is non-empty.

- [ ] **CHK-804** Verify no core-module or public-API touch. CHILmesh has **no** `_stages/*.py` "13 locked stage modules" — that is sibling-repo ADMESH's architecture (F-02). Real guard: `git diff development -- src/chilmesh/CHILmesh.py src/chilmesh/__init__.py` MUST be empty (no core-topology / `_peel` / `_build_adjacencies` / public-export changes); under `src/`, only `src/chilmesh/_vendor_admesh_truss.py` (additive capture block + the three kwargs) may show a diff.

- [ ] **CHK-805** Verify generator docstring = contract (T026 / F-06 / Constitution VIII): the `scripts/generate_hero_animation.py` module docstring and inline comments describe the SHIPPED design — reveal-based peel stage, equal-motion quantile pacing, the `snapshot_retriangulate` render-vs-physics note that explicitly warns against "restoring" the stale captured/return `t` and states the FEM input is `hist[-1]`, and the `_optimize_gif` post-process. No stale "single held peel frame" / "tiered hold" / "converges within ~15 iterations, each snapshot HELD" language remains.

- [ ] **CHK-806** Verify the auditor's own tests exist (T027 / F-12): `tests/test_audit_hero_gif.py` builds synthetic fixtures — one inverted triangle, one hole-spanning triangle, one `q=1.0` element — and asserts `audit_degenerate` flags the first two and bins the `q=1.0` element into the last bin (`HBINS-1`). Guards the auditor against silently self-certifying.

---

## 9. Process & Binary Push

All work committed to `development` branch via git CLI only; no MCP file tools used for binary upload.

- [ ] **CHK-901** Verify branch: run `git rev-parse --abbrev-ref HEAD`; assert output is `development` (not `main`, `daily-maintenance`, or a `claude/*` branch).

- [ ] **CHK-902** Verify commit format: run `git log --oneline -5` on `development`; inspect the 5 most recent commits:
  - Each follows `<type>: <imperative summary>` where type ∈ {fix, feat, docs, chore, refactor, test}
  - No `wip`, `fixup!`, `squash!`, `tmp` prefixes on commits bound for main

- [ ] **CHK-903** Verify git CLI binary push (from T024): `git add docs/gallery/readme_pipeline_annulus.gif` plus script/test files, `git commit` (with fallback `git -c commit.gpgsign=false commit` if signing fails), `git push origin development`.

- [ ] **CHK-904** Confirm binary integrity on remote: fetch the raw URL of the pushed GIF from GitHub and verify magic bytes (GIF89a or GIF87a header) and non-empty file size via a local tool (e.g., `curl -sI`), not an MCP read.

- [ ] **CHK-905** Verify seed-stage optimization: **only the `f_bfade`/`f_fall`/`f_settle` frames** (~50, dots-only animation) contain **no** histogram panel, keeping them ~10.8 KB each; do not accidentally add the histogram there (guardrail from plan.md § D). NOTE (F-05): `f_wire` + the wire-hold **intentionally DO** carry the fading histogram (~33 KB) — that is by design and must NOT be "optimized" away.

---

## Summary

**Pass Criteria**: Every checkbox `CHK-###` is green. All 9 functional requirements (FR-001…FR-009) and 7 success criteria (SC-001…SC-007) verified:

| Requirement | Checkboxes | Status |
|---|---|---|
| FR-001 (solver-side degeneracy suppression, default-off) | CHK-101 through CHK-107 | ✓ (solver kwargs + tests) |
| FR-002 (real captured node positions) | CHK-205 | ✓ (positions untouched) |
| FR-003 (uniform truss pacing, ≤3× band) | CHK-301 through CHK-304 | ✓ (quantile resampling audit) |
| FR-004 (≥1.5 s inter-stage pauses) | CHK-401 through CHK-405 | ✓ (schedule + manual frame count) |
| FR-005 (peel boundary-inward reveal, ≥0.8 s/layer) | CHK-401, CHK-402 | ✓ (peel schedule) |
| FR-006 (lock-step histogram, in-place, totals invariant) | CHK-501 through CHK-505 | ✓ (stacked-bar invariance audit) |
| FR-007 (green Median, red Mean, neutral Min) | CHK-601 through CHK-605 | ✓ (color-keyed metrics + render) |
| FR-008 (≤4.2 MB, 864×396, renders on GitHub) | CHK-701 through CHK-705 | ✓ (size gate + README check) |
| FR-009 (deterministic, fixed seed) | CHK-706 | ✓ (determinism audit) |
| SC-001 (zero degenerate frames) | CHK-202, CHK-203, CHK-204 | ✓ |
| SC-002 (thirds pacing ratio ≤3) | CHK-302 | ✓ |
| SC-003 (≥15 f boundary pauses) | CHK-405 | ✓ |
| SC-004 (per-bin totals invariant) | CHK-503 | ✓ |
| SC-005 (color-keyed metric texts) | CHK-604, CHK-605 | ✓ |
| SC-006 (GitHub render + size) | CHK-705 | ✓ |
| SC-007 (determinism) | CHK-704 | ✓ |

**Ready to merge**: All checkboxes ✓, binary pushed to `development`, README verification (CHK-705) passed on the pushed branch.
