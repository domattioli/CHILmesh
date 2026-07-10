---

description: "Task list for README Hero GIF Refinement"
---

# Tasks: README Hero GIF Refinement

**Input**: Design documents from `.specify/specs/002-hero-gif-refinement/`
**Prerequisites**: [plan.md](./plan.md) (required), [spec.md](./spec.md) (required for user stories), [research.md](./research.md), [data-model.md](./data-model.md), [quickstart.md](./quickstart.md)

**Tests**: Explicitly requested for the solver-library change only (Constitution III test-first mandate, plan.md § Test Plan items 1–3) plus one generator smoke test (plan.md § Test Plan item 4, Constitution III "visualization is smoke-test only"). The remaining visual/behavioral acceptance criteria (SC-001…SC-007) are covered by **scripted audits** (`scripts/audit_hero_gif.py`, built incrementally per story) and the manual render/size/README checks in `quickstart.md` — not unit tests, per Constitution III.

**Organization**: Tasks are grouped by user story so each story is independently implementable and testable. Phase order follows **spec priority** (P1 stories first: US1, US2, US4; then P2: US3, US5) — this does **not** match the spec's US1–US5 read order, since US4 (P1) is a higher priority than US3/US5 (P2). Story labels preserve the spec's own numbering.

## Format: `[T00N] [P?] [US#?] Description`

- **[P]**: Can run in parallel (touches a different file than sibling tasks at the same dependency tier, no ordering dependency between them)
- **[US#]**: Which user story this task belongs to (US1, US2, US3, US4, US5 — matches spec.md numbering)
- Every task names its exact file path(s)

## Path Conventions (this feature)

- `src/chilmesh/_vendor_admesh_truss.py` — solver library (additive-only change; not the CHILmesh public API)
- `scripts/generate_hero_animation.py` — the hero-GIF generator (rewritten)
- `scripts/audit_hero_gif.py` — **NEW**, scripted SC audits, built incrementally
- `tests/test_admesh_warmstart.py` — extended (solver-kwarg regression tests)
- `tests/test_generate_hero_animation.py` — **NEW**, generator headless smoke test
- `docs/gallery/readme_pipeline_annulus.gif` — regenerated binary output

---

## Phase 1: Setup

**Purpose**: Confirm a clean starting point before any edits, and capture a pre-change baseline for later regression diffing (SC-006, SC-007).

- [ ] [T001] Confirm the dev venv is active and green, and record the pre-change baseline: `source .venv/bin/activate` (or `pip install -e ".[dev]"` if missing); run `MPLBACKEND=Agg pytest tests/test_admesh_warmstart.py -v` and confirm all pre-existing tests pass; note the current `docs/gallery/readme_pipeline_annulus.gif` file size (`ls -l`), dimensions (`python -c "from PIL import Image; print(Image.open('docs/gallery/readme_pipeline_annulus.gif').size)"`), and stored-frame count — these are the "current holds" / "current 3.8 MB" baselines FR-004 and FR-008 reference.

**Checkpoint**: Baseline captured; safe to start Foundational work.

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: The one piece of shared "library" infrastructure every later stage's frames render on top of — the additive solver kwargs that fix the truss-frame degeneracy at the source (FR-001/FR-002). Per Constitution III, written test-first. Sequencing matters here specifically because research.md Probe 2's pacing measurements (US2) must be re-validated against the **post-fix** solver dynamics (plan.md § Complexity Tracking: "Re-run Probe 2 after the solver change lands"), so this must land before US2's resampling is finalized.

**⚠️ CRITICAL**: No user-story work in Phase 3+ begins until T002–T004 are complete (T005 may run in parallel — see below).

- [ ] [T002] Write failing tests in `tests/test_admesh_warmstart.py` for the not-yet-existing solver kwargs, per plan.md § Test Plan items 1–3:
  1. Extend `TestHistoryCapture` (class starts at the existing `test_history_none_is_default_noop`, line ~470) with a default-path byte-identity assert: call `distmesh2d_warmstart` twice on identical inputs — once with no history hook, once with `history_out=[]` and the new kwargs left at their defaults (and once more with `snapshot_retriangulate=False` passed explicitly) — assert `p`/`t` equal in all cases.
  2. Add a new `TestSnapshotRetriangulate` class: with `snapshot_retriangulate=True, history_every=1` on a small seeded domain, assert every captured `(p_i, t_i)` has `t_i` a valid Delaunay of `p_i` (each triangle centroid satisfies `fd(centroid) < -geps`), assert **zero** inverted triangles (signed-area sign matches the snapshot majority), and assert boundary rows stay pinned (mirrors the existing capture-test pattern).
  3. In the same class, assert return-value invariance: the `(p_out, t_out)` returned with `snapshot_retriangulate=True` equals the return with the hook off entirely (physics unchanged) — this is Probe 1's `p_identical/t_identical` result as a regression gate.
  Run `MPLBACKEND=Agg pytest tests/test_admesh_warmstart.py -v` and confirm the new tests **FAIL** (kwargs don't exist yet — `TypeError: unexpected keyword argument`).

- [ ] [T003] Implement the three additive keyword-only kwargs on `distmesh2d_warmstart` in `src/chilmesh/_vendor_admesh_truss.py`: insert `snapshot_retriangulate: bool = False,`, `snapshot_strict_interior: bool = False,`, `snapshot_hole_tau: float = 0.02,` immediately after the existing `history_every: int = 1,` (line 75). Replace **only** the capture block at lines 285–286 (`if history_out is not None and ...: history_out.append((p.copy(), t.copy()))`) with the branching logic from plan.md § Implementation Architecture A: when `snapshot_retriangulate`, run a fresh `Delaunay(p)` (already-imported `Delaunay`), keep triangles whose centroid satisfies `fd(centroid) < -geps` (already-defined `fd`, `geps`), and when `snapshot_strict_interior` also require all three edge-midpoints to satisfy `fd(midpoint) < snapshot_hole_tau`; otherwise keep the unchanged `t.copy()` default path. Add a Parameters docstring entry for each new kwarg stating default-off + capture-side-only semantics (Constitution VIII). Do not touch any other line. Depends on T002 (tests must exist and fail first).

- [ ] [T004] Run `MPLBACKEND=Agg pytest tests/test_admesh_warmstart.py -v`; confirm every test in the file passes, including the new T002 tests, **and** every pre-existing test passes unmodified (backward-compat gate — no existing test may be edited to make this pass). Depends on T003.

- [ ] [T005] [P] Add `tests/test_generate_hero_animation.py` (**NEW** file) with one fast headless smoke test per plan.md § Test Plan item 4 / Constitution III ("visualization is smoke-test only"): monkeypatch or parametrize `scripts/generate_hero_animation.py` to a low frame/iteration count, run its `main()` (or `_stage_data()` + a truncated schedule) under `MPLBACKEND=Agg`, and assert it completes headless and writes a non-empty GIF file. This is a regression safety net for the whole generator, independent of any specific story's visual behavior — does not depend on T002–T004 and may be written in parallel with them.

**Checkpoint**: Foundational ready — solver kwargs shipped, tested, and documented; the generator has a smoke-test safety net. User-story implementation can now begin.

---

## Phase 3: User Story 1 - Clean truss playback (Priority: P1) 🎯 MVP

**Goal**: No truss-stage frame renders a hole-spanning, grossly-oversized, or overlapping ("sail") element; every rendered node position remains a real captured solver iterate.

**Independent Test**: Regenerate the GIF; audit every one of the ~200 captured truss snapshots (not just the sampled/rendered frames) with the corrected degeneracy metrics; assert zero degenerate elements in any snapshot.

### Implementation for User Story 1

- [ ] [T006] [US1] In `scripts/generate_hero_animation.py::_stage_data()`: **(a)** pass `snapshot_retriangulate=True, snapshot_strict_interior=True` (leave `snapshot_hole_tau` at its 0.02 default) to the existing `distmesh2d_warmstart(...)` call (currently lines 108–112). **(b) F-03 — source the converged-truss frame + FEM input from `hist[-1]`, never from the solver return.** The return `(p_out, t_out)` keeps the pre-move stale `t` that `snapshot_retriangulate` does not fix, so: do NOT append `(p_out, t_out)` to `snapshots` (remove the `+ [(p_out, t_out)]`; `sel[-1] == len(hist)-1` already selects the final iterate); set `p_fin, t_fin = hist[-1]`; compute `q_truss = _quality(p_fin, t_fin)`; build the FEM input as `Mesh(connectivity=np.c_[t_fin, t_fin[:, 0]], points=np.c_[p_fin, np.zeros(len(p_fin))], ...)` and compute `q_fem` on `t_fin`. **(c) F-11 — compute the seed triangulation `T_seed` as a fresh `Delaunay(P_seed)` filtered by the centroid-in-domain keep test, NOT `snapshots[0][1]`** (which under `snapshot_retriangulate=True` is post-iteration-0-move), so the seed stage shows the true pre-solver scattered state; `q_seed = _quality(P_seed, T_seed)`. Depends on T003 (kwargs must exist on the solver).

- [ ] [T007] [P] [US1] Create `scripts/audit_hero_gif.py` (**NEW** file) with `audit_degenerate(hist, hole_tau=0.02, oversize_mult=8.0) -> dict` implementing the **corrected** metrics from plan.md § Implementation Architecture A / research.md Probe 1: `hole-span = edge-midpoint fd > +hole_tau` (NOT `fd > -geps`, which false-positives on legitimate boundary tris), `oversize = triangle area > oversize_mult × per-snapshot-median-area` (NOT 4×, since a clean annulus mesh's natural max/median area spread is ~4.3×), `inverted = signed-area sign flip vs. that snapshot's majority sign`. Return per-snapshot and total counts. Add a `--check degenerate` CLI entry point that imports `_stage_data()` and runs the audit over every captured snapshot in `hist`, **and additionally asserts by name that the FEM-input connectivity `hist[-1]` is degenerate-free** (F-03 — the FEM mesh is built from `hist[-1]`, NOT the solver return `(p_out, t_out)`; `hist[-1]` is already inside the 0/200 sweep but is asserted explicitly so a FEM-path regression can't slip through). Independent of T006 (different file) — may run in parallel.

- [ ] [T008] [US1] Render gate + verify: `MPLBACKEND=Agg .venv/bin/python scripts/generate_hero_animation.py` (~2 min), then `.venv/bin/python scripts/audit_hero_gif.py --check degenerate`; assert TOTAL hole=0, oversize=0, inverted=0 across all captured snapshots, MAX-per-snapshot = 0, **and** the explicit `hist[-1]` FEM-input connectivity assertion passes (F-03) (SC-001). Depends on T006, T007.

**Checkpoint**: US1 independently functional and testable — the truss stage never shows a degenerate element, and node positions are still real solver iterates (FR-002 untouched by this story).

---

## Phase 4: User Story 2 - Uniform convergence pacing (Priority: P1)

**Goal**: Truss-stage playback speed reads as approximately constant from start to finish — no visible acceleration toward convergence.

**Independent Test**: Measure mean per-rendered-frame node displacement across the truss stage; the ratio between the slowest and fastest thirds must be ≤ 3×.

### Implementation for User Story 2

- [ ] [T009] [US2] In `scripts/generate_hero_animation.py::_stage_data()`, replace the fixed-index `snap_idx` selection (currently `list(range(0, 16)) + list(range(16, 40, 2)) + list(range(40, len(hist), 15))` at line 119) with **equal-motion quantile resampling** per plan.md § Implementation Architecture B: compute per-iteration mean node displacement `d_iter`, its cumulative sum `Dcum`, then select `N_SNAP=32` snapshots at equal quantiles of `Dcum` via `np.searchsorted`, with a zero-motion uniform fallback and a min-count floor that tops up via uniform indices if `np.unique()` collapses selections. Depends on T006 (same function, sequential edit).

- [ ] [T010] [US2] In `scripts/generate_hero_animation.py::main()`, replace the tiered-hold `play` list build (currently `play += [snap] * (3 if k < 12 else 2 if k < 24 else 1)`, lines ~245–248) with a **uniform hold**: `TRUSS_HOLD = 2`, truss stage = 3 preroll seed frames + `TRUSS_HOLD × len(snapshots)`. Depends on T009 (same file/region, sequential).

- [ ] [T011] [P] [US2] Add `audit_pacing(hist, snap_idx, hold) -> float` to `scripts/audit_hero_gif.py` computing the thirds mean-per-rendered-frame-displacement ratio over the truss playback schedule (mirrors research.md Probe 2's measurement). Add a `--check pacing` CLI entry. Independent of T009/T010 (different file) — may run in parallel.

- [ ] [T012] [US2] Render gate + verify: regenerate the GIF, then `.venv/bin/python scripts/audit_hero_gif.py --check pacing`; assert the thirds displacement ratio ≤ 3 (target ≈ 1.26 per Probe 2), and confirm the truss stage still occupies 50–80 frames / 5–8 s (edge case: minimum-legible-duration). Depends on T009, T010, T011.

**Checkpoint**: US1 + US2 both independently functional — clean, evenly-paced truss playback.

---

## Phase 5: User Story 4 - Layer-by-layer peel reveal with live histogram (Priority: P1)

**Goal**: The peel stage starts from the final smoothed mesh and reveals `peel_layers()` layers one at a time, boundary-inward, with the quality histogram converting to matching stacked-bar segments in lock-step; per-bin totals never change.

**Independent Test**: Step through the peel-stage frames: layers recolor one at a time starting at the domain boundary and moving inward with a visible pause between each; the histogram's per-bin totals are bit-identical across every peel-stage frame while segment colors convert.

### Implementation for User Story 4

- [ ] [T013] [US4] In `scripts/generate_hero_animation.py`, add three new functions per plan.md § Implementation Architecture C: `_layer_color(li, n_layers)` (single-layer `viridis` color via `Normalize(0, n_layers - 1)`), `_peel_facecolors(q, elem_layer, n_layers, k)` (mesh facecolors: start from `cool_r(q)`, overwrite elements with `elem_layer < k` to their `_layer_color`), and `_draw_peel_hist(ax, q, elem_layer, n_layers, k, D)` (bottom-up stacked bars over the shared `HBINS=40` edges: per-element bin index via `np.clip(np.floor(q * HBINS).astype(int), 0, HBINS - 1)` so a `q=1.0` element lands in the last bin, F-10; revealed layers `0..k-1` each in `viridis(li)`; remainder — layers `>= k` — colored **per bin** by `cool_r(Normalize(0, 1)(bin_midpoints))` pinned EXACTLY equal to `_draw_hist`'s `bar_colors` (F-04), so at `k=0` the histogram is pixel-for-pixel identical to the FEM-hold histogram per FR-006's convert-in-place clarification; per-bin totals invariant across all `k`). Depends on T010 (same file, sequential edit after US2's changes land).

- [ ] [T014] [US4] In `scripts/generate_hero_animation.py::main()`, replace the single `add(40, f_peel_hold)` schedule entry (line 360) with the boundary-inward reveal schedule from plan.md § Implementation Architecture C: `k=0` hold **18 frames / 1.8 s** (F-13 — all inter-stage holds standardized to 18 f; this doubles as the FR-004 FEM→Peel inter-stage pause) → `k=1..n_layers-1` hold 8 frames / 0.8 s each → `k=n_layers` final hold 20 frames / 2 s (FR-005's held full-layer end; loop ends here, no re-peel). Recolor is **instant** at each reveal boundary (no fade), using `_peel_facecolors`/`_draw_peel_hist` from T013 and the boundary-inward layer ordering already computed into `D["elem_layer"]`/`D["n_layers"]` from `m.Layers`. In the same task, add a cheap rim-distance-monotonicity assert in `_stage_data()` (mean element distance to nearest domain rim must be strictly increasing from layer 0 to layer `n_layers - 1`) guarding this boundary-inward assumption against a future `chilmesh` layer-reorder (plan.md § Complexity Tracking risk: reveal order must NOT be verified by raw radius, which is non-monotonic on the annulus). Depends on T013.

- [ ] [T015] [P] [US4] Add `audit_peel_invariance(q_fem, elem_layer, n_layers) -> bool` to `scripts/audit_hero_gif.py`, asserting the per-bin stacked histogram totals are `np.array_equal` across `k=0..n_layers` and each equals a plain `np.histogram(q_fem, bins=HBINS, range=(0.0, 1.0))` computed **independently** — NOT by reusing the render binning (F-10), so a binning bug can't self-certify. **Additionally (F-04):** render `_draw_peel_hist(..., k=0)` and `_draw_hist(q_fem)` each to a pixel array via the Agg `fig.canvas` buffer (`np.asarray(fig.canvas.buffer_rgba())`) and assert `np.array_equal` on the two buffers — the byte-identical convert-in-place entry proof. Add a `--check peel-invariance` CLI entry. Independent of T013/T014 (different file) — may run in parallel.

- [ ] [T016] [US4] Render gate + verify: regenerate the GIF, then `.venv/bin/python scripts/audit_hero_gif.py --check peel-invariance`; step the rendered peel frames and confirm layers reveal boundary-first, ≥ 0.8 s (≥ 8 held frames) between consecutive reveals, and the stage ends on a held full-layer view (SC-004). Depends on T014, T015.

**Checkpoint**: US1 + US2 + US4 (all three P1 stories) independently functional — this is the MVP candidate (see Implementation Strategy below).

---

## Phase 6: User Story 3 - Distinct stage pauses (Priority: P2)

**Goal**: Every stage boundary (Seed→Truss, Truss→FEM, FEM→Peel) and the pre-loop-end hold pause for at least 1.5 s, noticeably longer than the pre-feature holds.

**Independent Test**: Count held (pixel-identical) frames at each of the 4 boundary points in the rendered GIF; each must meet ≥ 15 frames @ 10 fps.

### Implementation for User Story 3

- [ ] [T017] [US3] In `scripts/generate_hero_animation.py::main()`, set the Seed→Truss pre-truss static hold (currently `add(8, lambda i, n: f_wire(n - 1, n))`, 0.8 s) to **18 frames** and the Truss→FEM hold (currently `add(16, f_truss_hold)`, 1.6 s) to **18 frames** — F-13 standardizes ALL inter-stage holds to 18 f / 1.8 s, which clears the ≥ 15-frame / 1.5 s FR-004 floor and is noticeably longer than every pre-feature value. Note: the FEM→Peel boundary pause and the pre-loop-end hold are **already satisfied** by T014's peel schedule (`k=0` hold = **18 f**, `k=n_layers` final hold = 20 f) — do not re-implement those here. Depends on T014 (US4's peel schedule must exist so this task's verification in T018 can audit all 4 boundaries against a complete schedule).

- [ ] [T018] [US3] Add a schedule-introspection check (count consecutive identical-`(fn, i, n)`-effect runs at each of the 4 boundaries: Seed→Truss, Truss→FEM, FEM→Peel, pre-loop-end) and run it against the built `schedule` list; confirm each boundary holds ≥ 15 frames (SC-003). Render gate to confirm in the actual rendered GIF. Depends on T017.

**Checkpoint**: US1 + US2 + US3 + US4 independently functional.

---

## Phase 7: User Story 5 - Color-keyed metric annotations (Priority: P2)

**Goal**: In every histogram-bearing frame (all 4 stages), the Median is marked by a green dotted line + green text, the Mean by a red dotted line + red text, and the Min stays neutral text.

**Independent Test**: Inspect any histogram-bearing frame from any stage: the green-line metric's text is green, the red-line metric's text is red, the remaining metric's text is neutral.

### Implementation for User Story 5

- [ ] [T019] [US5] In `scripts/generate_hero_animation.py`: add the **one** new module constant `RED = "#ff5555"` (Mean); **reuse** the existing `GOOD = "#2ecc71"` for the green Median and the existing `TEXT = "#e8e8ee"` for the neutral Min (F-08 — do NOT add `GREEN`/`NEUTRAL` duplicates). Add `_draw_metrics(ax, q, alpha=1.0)` per plan.md § Implementation Architecture C: green (`GOOD`) dotted Median `axvline` + green `"Median: {:.3f}"` text, red (`RED`) dotted Mean `axvline` + red `"Mean: {:.3f}"` text, neutral (`TEXT`) `"Min: {:.3f}"` text — at axes-fraction `x = 0.02 / 0.40 / 0.76`, `y = 1.06`, `zorder=20` so the two lines/texts stay individually distinguishable even when median ≈ mean; every line/text scales by the `alpha` arg (F-09) so `f_wire` fades the metrics in with its bars. Call `_draw_metrics(ax_hist, q, alpha)` from `_draw_hist()` (lines ~196–208) **in place of** the current single `ax_hist.axvline(med, ...)` + `ax_hist.set_title(...)` block — this reaches the **seed/truss/FEM** stages, since `f_wire`, `f_truss`, `f_truss_hold`, `f_fem`, `f_fem_hold` all call `_draw_hist()`. **F-01 — the peel-reveal frames do NOT call `_draw_hist()`; they call `_draw_peel_hist()`** (T013), so the peel frame function (T014) must additionally call `_draw_metrics(ax_hist, q_fem)` explicitly, immediately after `_draw_peel_hist(...)` (peel `q = q_fem`, constant across `k`). The prior claim that "peel-reveal frames already call `_draw_hist()`" was FALSE and is corrected here. Depends on T017 (same file, sequential after US3's schedule edits).

- [ ] [T020] [US5] Render gate + visual frame-dump review: save PNG stills from each of the 4 stages (a seed/wire frame, a truss frame, a FEM frame, a peel frame) and confirm every one carries both color-keyed lines and matching colored texts, individually distinguishable in the peel stage where median/mean are ~2.4 bins apart (SC-005). Depends on T019.

**Checkpoint**: All 5 user stories independently functional and verified.

---

## Phase 8: Polish & Cross-Cutting Concerns

**Purpose**: Architecture D (frame/size budget) is cross-cutting by design — it applies to the whole regenerated GIF, not one story — plus the final combined verification and binary push.

- [ ] [T021] [P] In `scripts/generate_hero_animation.py`, add `_optimize_gif(path, colors=256)` per plan.md § Implementation Architecture D: open the saved GIF, build a single global `ADAPTIVE` palette from a sampled stack of frames, re-quantize every frame against it with `dither=Image.NONE`, and re-save with `optimize=True` and the original per-frame durations preserved; call it immediately after `anim.save(...)` in `main()` (line 372). In the same task, verify **only the `f_bfade`/`f_fall`/`f_settle` seed frames** (~50) stay dots-only with **no** histogram-panel content (research.md Probe 4 guardrail — an added histogram there costs +~1.1 MB). NOTE (F-05): `f_wire` + the wire-hold **intentionally carry the fading histogram** (~33 KB/frame, +~220 KB) — that is by design and is NOT part of this guardrail. Independent of T022 (different file) — may run in parallel.

- [ ] [T022] [P] Add `audit_determinism()` to `scripts/audit_hero_gif.py`: run `_stage_data()` plus the `main()` schedule-build twice under the fixed RNG seed, diff the resulting frame count and the `(fn.__name__, i, n)` schedule sequence for equality. Add a `--check determinism` CLI entry, and a `--check all` entry that runs `degenerate` + `pacing` + `peel-invariance` + `determinism` together (SC-007 / FR-009). Independent of T021 (different file) — may run in parallel.

- [ ] [T023] Full combined gate: render gate (`MPLBACKEND=Agg .venv/bin/python scripts/generate_hero_animation.py`, ~2 min) → size gate (`ls -l docs/gallery/readme_pipeline_annulus.gif` ≤ 4.2 MB; confirm dims still 864×396 via `python -c "from PIL import Image; print(Image.open('docs/gallery/readme_pipeline_annulus.gif').size)"`) → `.venv/bin/python scripts/audit_hero_gif.py --check all` (degenerate + pacing + peel-invariance + determinism, all green) → re-run `MPLBACKEND=Agg pytest tests/test_admesh_warmstart.py tests/test_generate_hero_animation.py tests/test_audit_hero_gif.py -v` for a final regression confirmation (all green, no pre-existing test modified). Depends on T021, T022, T026, T027, and the checkpoints of every prior story (T008, T012, T016, T018, T020).

- [ ] [T024] Binary push (**git CLI ONLY — never the GitHub MCP file tools, which double-encode base64 and corrupt binaries**; `development` branch ONLY): `git checkout development`; `git add docs/gallery/readme_pipeline_annulus.gif scripts/generate_hero_animation.py scripts/audit_hero_gif.py src/chilmesh/_vendor_admesh_truss.py tests/test_admesh_warmstart.py tests/test_generate_hero_animation.py tests/test_audit_hero_gif.py`; commit (`git -c commit.gpgsign=false commit` fallback if signing fails, per repo CLAUDE.md); `git push origin development`. Confirm the pushed binary is intact (non-empty, GIF magic bytes) via the raw URL, not an MCP read. Depends on T023 (all gates green).

- [ ] [T025] Manual README verification on the pushed `development` branch: confirm the regenerated GIF actually animates in GitHub's rendered README view (SC-006's manual leg per `quickstart.md`) — this must run **after** T024's push, never before. Depends on T024.

- [ ] [T026] [P] **(F-06 — docstring = contract, Constitution VIII)** Rewrite the `scripts/generate_hero_animation.py` module docstring (lines ~1–23) and the stale inline comments at lines ~116–118 (the "solver converges within ~15 iterations… each early snapshot is HELD" note — superseded by equal-motion quantile pacing) and ~245 (the "hold each early snapshot 3 frames, mid 2, late 1" tiered-hold note — now a uniform `TRUSS_HOLD`) so the module contract matches the shipped design: (1) **reveal-based peel stage** (boundary-inward `k`-layer reveal + lock-step stacked histogram), NOT "a single held frame, no reveal"; (2) **equal-motion quantile pacing** (not fixed-index HELD sampling); (3) a **`snapshot_retriangulate` render-vs-physics connectivity note** that explicitly WARNS future sessions NOT to "restore" the stale captured/return `t` (F-03) and states the FEM input comes from `hist[-1]`; (4) the mandatory **`_optimize_gif`** post-process. Comments/docstring only — no behavioral code change. **⚠️ This task edits a `.py` file's prose; per repo coding-dispatch policy it is a code edit and must be dispatched to a Haiku subagent, then verified.** Depends on T014, T021 (the shipped behaviors the docstring describes must exist first).

- [ ] [T027] [P] **(F-12 — guard the auditor)** Add `tests/test_audit_hero_gif.py` (**NEW** file): unit-test `scripts/audit_hero_gif.py`'s `audit_degenerate` (and its binning helper) against **hand-built synthetic fixtures** the auditor MUST flag / bin correctly — (a) one inverted triangle (signed-area sign flipped vs. majority → `inverted` count ≥ 1), (b) one hole-spanning triangle (edge-midpoint `fd > +hole_tau` → `hole` count ≥ 1), (c) one `q = 1.0` element (must land in the last histogram bin, index `HBINS-1`, not overflow). Guards the auditor so a broken audit can't silently self-certify the GIF (pairs F-10's independent-binning rule). Independent of T021/T022/T026 (different file) — may run in parallel. Depends on T007 (audit script must exist).

**Checkpoint**: Feature complete — all 9 functional requirements (FR-001…FR-009) and all 7 success criteria (SC-001…SC-007) verified; binary pushed to `development`.

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — starts immediately.
- **Foundational (Phase 2)**: Depends on Setup. T002→T003→T004 is a strict TDD chain; T005 (generator smoke test) is independent and may run in parallel with T002–T004. **Blocks all user stories.**
- **User Stories (Phase 3–7)**: All depend on Foundational completion. Phase order = spec priority (P1: US1→US2→US4, then P2: US3→US5), **not** spec read-order (US1..US5), since US4 is a higher priority than US3/US5.
- **Polish (Phase 8)**: Depends on all 5 user stories being complete.

### User Story Dependencies

- **US1 (P1)**: Depends only on Foundational (needs the solver kwargs from T003). No dependency on other stories.
- **US2 (P1)**: Code depends on US1's `_stage_data()` edit landing first (same function, sequential in the same file) — not a *behavioral* dependency (US2's resampling logic doesn't require US1's fix to be conceptually correct), but the two edits share code so must be sequenced.
- **US4 (P1)**: Code depends on US2 having landed (same file, sequential edits to `main()`/`_stage_data()`). Behaviorally independent of US1/US2's specific fixes.
- **US3 (P2)**: Its own two hold-length edits have no *behavioral* dependency on US4, but its **verification** (SC-003 must audit all 4 boundaries, including FEM→Peel and pre-loop-end) requires US4's peel schedule to exist — hence placed after US4.
- **US5 (P2)**: Single choke-point edit (`_draw_hist()`) reaches every stage already built by US1–US4; sequenced last among stories so its render/verify checkpoint reflects the final combined visual.

### Parallel Opportunities

- T005 (generator smoke test) run alongside T002–T004 (solver TDD chain) — different files, no shared dependency.
- Within each user-story phase, the `scripts/audit_hero_gif.py` task (T007, T011, T015) runs in parallel with that story's `scripts/generate_hero_animation.py` implementation task(s) — different files.
- T021 (`_optimize_gif` in the generator) and T022 (`audit_determinism` in the audit script) run in parallel in Polish — different files. T027 (`tests/test_audit_hero_gif.py`, F-12) is also parallel-safe there; T026 (F-06 docstring rewrite) touches the generator so it sequences after T021's `_optimize_gif` edit lands.

---

## Parallel Example: Foundational + User Story 1

```bash
# Can run simultaneously (different files, no shared dependency):
Task: "Add tests/test_generate_hero_animation.py generator smoke test" (T005)
Task: "Write failing solver-kwarg tests in tests/test_admesh_warmstart.py" (T002)

# Later, within US1 (after T003/T006 land):
Task: "Create scripts/audit_hero_gif.py with audit_degenerate()" (T007)
# ...runs independently of T006's edit to scripts/generate_hero_animation.py
```

---

## Implementation Strategy

### MVP First: US1 + US2 + US4 (all three P1 stories)

The spec marks **three** stories P1 — US1 (clean truss playback), US2 (uniform pacing), and US4 (peel reveal + live histogram) — not two; US3 and US5 are both P2. The MVP is therefore all three P1 stories together:

1. Complete Phase 1 (Setup) + Phase 2 (Foundational) — solver-kwarg fix + smoke test.
2. Complete Phase 3 (US1) → Phase 4 (US2) → Phase 5 (US4), in that order (each has a same-file sequential dependency on the previous, per Dependencies above).
3. **STOP and VALIDATE**: run T008, T012, T016's render-gate + audit checks independently; this is a shippable improvement even without US3/US5.
4. Add Phase 6 (US3) and Phase 7 (US5) as fast-follow — both are cheap, additive, single-choke-point edits.
5. Phase 8 (Polish) — the frame/size budget post-process and final combined gate apply once, after all 5 stories land, then binary push.

### Incremental Delivery

Each story phase ends on a render-gate + scripted-audit checkpoint (T008 / T012 / T016 / T018 / T020) that is independently verifiable without the later stories existing — so a session that runs out of budget after Phase 5 (US4) has already shipped 3 of the 6 operator-identified defects, verified end-to-end, and can safely stop at that checkpoint rather than mid-story.

---

## Notes

- [P] tasks touch a different file than their same-tier siblings and have no ordering dependency between them.
- [US#] labels map each task to its spec.md user story for traceability; tasks without a [US#] tag are Setup/Foundational/Polish (shared infrastructure, not story-specific).
- Every story phase ends on a render-gate + scripted-audit (or manual-review) checkpoint — do not consider a story "done" until that checkpoint task passes.
- No task in this list edits any core CHILmesh topology/skeletonization code or the public API. The concrete guard (F-02 — CHILmesh has **no** `_stages/*.py` "13 locked stage modules"; that is sibling-repo ADMESH's architecture): `git diff development -- src/chilmesh/CHILmesh.py src/chilmesh/__init__.py` MUST be empty (no core-topology / `_peel` / `_build_adjacencies` / public-export changes). Under `src/`, only `src/chilmesh/_vendor_admesh_truss.py` (the additive capture-block branch + the three keyword-only kwargs, T003) may change — it is underscore-prefixed and not re-exported. Every other task touches only the generator script, the new audit script, or test files.
- T024 (binary push) MUST precede T025 (README verification) — the verification target is the pushed `development` branch's rendered README, which doesn't exist until after the push.
- Commit after each task or logical group per repo convention; all writes to `docs/gallery/readme_pipeline_annulus.gif` go through the git CLI, never the GitHub MCP file tools (binary corruption risk).
