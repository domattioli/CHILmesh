# Implementation Plan: README Hero GIF Refinement

**Branch**: `development` (repo policy: single AI staging branch; no feature branch) | **Date**: 2026-07-10 | **Spec**: [`spec.md`](./spec.md)
**Input**: Feature specification from `.specify/specs/002-hero-gif-refinement/spec.md`

## Summary

Refine the README hero animation (`scripts/generate_hero_animation.py` →
`docs/gallery/readme_pipeline_annulus.gif`) along five operator-identified axes while
preserving the binding real-solver-playback mandate (#179): (1) suppress massively
degenerate/overlapping elements in truss-stage frames **at the solver source** via three
additive, default-off snapshot kwargs on `distmesh2d_warmstart`; (2) make truss pacing
perceptually uniform by resampling snapshots at equal quantiles of cumulative node
displacement with a uniform per-snapshot hold; (3) lengthen inter-stage pauses to ≥ 1.5 s;
(4) redesign the peel stage into a boundary-inward layer-by-layer reveal with a lock-step
stacked-bar histogram conversion; (5) color-key the histogram metric annotations (green
Median line + green text, red Mean line + red text, neutral Min). A mandatory
post-process GIF re-encode (single global adaptive palette + `optimize=True`) keeps the
output under the ~4.2 MB budget with no visual loss.

The work is validated by four completed probes that replicated `_stage_data()` exactly
(seed 11, annulus fixture, `deltat=0.02 Fscale=1.2 niter=200 history_every=1`) and
measured every number cited below. The technical approach is: **one additive solver change
(three keyword-only kwargs, default-off, capture-side only) + one rewritten generator
script**. No new runtime dependency. No change to the CHILmesh core
(`src/chilmesh/CHILmesh.py`) or the public API surface (`src/chilmesh/__init__.py`
exports); the only source file changing under `src/` is `_vendor_admesh_truss.py`
(capture block + additive kwargs).

## Technical Context

**Language/Version**: Python 3.10+ (repo floor; matches `pyproject.toml`).
**Primary Dependencies**: `numpy`, `scipy` (`scipy.spatial.Delaunay` — already imported in
the solver module), `matplotlib` (Agg backend, `animation.PillowWriter`), `Pillow`
(`PIL.Image` / `ImageSequence` — already a transitive matplotlib dependency; used
directly for the post-process GIF re-encode). **No new dependency is added.**
**Storage**: Files only — the generated binary GIF at
`docs/gallery/readme_pipeline_annulus.gif`; no database.
**Testing**: `pytest`. New solver-kwarg regression tests live in
`tests/test_admesh_warmstart.py::TestHistoryCapture` (extending the existing class), plus a
new `TestSnapshotRetriangulate` class. The generator itself is a matplotlib rendering
script (smoke-tested / manually verified per Constitution III — "visualization is
smoke-test only").
**Target Platform**: Linux (CI) / macOS / Windows; render is headless (`MPLBACKEND=Agg`).
**Project Type**: Single library + a repo-local generator script. Not web/mobile.
**Performance Goals**: Full GIF regeneration ≤ ~2 min locally (probe-measured within
budget with 200 extra Delaunay calls only on the hero run); default solver path byte-
identical and zero-cost when the new kwargs are unused.
**Constraints**: Output GIF ≤ 4.2 MB at fixed 864×396 px, 10 fps, single loop ending on
the peel result (FR-008); deterministic run-to-run under RNG seed 11 (FR-009 / SC-007);
truss node positions remain real captured solver iterations (FR-002); solver library
change additive and default-off, existing behavior byte-identical when unused (FR-001).
**Scale/Scope**: One ~380-node annulus mesh; ~200 captured solver iterations; ~110–120 unique
animation frames, ~265 total played frames after held-frame expansion (held frames are free —
see research Probe 4; the earlier "~180 unique" was an overcount, F-07).

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.* Gates evaluated
against `.specify/memory/constitution.md` (v1.1).

| Principle | Gate | Verdict |
|---|---|---|
| **III. Test-First + Regression-Gated** | New solver kwargs need TDD + a regression test that would have caught the staleness defect. | **PASS (with mandate).** The plan requires the solver-kwarg tests to be written before/with the implementation, mirroring the existing `TestHistoryCapture` pattern: (a) a default-path byte-identity assert (`history_out=[]` with the new kwargs absent → output identical to a no-hook run), and (b) an opt-in behavior assert (`snapshot_retriangulate=True` → every captured `t` is a valid Delaunay of its own `p`, zero inverted triangles). This directly regression-guards the one-Euler-step `(p_new, t_old)` staleness Probe 1 found. |
| **VII. Public API Stability** | Does the solver change break a public import or signature? | **PASS.** `distmesh2d_warmstart` is a `_vendor_admesh_truss` internal (underscore-prefixed module, not re-exported in `__init__.py`). The three new params are **keyword-only** (after the existing `*`), each defaults to the current behavior (`snapshot_retriangulate=False`, `snapshot_strict_interior=False`, `snapshot_hole_tau=0.02`), so no existing caller changes. Purely additive — same shape as the shipped `history_out=` / `history_every=` precedent. |
| **Backward compatibility / Zero-regression** | Default solver behavior byte-identical when new options unused (FR-001). | **PASS — probe-proven.** Probe 1 verified `p_identical=True, t_identical=True, maxabs dp=0.0` between the new-code default path and the original off-hook solver. Physics/forces/return path is untouched; the new branch executes only inside the `history_out is not None` capture block and only when `snapshot_retriangulate=True`. |
| **II. Mesh Immutability** | Does the change mutate mesh topology or the smoother? | **PASS.** No mesh mutation. The generator composes existing calls (`Mesh(...)`, `smooth_mesh(method="fem")`, `peel_layers`/`Layers`); the solver change only re-triangulates a **copy** of captured points for rendering. |
| **X / IV. Correctness over performance / MATLAB parity** | Does re-triangulating the captured snapshot alter solver semantics? | **PASS.** The captured `t` is a diagnostic/render artifact that the solver itself discards on the next iteration; FR-002 governs node **positions** (untouched). A fresh Delaunay of the real points is more faithful than the 1-step-stale connectivity, not less. Documented explicitly in the generator + solver docstrings so a future session does not "restore" the stale capture. |
| **VIII. Documentation = Contract** | New kwargs documented. | **PASS (mandate).** The solver docstring gains a Parameters entry for each new kwarg stating default-off + capture-side-only semantics; the generator gains a comment block explaining the render-vs-physics connectivity distinction. |

**Result: PASS.** No violations; Complexity Tracking is used only for the open risks
below, not for justifying a constitution deviation.

## Project Structure

### Documentation (this feature)

```text
.specify/specs/002-hero-gif-refinement/
├── plan.md              # This file
├── research.md          # Phase 0 — one section per probe (decisions + evidence)
├── data-model.md        # Phase 1 — entities + invariants
├── quickstart.md        # Phase 1 — regenerate + verify + SC checklist
└── tasks.md             # Phase 2 — NOT created here (/speckit-tasks)
```

No `contracts/` directory: this feature exposes no network/service API; the sole "contract"
is the additive solver-kwarg signature, captured inline in this plan and in the solver
docstring, and pinned by the new regression tests. **Agent-context update script skipped:**
no `.specify/scripts/bash/` agent-context updater exists in this repo; noted here per the
plan template rather than run.

### Source Code (repository root)

```text
src/chilmesh/
└── _vendor_admesh_truss.py      # MODIFIED — 3 additive keyword-only kwargs on
                                 #   distmesh2d_warmstart; capture-block branch only.
                                 #   Default path byte-identical.

scripts/
└── generate_hero_animation.py   # REWRITTEN — snapshot resampling (pacing), peel-reveal
                                 #   + stacked histogram, color-keyed metrics, longer
                                 #   stage pauses, post-process GIF optimize, passes the
                                 #   new solver kwargs on the hero run.

tests/
└── test_admesh_warmstart.py     # EXTENDED — TestHistoryCapture gains default-identity +
                                 #   opt-in-behavior asserts for the new kwargs
                                 #   (new class TestSnapshotRetriangulate).

docs/gallery/
└── readme_pipeline_annulus.gif  # REGENERATED binary output (git CLI push only, never MCP).
```

**Structure Decision**: Single-project library layout (`src/chilmesh/` + `tests/` +
`scripts/`), already in place. The feature touches exactly one library module (additive),
one script (rewrite), and one test file (extend) — no new package, directory, or module.

## Implementation Architecture

### A. Solver kwargs (from Probe 1 — degeneracy suppression)

Add three **keyword-only** parameters to `distmesh2d_warmstart`, inserted immediately after
the existing `history_every: int = 1,` (line 75), all defaulting to current behavior:

```python
    snapshot_retriangulate: bool = False,
    snapshot_strict_interior: bool = False,
    snapshot_hole_tau: float = 0.02,
```

Replace **only** the capture block at lines 285–286 (leave every other line untouched):

```python
if history_out is not None and (iteration % max(1, history_every) == 0):
    if snapshot_retriangulate:
        # ADDITIVE, default-OFF: capture a triangulation CONSISTENT with the
        # post-move points p. Fixes the 1-Euler-step (p_new, t_old) staleness that
        # renders inverted "overlapping sail" triangles. Physics/forces/return are
        # untouched -> convergence byte-identical to the default path. The captured
        # connectivity is a RENDER artifact; the solver discards its own t next
        # iteration anyway (do NOT "restore" the stale capture in a future edit).
        _tri = Delaunay(p)
        _ta = _tri.simplices
        _cen = p[_ta].mean(axis=1)
        _keep = fd(_cen) < -geps
        if snapshot_strict_interior:
            _p0, _p1, _p2 = p[_ta[:, 0]], p[_ta[:, 1]], p[_ta[:, 2]]
            _m01 = 0.5 * (_p0 + _p1); _m12 = 0.5 * (_p1 + _p2); _m20 = 0.5 * (_p2 + _p0)
            _keep = _keep & (fd(_m01) < snapshot_hole_tau) \
                          & (fd(_m12) < snapshot_hole_tau) \
                          & (fd(_m20) < snapshot_hole_tau)
        history_out.append((p.copy(), _ta[_keep].copy()))
    else:
        history_out.append((p.copy(), t.copy()))   # unchanged default path
```

Dependencies are already in scope: `Delaunay` (module import), `fd` (passed callable),
`geps` (defined line 146). No new imports. The hero run enables
`snapshot_retriangulate=True, snapshot_strict_interior=True` (tau left at 0.02). Probe 1
proved: over all 200 captured snapshots with corrected metrics, **hole=0, oversize=0,
inverted=0**, every-frame max 0; and the returned `(p_out, t_out)` is byte-identical to the
off-hook run. `snapshot_strict_interior` is a measured no-op on this annulus (no hole-
spanners) — the load-bearing lever is `snapshot_retriangulate` alone; strict is kept ON for
robustness on non-annulus domains.

**Return vs. capture (F-03).** The `snapshot_retriangulate` branch rewrites only the captured
`history_out` entries — the solver's *return* `(p_out, t_out)` still carries the pre-move stale
connectivity from its last Euler step. "Byte-identical to the off-hook run" is a *physics/node-
positions* invariance proof, **not** a promise that `t_out` is degeneracy-free. The generator
therefore sources both the converged-truss frame and the FEM input from `hist[-1]` (re-
triangulated), never from the return — see §B.

**Metric-definition correction (generator + spec audit):** any degenerate audit must use
`hole-span = edge-midpoint fd > +0.02` (NOT `> -geps`, which flags all legitimate boundary
tris) and `oversize = area > 8× snapshot-median` (NOT 4×, since a valid annulus mesh's
natural max/median area spread is ~4.3×). `inverted = signed-area sign flip vs snapshot
majority`. Hardcoding 4× would false-positive on a clean GIF.

### B. Pacing (from Probe 2 — uniform perceived motion)

Replace fixed-index snapshot selection + variable hold with **equal-motion quantile
resampling + uniform hold**. In `_stage_data()`, after `hist` is captured, compute per-
iteration mean node displacement, its cumulative sum `Dcum`, then select `N_SNAP=32`
snapshots at equal quantiles of `Dcum` (with a zero-motion uniform fallback and a min-count
floor that tops up via uniform indices if `unique()` collapses selections):

```python
Pn = [np.asarray(p)[:, :2] for (p, _t) in hist]
d_iter = np.zeros(len(Pn))
for i in range(1, len(Pn)):
    d_iter[i] = np.mean(np.linalg.norm(Pn[i] - Pn[i - 1], axis=1))
Dcum = np.cumsum(d_iter)
N_SNAP = 32
if Dcum[-1] <= 0:
    sel = np.unique(np.linspace(0, len(hist) - 1, N_SNAP).astype(int))
else:
    targets = np.linspace(0.0, Dcum[-1], N_SNAP)
    sel = np.clip(np.searchsorted(Dcum, targets), 0, len(hist) - 1)
    sel[0] = 0; sel[-1] = len(hist) - 1
    sel = np.unique(sel)
    if len(sel) < N_SNAP:
        sel = np.unique(np.concatenate(
            [sel, np.linspace(0, len(hist) - 1, N_SNAP).astype(int)]))
snap_idx = sel.tolist()
# F-03: do NOT append (p_out, t_out) — the solver RETURN keeps the pre-move stale t that
# snapshot_retriangulate never fixes. sel[-1] == len(hist)-1 already selects the final
# iterate (re-triangulated) as the converged-truss frame.
snapshots = [hist[i] for i in snap_idx]
```

Then in `main()`, replace the tiered-hold play-list build with a **uniform hold**
(`TRUSS_HOLD=2`): 3 preroll seed frames + `TRUSS_HOLD × len(snapshots)`. Probe 2 measured
SC-002 thirds ratio **1.261** (vs baseline 14.888 FAIL), 67 truss frames (3-frame preroll + 2×32 snapshots; probe measured 69 with the since-dropped appended final, F-03 corrigendum) = 6.7 s @ 10 fps
(inside the 5–8 s target), per-transition motion mean ≈ 10 px @ 864-wide. The algorithm is
quantile-based, so it self-adapts if Probe 1's dynamics shift.

**Converged-truss frame + FEM input (F-03 corrigendum).** `sel[-1] == len(hist) - 1` already
selects the last iterate as the final rendered truss frame, so **no** extra snapshot is
appended. Because the solver *return* `(p_out, t_out)` keeps the pre-move stale `t` (§A),
`_stage_data()` must derive the converged-truss mesh **and** the FEM input from `hist[-1]`,
not the return: `p_fin, t_fin = hist[-1]`; compute `q_truss = _quality(p_fin, t_fin)`; build
the smoother input `Mesh(connectivity=np.c_[t_fin, t_fin[:, 0]], points=np.c_[p_fin,
np.zeros(len(p_fin))], ...)` and compute `q_fem` on `t_fin`. `audit_hero_gif` asserts this
FEM-input connectivity (`hist[-1]`) passes the degenerate audit **by name** — it is already
covered by the 0/200 snapshot sweep, but is asserted explicitly so a FEM-path regression
cannot slip through.

**Seed triangulation (F-11).** The seed-stage wireframe/histogram must show the true *pre-
solver* scattered state, so `T_seed` is computed in-generator as a fresh `Delaunay(P_seed)`
filtered by the centroid-in-domain keep test — **not** `snapshots[0][1]`, which under
`snapshot_retriangulate=True` is the triangulation *after* iteration 0's move. `q_seed =
_quality(P_seed, T_seed)` then reflects the genuine unrelaxed seed.

### C. Peel + histogram + metrics (from Probe 3)

`n_layers = 4` for the annulus. Reveal order is **boundary(s)→inward** verified by
strictly-monotonic mean rim-distance (L0 0.0295 → L3 0.3225); do **not** order/verify by
raw radius (non-monotonic — L0 hugs both rims). State variable `k` = number of fully-
revealed layers (0..n_layers); `k=0` is byte-identical to the FEM-hold appearance
(satisfies convert-in-place FR-006). Function signatures + constants:

```python
HBINS = 40; RED = "#ff5555"          # F-08: RED is the ONE new constant. Reuse the existing
                                     # GOOD ("#2ecc71") for the green Median and TEXT
                                     # ("#e8e8ee") for the neutral Min — no GREEN/NEUTRAL dups.

def _layer_color(li, n_layers):                          # single-layer viridis color
    return matplotlib.colormaps["viridis"](Normalize(0, max(1, n_layers - 1))(li))

def _peel_facecolors(q, elem_layer, n_layers, k):        # mesh facecolors for state k
    fc = matplotlib.colormaps["cool_r"](Normalize(0.0, 1.0)(q))   # all quality-colored
    for li in range(k):
        fc[elem_layer == li] = _layer_color(li, n_layers)         # overwrite revealed
    return fc

def _draw_peel_hist(ax, q, elem_layer, n_layers, k, D):  # stacked bars, totals invariant
    # Per-element bin index (F-10): shared HBINS edges over [0,1]; a q=1.0 element lands in
    # the last bin, never overflows ->  b = np.clip(np.floor(q * HBINS).astype(int), 0, HBINS-1)
    # Bottom-up stack: revealed layers 0..k-1 each in viridis(li); the REMAINDER (elements
    # with elem_layer >= k) is colored PER BIN, pinned EXACTLY (F-04) to
    #   cool_r(Normalize(0, 1)(bin_midpoints))
    # i.e. byte-identical to _draw_hist's `bar_colors` (mids = hedges[:-1] + widths/2), so at
    # k=0 the whole histogram equals the FEM-hold histogram pixel-for-pixel (convert-in-place, FR-006).
    ...  # counts via np.add.at into HBINS; bottom-accumulate the stacked segments

def _draw_metrics(ax, q, alpha=1.0):                     # REPLACES ax.set_title (F-09: alpha)
    med, mn, mean = float(np.median(q)), float(np.min(q)), float(np.mean(q))
    ax.axvline(med,  color=GOOD, linestyle=":", linewidth=1.8, alpha=0.9 * alpha, zorder=20)
    ax.axvline(mean, color=RED,  linestyle=":", linewidth=1.8, alpha=0.9 * alpha, zorder=20)
    ax.text(0.02, 1.06, f"Median: {med:.3f}", transform=ax.transAxes, color=GOOD,
            fontsize=11, fontweight="bold", ha="left", va="bottom", alpha=alpha)
    ax.text(0.40, 1.06, f"Mean: {mean:.3f}",  transform=ax.transAxes, color=RED,
            fontsize=11, fontweight="bold", ha="left", va="bottom", alpha=alpha)
    ax.text(0.76, 1.06, f"Min: {mn:.3f}",     transform=ax.transAxes, color=TEXT,
            fontsize=11, ha="left", va="bottom", alpha=alpha)
```

`_draw_metrics` reaches the **seed/truss/FEM** stages by replacing the single `ax.set_title(...)`
in `_draw_hist` (those frames call `_draw_hist`). **The peel-reveal frames do NOT call
`_draw_hist` — they call `_draw_peel_hist`** (F-01), so the peel frame function must call
`_draw_metrics(ax, q_fem)` explicitly, immediately after `_draw_peel_hist`. Peel `q` is `q_fem`
(constant across `k` — peeling recolors, it does not re-smooth), so the Median/Mean/Min lines and
texts are identical on every peel frame. Probe 3 verified per-bin stacked totals are
`np.array_equal` across k=0..4 and equal to the plain histogram (FR-006 / SC-004 invariant), and
that median (0.930) / mean (0.870) lines are 2.4 bins apart — distinguishable, color+position keyed.

**Peel frame schedule (10 fps):** k=0 hold **18 f (1.8 s)** — doubles as the FR-004 FEM→Peel
inter-stage pause (F-13: all inter-stage holds standardized to 18 f) → k=1/k=2/k=3 hold 8 f
each → k=4 final hold 20 f (2 s, FR-005 held full-layer). Total peel = **62 frames**. Recolor
is **instant** at each reveal boundary, not a fade.

### D. Frame + size budget (from Probe 4)

**Held frames are free**: matplotlib's `PillowWriter` → PIL already merges runs of pixel-
identical consecutive frames into one stored frame with extended duration (proven byte-
identical: 30-unique vs 30-unique-held-4×). Budget by **unique** frames only. Cost per
unique frame: seed dots-only ≈ 10.8 KB, full mesh+hist ≈ 33.3 KB, peel ≈ 31.5 KB (raw
`PillowWriter`, no optimize).

**Total played frames after the F-13 hold changes ≈ 265 (26.5 s @ 10 fps):** seed 8+26+16 +
wire 10 + wire-hold **18** + truss (3 + 2×32 =) 67 + truss-hold **18** + FEM 24 + FEM-hold 16
+ peel **62**. **Unique** frames ≈ **110–120** (F-07: probe-4's conservative 50–80-truss-unique
budget was for a coarser snapshot set; the shipped 32-snapshot config lands here — the earlier
"~179 unique" was an overcount, no action needed). **Seed-frame budget correction (F-05):** only
`f_bfade`/`f_fall`/`f_settle` (~50 frames) are histogram-free ≈ 10.8 KB; `f_wire` (~10 unique)
**and** the wire-hold **intentionally carry the fading histogram** at ≈ 33 KB/frame (kept by
design), which adds ≈ **+220 KB** over a naive all-seed-dots-only assumption. Net raw projection
≈ **3.0–3.4 MB** — tight but under 4.2 MB even before optimize. Therefore a **mandatory
post-process re-encode** is appended after `anim.save(...)` in `main()`:

```python
from PIL import Image, ImageSequence
def _optimize_gif(path, colors=256):
    im = Image.open(path)
    frames = [f.convert("RGB").copy() for f in ImageSequence.Iterator(im)]
    durs   = [f.info.get("duration", 100)
              for f in ImageSequence.Iterator(Image.open(path))]
    stack = np.concatenate([np.asarray(frames[i])
                            for i in range(0, len(frames), max(1, len(frames)//8))], axis=0)
    pal = Image.fromarray(stack).convert("P", palette=Image.ADAPTIVE, colors=colors)
    q = [f.quantize(palette=pal, dither=Image.NONE) for f in frames]
    q[0].save(path, save_all=True, append_images=q[1:],
              duration=durs, loop=0, optimize=True, disposal=1)

anim.save(str(OUT_PATH), writer=animation.PillowWriter(fps=10), dpi=72)
_optimize_gif(str(OUT_PATH), colors=256)   # dims 864x396 + frame count + durations preserved
```

Probe 4 measured this at **−31.6%** on the real shipped GIF (3.77 → 2.58 MB) with identical
dimensions and frame count, moving the ≈ 3.0–3.4 MB raw projection to ≈ **2.1–2.3 MB** —
comfortably under 4.2 MB. `colors=192` (−36.6%) / `128` (−41.6%) are fallbacks if a hard
margin is needed; `dpi<72` shrinks below 864×396 and is **rejected** (FR-008 violation). Two
budget guardrails: keep **only the `f_bfade`/`f_fall`/`f_settle` seed frames** (~50) dots-only
with an empty histogram panel (so they stay ~10.8 KB, not ~33 KB) — the `f_wire` + wire-hold
histogram is intentional and **NOT** subject to this guardrail; if size still regresses, trim
truss holds before touching uniques.

### FR → design traceability

| FR | Design element |
|---|---|
| FR-001 (solver-side degeneracy suppression, additive/default-off) | A — `snapshot_retriangulate`/`snapshot_strict_interior`/`snapshot_hole_tau`; default path byte-identical (probe-proven). |
| FR-002 (real captured node positions, no interpolation) | A — only the render **connectivity** is re-triangulated; captured `p` positions are the untouched solver iterates. |
| FR-003 (uniform truss pacing, ≤3× band) | B — equal-motion quantile resampling + uniform hold (ratio 1.261). |
| FR-004 (≥1.5 s inter-stage pauses) | B/C schedule — all three inter-stage holds standardized to 18 f = 1.8 s (Seed→Truss wire-hold, Truss→FEM, FEM→Peel k=0), F-13; held frames are free (D). |
| FR-005 (peel boundary-inward reveal, ≥0.8 s/layer, held full-layer end) | C — reveal k=1..4 at 8 f (0.8 s) each, final hold 20 f, instant recolor. |
| FR-006 (lock-step stacked histogram, in-place, totals invariant) | C — `_draw_peel_hist` bottom-up stack; k=0 == quality-colored; per-bin totals `np.array_equal` across k. |
| FR-007 (green Median line+text, red Mean line+text, neutral Min) | C — `_draw_metrics` (green=GOOD, red=RED, neutral=TEXT; F-08) from `_draw_hist` for seed/truss/FEM, and called explicitly after `_draw_peel_hist` on every peel frame (F-01). |
| FR-008 (≤4.2 MB, 864×396, single loop) | D — post-process 256-color + optimize re-encode; dpi pinned at 72; loop ends on peel k=4 hold. |
| FR-009 (deterministic, fixed seed) | A–D — seed 11 unchanged; resampling/re-encode add no RNG; Delaunay + ADAPTIVE palette deterministic on fixed input. |

## Test Plan

Mirrors the existing `TestHistoryCapture` precedent in `tests/test_admesh_warmstart.py`.
Written **with** the solver change (Constitution III):

1. **Default-path byte-identity** (regression guard for FR-001). Run
   `distmesh2d_warmstart` twice on the same inputs — once with no history hook, once with
   `history_out=[]` and the new kwargs left at defaults — and assert
   `np.testing.assert_allclose(p1, p2)` / `t` equality. Extends the existing
   `test_history_none_is_default_noop`. Also assert a call passing
   `snapshot_retriangulate=False` explicitly returns output identical to the no-kwarg call.
2. **Opt-in behavior** (new `TestSnapshotRetriangulate`). With `snapshot_retriangulate=True,
   history_every=1` on a small annulus seed: for every captured `(p_i, t_i)`, assert `t_i`
   is a valid Delaunay of `p_i` (each triangle's centroid satisfies `fd(centroid) < 0`) and
   that **zero** triangles are inverted (signed-area sign matches the snapshot majority) —
   the exact defect Probe 1 found in the stale-`t` default capture. Assert boundary rows
   stay pinned (as the existing capture test does).
3. **Return-value invariance under the hook**: assert the `(p_out, t_out)` returned with
   `snapshot_retriangulate=True` equals the return with the hook off (physics unchanged) —
   Probe 1's `p_identical/t_identical` result as a regression gate.
4. **Generator smoke** (Constitution III — viz is smoke-only): a fast, low-`niter` or
   monkeypatched invocation asserting `main()` runs headless under `MPLBACKEND=Agg` and
   writes a non-empty GIF; full visual acceptance (SC-001, SC-006) is the manual
   frame-audit + size check in `quickstart.md`, not a unit assertion.

No change to any existing test is required (backward-compat gate). All prior
`test_admesh_warmstart.py` tests must still pass unmodified.

## Complexity Tracking

> No Constitution violations to justify. This table tracks **open risks** carried from the
> probes (all four probes succeeded; these are residual risks, not failures).

| Risk | Impact | Mitigation / trigger |
|---|---|---|
| Extra Delaunay per captured iteration (hero run only, ~200 calls) | Slower hero render | Measured within ~2 min budget; zero cost on default path. No action unless render exceeds budget. |
| Rendered connectivity ≠ physics connectivity | Fidelity question (FR-002) | FR-002 governs positions (untouched); captured `t` is a valid Delaunay of the real points and more honest than 1-step-stale `t`. Documented in solver + generator docstrings so a future edit does not revert it. |
| `snapshot_hole_tau` floor | Too-low tau culls inner-ring boundary tris → empty-ring frames | Keep tau ≥ ~0.01 (inner-ring chord midpoints dip ~0.004 into the hole); 0.02 verified safe at h0≈0.095. For a different h0, re-verify. |
| Probe-1 dynamics could shorten `len(hist)` | Fewer than 32 snapshots → shorter truss stage | Probe 2 floor + zero-motion fallback present; if `len(hist)` routinely < 40 after the solver change, lower `N_SNAP` or raise `TRUSS_HOLD` to keep the stage ≥ 50 frames (spec min-stage edge case). **Re-run Probe 2 after the solver change lands** to re-confirm the SC-002 ratio. |
| Raw GIF budget sensitive to truss unique count | > 4.2 MB at high truss-unique end (raw path) | Mandatory 256c + optimize post-process (−31.6%) gives margin; if omitted, hard-cap truss uniques ≤ 60. Flag for the render-and-measure gate (SC-006). |
| `f_bfade`/`f_fall`/`f_settle` accidentally carrying the histogram | +~1.1 MB (~50 frames × 22 KB) | Keep only those ~50 dots-only seed frames with an empty histogram panel. NOTE (F-05): `f_wire` + the wire-hold intentionally DO carry the fading histogram (~33 KB) — by design, not a regression. |
| Metric lines coincide in some non-peel stage (median≈mean) | Two dotted lines overlap | Distinct colors + zorder(20) keep both readable; low risk (peel data is 2.4 bins apart). |
| `elem_layer` ordering depends on `m.Layers`/`m.n_layers` | A future chilmesh layer-reordering could break boundary-inward assumption | Add a cheap rim-distance-monotonicity assert in the generator so a reorder is caught at render time. |
