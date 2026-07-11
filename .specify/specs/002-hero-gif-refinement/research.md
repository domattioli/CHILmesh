# Phase 0 Research: README Hero GIF Refinement

Four live probes replicated `_stage_data()` exactly (RNG seed 11, annulus fixture from
`chilmesh.examples.annulus()`, `h0` = median fixture edge length ≈ 0.09531, rejection-
sampled interior in band `R_IN+0.04 .. R_OUT-0.04`, solver params
`deltat=0.02 Fscale=1.2 niter=200 track_best_quality=False history_every=1`,
`fd(p)=max(r-1.0, 0.3-r)`, `bbox=(-1,-1,1,1)`). Each section below records the decision,
rationale, alternatives rejected, and the verbatim measured evidence.

---

## Probe 1 — Degeneracy suppression (FR-001, FR-002)

**Decision.** Fix the degenerate/overlapping truss-frame elements at the **solver source**
by adding three additive, default-off keyword-only kwargs to `distmesh2d_warmstart`:
`snapshot_retriangulate=False`, `snapshot_strict_interior=False`, `snapshot_hole_tau=0.02`.
The hero run enables `snapshot_retriangulate=True` (and, harmlessly, `snapshot_strict_
interior=True`). The generator's degenerate-audit metric definitions are corrected in the
same pass.

**Rationale.** Two root-cause hypotheses were tested; only a subtle one is real. The visible
red "overlapping-sail" glitch is caused by the history hook appending
`(p.copy(), t.copy())` **after** `p = p_new` while `t` is still the pre-move triangulation
(computed from pre-move points). Every captured frame therefore pairs post-move points with
one-Euler-step-stale connectivity; points that crossed an edge within one `deltat` step
produce inverted overlapping triangles. Re-triangulating each captured snapshot's own post-
move points before appending fixes it — purely capture-side, physics/forces/return byte-
untouched (the solver discards its own `t` next iteration anyway).

**Alternatives considered / rejected.**
- *Display-side culling / snapshot-skipping* (spec clarification Q1 options A/B): **rejected
  by the operator** in favor of option C (solver-side prevention) — degenerates must not
  exist in any captured snapshot, not be hidden at render time.
- *Hypothesis (a): hole-spanning triangles surviving the centroid keep-test* — **FALSE.**
  Genuine hole-span (any edge-midpoint `fd > 0.02`) = 0 in every frame; the solver's
  centroid keep-test already culls true hole-spanners. The "60/61 flagged at `fd>-geps`"
  are legitimate inner-ring boundary tris whose chord midpoints dip only ~0.004 into the
  hole.
- *Hypothesis (b): stale triangulation from `ttol` skipping re-Delaunay* — **FALSE as
  stated** (`t` changes essentially every iteration; retriangulation is effectively 100%).
  The real defect is the subtler post-`p=p_new`/pre-move-`t` pairing above.
- *`snapshot_strict_interior` as the primary lever* — measured **no-op** on the annulus (no
  hole-spanners: 54→54 oversize at 4×, 0→0 at 8×). Kept ON only as belt-and-suspenders for
  non-annulus domains; the load-bearing lever is `snapshot_retriangulate` alone (fresh
  Delaunay already zeros inversions).
- *Prior dead ends (context):* default `deltat=0.2` diverges on a cold scattered seed (hence
  the `deltat=0.02` converging regime); interpolated/"lerp" fake frames were rejected by the
  operator (#179 real-frames mandate).
- *Generator metric bars:* `hole-span = fd > -geps` and `oversize = 4× median` were rejected
  as false-positive-prone — corrected to `fd > +0.02` and `8× median` (a valid annulus
  mesh's natural max/median area spread is ~4.3×, so 4× flags legitimate sparse-region tris).

**Measured evidence (verbatim).**
- BEFORE (baseline, stale-`t`), over 39 animation-sampled snapshots + all 200: inverted
  tris 30 total on sampled (20 @ iter0, 8 @ iter1, 2 @ iter3, sporadic), MAX 20/frame;
  oversize(4×) 54 sampled / 312 over all 200; oversize(8×) 0 with fresh-`t` but present with
  stale; genuine hole-span(>0.02) = 0; maxRel area 14.8 @ iter0. Quality traj (median
  script-q): start 0.6977, iter15 0.8841, final 0.8764.
- AFTER (winning: `snapshot_retriangulate=True, snapshot_strict_interior=True, tau=0.02`),
  audited over ALL 200 snapshots with corrected metrics (`hole_tau=0.02, over_bar=8×`):
  **TOTAL hole=0, oversize=0, inverted=0; MAX-per-frame all 0; frames with ANY degenerate
  0/200.** Quality traj: start 0.7580, iter15 0.8883, final 0.8847 (min 0.758, max 0.900,
  last-10 flat 0.874–0.894 — no oscillation, final ≥ 0.85). Slightly higher than baseline
  because inverted q=0 tris no longer drag the median.
- CONVERGENCE-IDENTITY: returned `p_identical=True, t_identical=True, maxabs dp = 0.00e+00`
  vs the original off-hook solver → physics provably unchanged.
- Fresh-Delaunay alone zeroed inversions (nInv 30→0); adding `strict_interior` changed
  nothing (confirming no hole-spanners exist).
- Spot-checks `spot_iter{0,1,3}.png`: baseline red inverted sails at the hole → gone in the
  winning config; visually clean watertight rings, no empty regions.

---

## Probe 2 — Uniform convergence pacing (FR-003, SC-002)

**Decision.** Replace fixed-index snapshot selection + variable hold with **equal-motion
quantile resampling + uniform hold**: select `N_SNAP=32` snapshots at equal quantiles of
cumulative node displacement, hold each `TRUSS_HOLD=2` frames. Truss stage = 3 preroll +
2×33 = 69 frames = 6.9 s @ 10 fps.

**Rationale.** Per-iteration mean displacement is near-uniform (~0.010, peak 0.0166 @ iter1)
and total motion is spread evenly (50% by iter99, 90% by iter179) — motion is **not** front-
loaded, so quantile resampling maps cleanly. The current `snap_idx` (early gaps=1 iter/snap
held 3 frames; late gaps=15 iters/snap held 1 frame) produces the exact "accelerates toward
end" defect. Quantile-based selection is dynamics-agnostic, so it self-adapts if Probe 1
changes solver dynamics; a min-count floor + zero-motion fallback keep it robust.

**Alternatives considered / rejected.**
- *Keep fixed-index selection with tuned tiers* — rejected: it is exactly what produced the
  14.9× ratio; any hand-tuned tier set is brittle to the solver dynamics.
- *Tuning sweep configs* (all pass ≤3, listed for the record): N24/H2=53f/ratio1.344;
  N24/H3=78f/1.286; N28/H2=61f/1.193; N32/H3=102f/1.098 (too long); N40/H2=85f/1.124
  (slightly over the 80-frame ceiling); N40/H3=126f/1.119 (too long); N32/H3=102f (too
  long). **N32/H2** chosen as the best duration/ratio balance inside 50–80 frames / 5–8 s.
- *Interpolated intermediate frames to smooth motion* — rejected by the #179 real-frames
  mandate (FR-002).

**Measured evidence (verbatim).**
- BASELINE (current): 39 snapshots, 79 frames, 7.9 s; thirds mean-per-frame displ ×1e3 =
  [4.24, 4.92, 63.11]; SC-002 ratio = **14.888 (FAIL)**.
- RECOMMENDED (N=32, HOLD=2): 32 unique snaps (+1 appended = 33 shown), 69 frames, 6.9 s;
  thirds ×1e3 = [26.25, 33.10, 27.72]; SC-002 ratio = **1.261 (PASS)** (1.098 preroll-
  excluded). Chosen hist iters = [0,5,12,18,25,32,38,45,51,57,63,70,76,82,89,95,102,109,
  116,122,128,134,141,147,154,161,167,173,180,187,193,199]. Per-transition displ (u): mean
  0.0626 / min 0.0250 / max 0.0713 → px@864: mean 10.0 / min 4.0 / max 11.4. Seed→hist[0]
  first-motion jump = 0.0250 u = 4.0 px.
- Pixel scale: fig 12×5.5 in @ 72 dpi = 864×396; mesh axis frac 0.4162 → 359.6 px; data
  range xlim(−1.12,1.12) = 2.24 u → 160.5 px/unit.
- Solver: 200 iters captured (no early dptol break), node count constant 380, `d_i` max
  0.0166 @ iter1, `D_final` = 1.978, motion evenly distributed.

---

## Probe 3 — Peel reveal + stacked histogram + metrics (FR-005, FR-006, FR-007)

**Decision.** `n_layers=4`; reveal **boundary(s)→inward** with state `k` = number of fully-
revealed layers; mesh facecolors = `cool_r(q)` overwritten by `viridis(li)` for elems with
`elem_layer < k`; histogram = bottom-up stacked bars (revealed layers in viridis, remainder
in per-bin `cool_r`) with per-bin totals invariant; every histogram-bearing frame draws a
green dotted Median line + green text, red dotted Mean line + red text, neutral Min text
(replacing the current single `ax.set_title`). Peel schedule: k=0 15f (1.5 s, doubles as
FR-004 inter-stage pause) → k=1/2/3 8f each → k=4 20f; instant recolor per reveal.

**Rationale.** Reveal direction must be verified by **rim-distance**, not raw radius: raw
mean radius is non-monotonic because layer 0 hugs BOTH the outer (R_OUT=1.0) and inner
(R_IN=0.3) rims, so its mean radius (0.764) sits mid-range; mean rim-distance IS strictly
monotonic (L0 0.0295 → L3 0.3225), giving a well-defined boundary→mid-band-core order.
Convert-in-place (spec clarification Q1→A) means k=0 is byte-identical to the FEM-hold
quality coloring, so the stage enters seamlessly and only converts as layers reveal.

**Alternatives considered / rejected.**
- *Order/verify layers by raw radius* — rejected: non-monotonic on the annulus (L0 0.764 >
  L2 0.657) would mis-order the reveal.
- *Fade recolor between reveals* — rejected in favor of instant recolor at each reveal
  boundary (clearer layer↔quality mapping; also cheaper).
- *Histogram restart / redraw fresh per k (not in-place)* — rejected by spec clarification
  Q1→A (convert-in-place); un-revealed portions must keep the quality appearance until their
  layer reveals.
- *Min or Median taking the red line* — rejected by spec clarification Q1→A: **Mean** takes
  red, Median takes green, Min stays neutral.
- *Offsetting the viridis Normalize domain to brighten layer 0* — deferred: spec says
  `Normalize(0, n_layers-1)`; wireframe edges keep dark-purple layer-0 elems legible. Only
  revisit if the operator flags contrast.

**Measured evidence (verbatim).**
- `n_layers=4`; `n_elems=580`; ymax(pinned, q_fem)=140. Per-layer (layer, nelems, mean_rad,
  min_rad, max_rad, mean_rim_dist): L0 257, 0.7641, 0.3132, 0.9900, **0.0295** (boundary-
  adjacent, both rims); L1 153, 0.7130, 0.3878, 0.9426, **0.1333**; L2 136, 0.6566, 0.4993,
  0.8797, **0.2552**; L3 34, 0.6511, 0.6014, 0.7333, **0.3225** (mid-band core, deepest).
  mean_rim_dist strictly monotonic ↑ → reveal direction verified boundary→inward; raw
  mean_rad non-monotonic → do NOT use radius.
- Metrics on q_fem: median=0.9301, mean=0.8700, min=0.361; |median−mean| = 0.0601 = 2.40
  bins (HBINS=40) → distinguishable, not coincident.
- Invariance: per-bin stacked totals `np.array_equal` across k=0..4 = **True**; == plain
  histogram = **True**; rendered-frame totals equal across k=[0,1,2,4] = **True**.
- Figure: figsize (12,5.5) dpi 72 = 864×396 px (matches FR-008).
- 4 PNGs rendered + inspected: stacked viridis colors distinguishable, metric texts colored
  + non-overlapping at 864×396, revealed layers viridis while unrevealed stay quality-colored
  (cool_r), legible over the dark bg via wireframe edges.

---

## Probe 4 — Frame + size budget (FR-008, SC-006)

**Decision.** GO on the redesign schedule, budgeting by **unique frames only** (held frames
are free), and ship with a **mandatory post-process re-encode** (single global 256-color
adaptive palette, `dither=NONE`, `optimize=True`, per-frame durations preserved) appended
after `anim.save(...)`.

**Rationale.** matplotlib `PillowWriter` → PIL already merges runs of pixel-identical
consecutive frames into one stored frame with extended duration, so every inter-stage hold,
peel per-layer hold, and final hold costs ~0 bytes — only unique frames cost bytes. The raw
path is budget-tight/marginally-over at the high truss-unique end (~4.4 MB at 80 truss
uniques); the re-encode is a large, visually-free win (mostly lossless LZW `optimize` +
one global palette instead of matplotlib's per-frame default palette with no optimize).

**Alternatives considered / rejected.**
- *Lower dpi to shrink bytes* — **rejected**: dpi 68 shrinks output to 816×374, violating
  FR-008's pinned 864×396.
- *Lower fps* — rejected as a size lever: fps changes only duration metadata, not bytes (and
  held-free makes fps irrelevant to budget anyway).
- *128- or 192-color re-encode as default* — deferred to fallback only: 128c can band the
  cool_r/viridis gradients; 256c has no visible loss and is the safe default. Drop to
  192c (−36.6%) / 128c (−41.6%) only if a hard margin is required.
- *Rely on the raw path with no post-process* — rejected as too tight; only acceptable if
  truss uniques are hard-capped ≤ 60.

**Measured evidence (verbatim).**
- Calibration: shipped GIF 3,771,075 B, 113 stored frames, 864×396, durations
  [100,200,300,900,1600,4000] ms, 24.3 s → 33.37 KB/unique; mimic full-frame 33.29 KB/unique
  (match).
- Held-free proof: 30-unique = 1,028,348 B; 30-unique-held-4× (120 nominal) = 1,028,348 B
  (identical). PIL path: A30 = B120 = 1,687,945 B.
- Per-unique cost (raw PillowWriter, no optimize): seed dots-only 10,758 B; full mesh+hist
  33,290 B; peel layer-colored 31,492 B; annotation (median+mean axvline + title) overhead
  1,591 B/unique.
- Redesign UNIQUE-frame projection (holds excluded): seed 70×10.76KB=753KB; truss 80×33.29KB
  =2,663KB (upper) | 50×=1,665KB (lower); fem 24×33.29KB=799KB; peel 5×31.49KB=158KB. RAW
  TOTAL = **4.37 MB (truss=80) | 3.37 MB (truss=50)**; break-even ~85 truss uniques.
- With post-process 256c+optimize (ratio 0.684): 4.37 MB → **2.99 MB** (upper) | 3.37 →
  **2.31 MB** (lower). Both < 4.2 MB. GO.
- Mitigation table (real GIF, dims+frames preserved): 256c+opt 3.77→2.58 MB (−31.6%);
  192c+opt →2.39 MB (−36.6%); 128c+opt →2.20 MB (−41.6%); dpi68 → 816×374 (FR-008 violation,
  rejected); fps → 0 byte effect.
- Guardrails: keep the 70 seed frames dots-only with an empty histogram panel (else they
  jump ~10.8 KB → ~33 KB, +1.6 MB); the post-process is deterministic (ADAPTIVE palette
  deterministic on identical pixel input; no RNG) — FR-009/SC-007 safe.
