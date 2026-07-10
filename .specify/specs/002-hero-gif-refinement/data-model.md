# Phase 1 Data Model: README Hero GIF Refinement

This feature has no persistent store or schema; its "entities" are the in-memory objects the
generator (`scripts/generate_hero_animation.py`) and the solver hook build per run. They are
captured here with attributes and invariants so the implementation and tests have a shared
vocabulary. All are derived deterministically from RNG seed 11 (FR-009).

---

## Solver snapshot

One captured `(points, triangulation)` pair from a real truss iteration; the atomic unit of
truss playback.

| Attribute | Type | Meaning |
|---|---|---|
| `p` | `(N, 2) float ndarray` | Node positions at this captured iteration (real solver iterate; N = 380 for the annulus, constant across iterations). |
| `t` | `(M, 3) int ndarray` | Triangle connectivity for this snapshot. |
| `iteration` | int | Solver iteration index it was captured at (0..199). |

**Invariants.**
- Positions are **real** solver iterates — never interpolated (FR-002).
- Boundary rows `p[:B]` are bit-exact-pinned equal to the input boundary in every snapshot.
- When captured with `snapshot_retriangulate=True`, `t` is a valid Delaunay of **this
  snapshot's own** `p`: every triangle centroid satisfies `fd(centroid) < -geps`, and no
  triangle is inverted (signed-area sign matches the snapshot majority). This is the FR-001
  fix for the one-Euler-step `(p_new, t_old)` staleness.
- When captured with the default (`snapshot_retriangulate=False`), behavior is byte-identical
  to the pre-feature hook (the `t.copy()` pre-move path) — the regression-guarded default.
- The captured `t` is a **render artifact**; the solver discards its own `t` next iteration.
  FR-002 governs positions, not this connectivity.

---

## Animation stage

One of Seed / Truss / FEM / Peel; owns a frame schedule and a boundary pause.

| Attribute | Type | Meaning |
|---|---|---|
| `name` | enum {Seed, Truss, FEM, Peel} | Stage identity (drives the on-canvas label). |
| `frames` | list of `(fn, i, n)` schedule entries | Ordered render callbacks. |
| `boundary_pause_frames` | int | Held settled frames at this stage's start/boundary. |

**Invariants.**
- Every stage boundary (Seed→Truss, Truss→FEM, FEM→Peel) holds ≥ 15 frames = 1.5 s @ 10 fps
  (FR-004 / SC-003). The peel k=0 hold (15 f) doubles as the FEM→Peel pause.
- Truss stage = 3 preroll + `TRUSS_HOLD × len(snapshots)` frames (uniform hold); its
  per-rendered-frame mean node displacement varies ≤ 3× across thirds (FR-003 / SC-002).
- Held frames are pixel-identical runs (merged free by PIL — do not budget them as bytes).
- The overall sequence is a single loop ending on the peel full-layer hold (no re-peel).

---

## Mesh layer

One concentric ring from `peel_layers()` / `m.Layers` (outer→inner ordering); owns a color
and an element set.

| Attribute | Type | Meaning |
|---|---|---|
| `index` (`li`) | int 0..n_layers−1 | Layer id; 0 = boundary-adjacent (outermost), n_layers−1 = mid-band core. |
| `elements` | int index set | Elements assigned to this layer (`elem_layer == li`). |
| `color` | RGBA | `viridis(Normalize(0, n_layers−1)(li))`. |
| `mean_rim_dist` | float | Mean element distance to nearest domain rim (the reveal-order key). |

**Invariants.**
- `n_layers = 4` for the annulus fixture; `elem_layer` partitions all 580 elements (every
  element belongs to exactly one layer).
- Reveal order is by **strictly increasing `mean_rim_dist`** (boundary→inward): L0 0.0295 <
  L1 0.1333 < L2 0.2552 < L3 0.3225. Raw radius is **not** monotonic and must not be used to
  order or verify the reveal.
- Layer assignment comes from `m.Layers` OE+IE exactly as `_stage_data()` computes it;
  a cheap rim-distance-monotonicity assert guards against a future chilmesh layer reorder.

---

## Quality histogram state

Per-frame element-quality distribution; in the peel stage it becomes a stacked composition
keyed by revealed layers.

| Attribute | Type | Meaning |
|---|---|---|
| `bins` | 40 fixed bins over [0, 1] | `HBINS = 40`; edges shared across all frames. |
| `revealed_k` | int 0..n_layers | Number of fully-revealed layers (peel-stage state variable). |
| `ymax` | int | Y-axis cap, pinned to the global max bin count across all stages. |
| `stack` | per-bin, per-layer counts + remainder | Bottom-up: layers 0..k−1 each in `viridis(li)`, then a remainder segment (layers ≥ k) in per-bin `cool_r`. |

**Invariants (FR-006 / SC-004).**
- `revealed_k = 0` renders byte-identically to the FEM-hold quality histogram (convert-in-
  place entry, spec clarification Q1→A).
- **Per-bin total height is invariant across all k**: `np.array_equal` holds for the stacked
  per-bin totals at k=0..n_layers, and each equals the plain `np.histogram(q_fem)` of the
  same data. Segments recolor; totals never change.
- A layer contributing zero elements to a bin yields an empty segment there; the bin total
  still matches.
- `ymax` is identical for every frame in the run (pinned once in `_stage_data()`).

---

## Metric annotation

A `(vertical line, header text)` pair, color-keyed per metric, drawn on **every** histogram-
bearing frame.

| Attribute | Type | Meaning |
|---|---|---|
| `metric` | enum {Median, Mean, Min} | Which statistic. |
| `line` | dotted axvline or none | Median → green dotted; Mean → red dotted; Min → none. |
| `text_color` | color | Median → green `#2ecc71`; Mean → red `#ff5555`; Min → neutral `#e8e8ee`. |
| `value` | float | The statistic's value for the frame's quality array. |

**Invariants (FR-007 / SC-005).**
- 100% of histogram-bearing frames (all stages, not just peel) carry the green Median line +
  green text and the red Mean line + red text; Min text stays neutral (no line).
- The two lines/texts remain individually distinguishable even if median ≈ mean (distinct
  colors + zorder 20); in the peel data they are 2.4 bins apart.
- Replaces the previous single `ax.set_title(...)` metric header — three separate colored
  `ax.text` objects at axes-fraction x = 0.02 / 0.40 / 0.76, y = 1.06.
