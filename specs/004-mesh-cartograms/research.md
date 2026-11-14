<!-- provenance: author=domattioli model=unknown effort=unknown date=2026-09-24 skill=speckit-pipeline repo=CHILmesh session=unknown -->
# research: size-controlled mesh cartograms (#219)

Advance-only research note for [#219](https://github.com/domattioli/CHILmesh/issues/219)
(`request: research` / `status: brainstorming` / `priority: someday`). Grounds the
proposed `mesh.cartogram(...)` API against the **current** CHILmesh plotting + layer
API so a later implementer starts from real call sites, not the issue's prose. No code
shipped; no label flip (advance only, per the brainstorming entitlement).

---

## 1. Problem (restated, one line)

An area-weighted `tripcolor` of a non-uniform mesh spends most of the canvas on the
few coarse elements and shrinks the many refined elements — where the scalar signal
concentrates — to sub-pixel threads. A cartogram equalizes *visual* weight per element
so a per-element (or per-vertex) scalar field is readable regardless of geometric size.

## 2. Existing substrate to compose (verified in-tree)

The cartogram is **additive** — it re-lays-out element geometry, then hands the new
polygons to the plotting path that already exists. No new render engine.

| Need | Existing API (verified) | File |
|---|---|---|
| per-element polygons from `(points, connectivity)` | `build_polygons(points, connectivity, ...)` | `src/chilmesh/chilplotting.py:100` |
| scalar → colormap fill of a polygon set | `plot_filled(points, connectivity, *, values, cmap, vmin, vmax, ...)` | `chilplotting.py:205` |
| length-check `len(values) == n_elems` | already enforced in `plot_filled` | `chilplotting.py:227` |
| per-layer boundary vertex **path ordering** | `paths_on_outer_vertices(mesh, layer_idx) -> list[np.ndarray]` | `src/chilmesh/layer_paths.py:51` |
| layer decomposition (OE/IE/OV/IV per layer) | `mesh.layers` dict of lists | `CHILmesh.py:307`, docstring `:53-56` |
| axis config / limits | `configure_axes`, `_new_ax` | `chilplotting.py:72,162` |

Consequence: `noncontig` and `hybrid` need only a **centroid + rescale** step, then
`PolyCollection` fill; `unrolled` is the only variant that consumes `mesh.layers` +
`paths_on_outer_vertices`, both of which already exist.

## 3. Proposed API

```python
fig, ax = mesh.cartogram(
    values,                 # (n_elems,) or (n_verts,) scalar; per-vert reduced to per-elem by mean
    kind="unrolled",        # "noncontig" | "unrolled" | "hybrid"
    reference="mean",       # hybrid only: mean|median|mode|max|min  → element extent ∝ |value-ref|+eps
    cmap="inferno",
    ax=None,
) -> tuple[Figure, Axes]
```

Thin instance method on `CHILmesh` delegating to a new `chilplotting.cartogram(mesh, values, kind, ...)`
free function (keeps the mixin thin, mirrors the `chilplotting.plot`/`axis_chilmesh(mesh)` pattern at
`chilplotting.py:299,304`). Returns Matplotlib handles like every other plot fn here.

## 4. Algorithm per `kind` (grounded)

- **`noncontig` — equal-area, geographic position preserved**

    I. compute each element centroid `c_e` from `points[connectivity[e]]` (mean of its verts).

    II. redraw every element as a fixed-area glyph about `c_e` (square or the element's own shape
        scaled to a constant target area `A = bbox_area / n_elems`). Position = true centroid → the
        map still reads geographically; only per-element *area* is equalized.

    III. hand the rescaled polygons + `values` to the `PolyCollection` fill (reuse `plot_filled`'s
         norm/cmap branch, `chilplotting.py:225-233`).

    Known cost: equal-area glyphs at true centroids **overlap** in refined bands and **gap** in coarse
    ones — acceptable for a first cut; the Dorling equal-circle / Gastner–Newman diffusion variants in
    the backlog fix overlap but are materially more code (deferred, §7).

- **`unrolled` — skeleton layers as equal-height bands**

    I. peel layers already available in `mesh.layers` (`_peel` output). Each layer → one horizontal band.

    II. within a band, order elements by the **OV path arc-length** from `paths_on_outer_vertices(mesh, i)`
        (true boundary order) — falls back to centroid-angle if a layer has no clean path (islands / multi-loop).

    III. each element → an equal-width cell in its band; colour by `values`. Band height = equal-per-layer
         (v1); count-weighted / log-count are backlog knobs (§7).

    This is the variant that exposes the #219 motivating result (deep small layers L21–L25 clustering) —
    it deliberately discards geography for a layer×angle readout.

- **`hybrid` — focus+context**

    element extent ∝ `|value - reference| + eps`, `reference ∈ {mean, median, mode, max, min}` of `values`.
    Same centroid-anchored redraw as `noncontig` but the per-element target area is signal-driven, so
    outliers vs the reference grow and on-reference elements shrink. Reuses the `noncontig` layout code
    with a variable area vector.

## 5. Dependencies / gating

- **matplotlib only** for v1 — no new extra. `PolyCollection` is already the render path
  (`plot_filled`), so cartograms inherit the same headless/`Agg` behavior and the same
  ~11.6 µs/element raster ceiling measured in #167. Large-mesh (10⁶+) batching/decimation is a
  backlog item, NOT a v1 blocker (the QuADMESH study meshes are ≤10⁵ elems).
- No GPU dependency — this is orthogonal to the #167 GPU-backend track (that issue is env-blocked;
  this one is not, it composes the existing matplotlib path).

## 6. Validation / test plan

- **API/shape tests** (headless `MPLBACKEND=Agg`, mirrors existing plotting tests): each `kind`
  returns `(Figure, Axes)`; `values` length mismatch raises (reuse the `plot_filled` guard); per-vertex
  input is reduced to per-element without error.
- **Invariant tests**: `noncontig`/`hybrid` preserve element **count** (n polygons out == n elems in);
  `unrolled` places every element in exactly one band and `sum(cells per band) == n_elems`.
- **Reference reproduction** (acceptance, #219): regenerate the QuADMESH pass-frequency cartograms
  (PR #96 / #97 `experiments/mc_layer_pass/`) directly from a CHILmesh mesh + a per-element scalar,
  asserting the fine-coastal-band structure is visible (non-background pixel fraction for the refined
  9.9%/4.9% element sets exceeds the plain-map baseline). Gallery doc page (MATLAB-help register docstrings).

## 7. Backlog (from the issue, tiered so v1 stays small)

Deferred past v1 (each is its own follow-up, not a v1 blocker): Dorling equal-circle + Gastner–Newman
contiguous diffusion (overlap/gap fix); count-weighted / log-count band heights; true OV arc-length vs
centroid-angle within-band ordering for multi-loop layers + islands; diverging colormaps + percentile-clip
+ colorblind-safe defaults; `PolyCollection` batching/decimation for 10⁶ elems; optional hover interactivity.

## 8. Open questions for the operator (toward `ready`)

1. **v1 scope** — ship all three `kind`s, or land `noncontig` + `unrolled` first and defer `hybrid`?
   (`hybrid` reuses `noncontig` layout, so marginal cost is low — leaning all-three.)
2. **`noncontig` glyph** — fixed square vs the element's own shape scaled to constant area? (square is
   simpler + reads as a cartogram; shape-preserving is prettier but adds per-element affine work.)
3. **`unrolled` band height** default — equal-per-layer (simplest, recommended v1) vs count-weighted?
4. **Priority** — #219 is `someday`; confirm it stays deferred, or promote if the QuADMESH pass-frequency
   visualization is now a near-term need.

A one-line answer to Q1–Q3 turns this into a `status: ready` spec (plan → tasks → implement). Until then
`mesh.cartogram()` is intentionally **not** added.

## 9. Operator rulings (grill, 2026-10-20) — supersede §3, §4 and §8 where they differ

Line refs re-verified on `fix/issues-293-290` (a2ac607): `build_polygons` is now `chilplotting.py:193`,
`plot_filled` `:298`, `plot_field` `:626` (per-vertex → per-element mean + NaN blanking; reuse it),
`plot_layer` `:535` (OE+IE per layer). Prototypes: QuADMESH `origin/monte-carlo-quality-analysis`,
`experiments/mc_layer_pass/cartogram.py`, `cartogram_hybrid.py`.

| # | Question | Ruling |
|---|---|---|
| 1 | v1 scope | All three views ship together. |
| 2 | Hybrid exposure | Not a third kind. `kind ∈ {"noncontig","unrolled"}`; `reference=None` (default) = equal size, `reference ∈ {"mean","median","mode","max","min"}` or a float = focus+context. Applies to both layouts (prototype renders `hybrid_geographic` and `hybrid_unrolled`). |
| 3 | noncontig glyph | Element's own shape, scaled about its centroid by `sqrt(A / area_e)`. |
| 4 | noncontig target area `A` | Median element area (prototype `cartogram.py:49`). |
| 5 | noncontig draw order | Ascending by value; high values drawn last, on top (prototype `:54`). |
| 6 | unrolled band contents | One band per layer = `OE[i] ∪ IE[i]` (partitions all elements; annulus 580/580 unique). Layer 0 at top (prototype `:103`). |
| 7 | unrolled band height | Equal per layer. |
| 8 | unrolled within-band order | Centroid angle `atan2` about mean of element centroids (prototype `:89`). OV arc-length → backlog (annulus has 2–3 OV paths per layer, so multi-loop is the norm). |
| 9 | per-elem vs per-vertex input | Infer from length; per-vertex reduced like `plot_field`. If `n_elems == n_verts`, `on="elem"|"vertex"` is required, else `ValueError` naming both. |
| 10 | `reference="mode"` | Exact-value mode (`np.unique` counts; ties → smallest; prototype rounds to 6 dp). Docstring: intended for discrete fields. |
| 11 | `eps` | Default `0.004 × (max − min)` of finite values; `eps=` overrides with an absolute value. |
| 12 | Geographic focus+context clamp | Per-element scale clipped to `scale_clip=(0.05, 10)` × base (prototype `cartogram_hybrid.py:114`). In unrolled, cell width ∝ `|v − ref| + eps` within the band. |
| 13 | NaN | Element blank (unrolled: slot kept, empty); excluded from reference, colour limits and eps. |
| 14 | Acceptance | CI tests on `examples.annulus` with synthetic fields + a gallery script rendering `WNAT_Hagen.14` with a CHILmesh-computed scalar (quality). Exact pass_freq reproduction stays QuADMESH-side (`results/*_pass.csv` is gitignored there). |
| 15 | Glossary | Terms (cartogram, layout, band, reference) defined in the docstring + gallery page; no `CONTEXT.md`. |
| 16 | Priority | Promoted to now. Branch `feat/219-cartograms` off `development` (`git.md:3`). |

Defaults carried from prototype without a question: `cmap="inferno"`, `edgecolor="none"`, returns `(Figure, Axes)`.

---

_[model: claude-opus-4-8, repo: CHILmesh, session: 013NhKuhHJwC9wmmysELCufc, refs: #219, #167, QuADMESH#96/#97]_
