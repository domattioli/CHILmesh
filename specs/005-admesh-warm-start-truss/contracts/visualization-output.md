# Contract: Visualization Output (New 4-Row Layout)

**Date**: 2026-05-02
**Phase**: 1 — Design & Contracts
**Type**: Output artifact contract (PNG image structure)
**Supersedes**: `specs/004-fix-readme-viz/contracts/visualization-output.md`

---

## Purpose

The demo script produces a single PNG asset consumed by the project README. This contract is the PNG's structure and content. Spec 004's contract is replaced by this one; the file path is preserved so README links keep working.

---

## Output File Contract

| Property | Value |
|----------|-------|
| Path | `tests/output/annulus_quickstart.png` (UNCHANGED from spec 004) |
| Format | PNG, 8-bit RGBA, non-interlaced |
| Resolution | 100 DPI |
| Dimensions | ~1500 × 1800 px (figsize 15×18 in @ 100 dpi) |
| File size | < 5 MB |

### Image structure (4×3 subplot grid)

The image MUST contain exactly twelve subplots arranged in a 4×3 grid, with a figure-level title at the top.

#### Figure title

```text
CHILmesh × ADMESH: Warm-Start Truss Optimization Pipeline
```

#### Row contracts (NEW — replaces spec 004)

| Row | Mesh source | Smoother applied to that source |
|-----|-------------|--------------------------------|
| 1 | `chilmesh.examples.annulus()` | None (raw Delaunay over random points) |
| 2 | **Row 1's mesh** | `chilmesh.optimize_with_admesh_truss(row1, ANNULUS_SDF, size_fn=None, seed=0)` |
| 3 | **Row 2's mesh** (NOT row 1, NOT a fresh ADMESH) | `mesh.smooth_mesh(method='fem', acknowledge_change=True)` |
| 4 | **Row 2's mesh** (NOT row 3) | `admesh.quad_prep.smooth_for_quadrangulation(row2.points, row2.triangles, ANNULUS_SDF.fd, h=h0)` |

**Critical invariant**: Rows 3 and 4 BOTH branch off Row 2 — they are siblings, not sequential. Row 4 is NOT "row 3 with right-iso applied"; it is "row 2 with right-iso applied."

**Geometric invariant within each row**: All three columns within the same row MUST display the same node positions and the same triangle connectivity. Only the per-element coloring differs.

**Boundary-preservation invariant**: All four rows MUST have identical boundary point coordinates (subject to FR-013 V_BND and V_BND_PROP). Row 2's boundary equals Row 1's bit-exactly. Row 3 and Row 4 boundaries equal Row 2's within each smoother's documented tolerance (FEM smoother: bit-exact for fixed boundary; right-isoceles: ADMESH's pfix tolerance).

#### Column contracts (preserved from spec 004)

| Col | Visualization | Colormap | Element coloring rule |
|-----|---------------|----------|----------------------|
| 1 | Mesh wireframe | None (black on white) | `triplot` only — NO fill |
| 2 | Layer skeletonization | MATLAB parula (64 stops, ListedColormap) | `parula_cmap(layer_idx / max(1, n_layers - 1))` |
| 3 | Element quality | matplotlib `cool_r` (red→blue, NOT plain `cool`) | `cool_r_cmap(quality_value)` normalized to [0,1] |

The `cool_r` (reversed) choice is from spec 004's clarify Q5=b: red = poor quality (0), blue = good quality (1). Inherited unchanged.

#### Subplot title contracts

Each subplot has a two-line title `"<Row Label>\n<Column Label>"`:

| Row | Row Label |
|-----|-----------|
| 1 | `Raw Delaunay` |
| 2 | `+ ADMESH Truss (warm-start)` |
| 3 | `Row 2 + FEM Smoother` |
| 4 | `Row 2 + Right-Isoceles Smoother` |

| Column | Column Label |
|--------|-------------|
| 1 | `Mesh` |
| 2 | `Layers` |
| 3 | `Quality` |

Example: top-left subplot title = `"Raw Delaunay\nMesh"`. Row 3 first column = `"Row 2 + FEM Smoother\nMesh"`.

#### Colorbar contracts (preserved from spec 004)

| Subplot | Colorbar |
|---------|----------|
| Col 1 (any row) | NONE |
| Col 2 (any row) | Discrete parula colorbar with one tick per layer (n_layers ticks). Label: `"Layer"` |
| Col 3 (any row) | Continuous cool_r colorbar in [0, 1]. Label: `"Quality"` |

---

## Programmatic verification (run before saving PNG)

The demo script enforces these assertions in order. Any failure raises `RuntimeError` and the PNG is NOT written.

| Check ID | What it verifies | Failure action |
|----------|------------------|---------------|
| **V_BND** | `np.array_equal(row2.points[boundary_indices], row1.points[boundary_indices])` — Row 2's boundary points are bit-exact equal to Row 1's | Raise `RuntimeError("V_BND failed: warm-start did not preserve boundary; max delta = ...")` |
| **V_BND_PROP** | Row 3 and Row 4 boundary points are equal to Row 2's, within each smoother's documented tolerance | Raise `RuntimeError("V_BND_PROP failed in Row {i}: max delta = ...")` |
| **V_QI** | `np.median(row2.elem_quality()[0]) > np.median(row1.elem_quality()[0])` — warm-start improved demo mesh quality | Raise `RuntimeError("V_QI failed: warm-start did not improve quality (in={...}, out={...})")` |
| **V_CONN** | All four rows: every triangle has positive area, no duplicate vertex indices | Raise `RuntimeError("V_CONN failed in Row {i}: {detail}")` |
| **V_CHAIN** | Row 3's input mesh is `id(row2)` (or hash-equal); Row 4's input mesh is `id(row2)` (or hash-equal); neither used Row 1 directly nor a fresh ADMESH | Raise `RuntimeError("V_CHAIN failed: Row {i} input was not Row 2")` |
| **V_TRUSS_INVOKED** | The vendored `distmesh2d_warmstart` was called exactly once for Row 2 (tracked via module-level flag) | Raise `RuntimeError("V_TRUSS_INVOKED failed: warm-start truss was not actually invoked for Row 2")` |

If any check fails, the PNG is NOT written — the README's existing image (the last successfully validated one) remains in place.

---

## Visual acceptance (human inspection)

A reviewer opening the rendered PNG should verify:

- [ ] Each row clearly demonstrates a different stage of the pipeline.
- [ ] Row 1 looks visibly worse (irregular triangles) than Row 2.
- [ ] Row 2's outer/inner boundaries are visually identical to Row 1's (the rings).
- [ ] Rows 3 and 4 both show further refinement of Row 2's mesh, but in different ways (FEM = continuous quality improvement; right-iso = elements pushed toward right-isosceles shapes).
- [ ] Left column has only black mesh edges on white — no color fills anywhere.
- [ ] Center column gradient is parula (deep blue at low end, yellow at high end).
- [ ] Right column gradient is cool_r (red at low quality, blue at high quality).
- [ ] Each subplot's title makes the row/column purpose immediately obvious.
- [ ] Row 2's title explicitly says "ADMESH Truss (warm-start)" — not just "ADMESH" or "Smoothed".
- [ ] Within each row, the geometry in all three columns is identical.
- [ ] Image renders correctly at GitHub's default README viewport width.

---

## Migration from spec 004

What changed:
- **Row 2** was "FEM smoother applied to Row 1" → now "ADMESH warm-start applied to Row 1".
- **Row 3** was "Fresh ADMESH from bbox + size function" → now "FEM smoother applied to Row 2".
- **Row 4** was "Right-isoceles applied to Row 3" → now "Right-isoceles applied to Row 2" (siblings with Row 3, not sequential).
- **Title** updated to reflect the new pipeline narrative.
- **V_TRUSS_INVOKED** added (was V3 in spec 004 for right-isoceles; now applies to warm-start).
- **V_BND_PROP** is new (spec 004 didn't need it because rows 3-4 came from a different mesh).
- **V_CHAIN** is new (specifically guards against the regression "rows 3-4 secretly used row 1 because someone restructured the script wrong").

What stayed the same:
- File path (`tests/output/annulus_quickstart.png`).
- 4×3 grid layout.
- Column conventions (mesh / layers / quality).
- Colormaps (parula 64-stop ListedColormap / cool_r continuous).
- DPI and figsize ratio.

---

## Non-contracts (out of scope)

- Element count per row is NOT fixed (varies with truss convergence).
- Specific RGB pixel values are NOT fixed (depend on matplotlib version anti-aliasing).
- Run-to-run determinism is enforced **only when** all dependency versions and the RNG seed are pinned (per the API contract's Determinism section). On a different machine with different scipy version, pixels may differ but structure is stable.
- Animation, interactivity, or alternate output formats (SVG, PDF) are NOT in scope.
- The `generate_4row_admesh.py` filename is NOT contractually fixed; rename is permitted as long as README is updated and the PNG path stays at `tests/output/annulus_quickstart.png`.
