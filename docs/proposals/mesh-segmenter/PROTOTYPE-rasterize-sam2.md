# Prototype Plan — rasterize → SAM2 field bootstrap

Goal: decide if SAM2 on a rasterized mesh field is worth building, with **zero model
training**. The prototype answers one question: **does SAM2 beat the cheap
deterministic baseline on element-IoU?** If no → drop SAM2, ship deterministic
`by_click(field_gradient)`.

This is throwaway/sandbox code. It does NOT belong in `mesh_segmenter`'s shipped API
or in chilmesh — keep it under `prototype/` with its own optional `[proto]` extra.

## Hypothesis

Rasterize a mesh scalar field → image → SAM2 point-click → mask → project back onto
elements = a usable `Selection`, no training.

## Pipeline

```
mesh + field ──rasterize──► field_raster  (HxW, multi-channel) ──► SAM2(click) ──► mask_raster (HxW bool)
mesh         ──rasterize──► elemid_raster (HxW int label)      ──project────────► Selection (element ids)
```

- **field_raster** — one channel per signal (bathymetry / curvature / size-fn),
  normalized to uint8, fed to SAM2 as a pseudo-image.
- **elemid_raster** — paint each element polygon with its element id
  (`skimage.draw.polygon`) on the SAME grid. Gives exact back-projection, no
  centroid-sampling loss.
- **project** — element selected iff ≥ 50% of its pixels fall inside the SAM2 mask.
- **click map** — mesh `(x, y)` → pixel via a bbox affine transform.

## Phases

- **P0 — CPU, no model.** Build rasterize + elemid-label + back-project. Stub model =
  `skimage` watershed / flood-fill from the click on `field_raster`. Proves plumbing
  and sets the **deterministic baseline IoU**. No GPU, no checkpoint.
- **P1 — real SAM2.** Swap stub → SAM2 image-predictor with a point prompt. Checkpoint
  via `huggingface_hub`. Same metrics.
- **P2 — decision.** Resolution sweep + click-jitter robustness + field-channel
  ablation → apply the gate.

## Fixtures

- **Synthetic 2-basin (primary)** — deform a structured grid + analytic depth (two
  gaussians). Ground-truth region = a known basin → exact IoU target. Fully
  controllable; build this first.
- chilmesh `annulus` / `donut` — plumbing sanity only.
- One real (WNAT + bathymetry) if reachable; skip on 403.

## Tests / metrics

| test | assertion |
|---|---|
| `test_roundtrip_recall` | paint → mask(ALL) → project recovers ≥ 0.95 of elements @ 512 (res adequacy) |
| `test_iou_vs_gt` | SAM2 click-in-basin: IoU(pred, true_basin) ≥ 0.70 on ≥ 2 cases |
| `test_sam2_beats_baseline` | IoU(SAM2) > IoU(P0 watershed) — **the real question** |
| `test_jitter_robust` | click ± 8 px, pairwise IoU ≥ 0.80 |
| `test_res_sweep` | IoU vs {256, 512, 1024} → min stable res |
| `test_channel_ablation` | bathy vs + curvature vs + size-fn → does an extra channel help |

## Success gate

All true → v2 rasterize-SAM2 demo is viable:

- roundtrip recall ≥ 0.95 @ 512
- SAM2 IoU ≥ 0.70 on ≥ 2 cases
- jitter IoU ≥ 0.80
- **SAM2 IoU > baseline IoU** — if this fails, SAM2 adds nothing; ship deterministic
  region-grow and kill the SAM2 track.

## Files

```
prototype/
  raster.py          mesh → field_raster + elemid_raster
  project.py         mask_raster → element Selection
  model_stub.py      P0 watershed baseline (same interface as SAM2 wrapper)
  model_sam2.py      P1 SAM2 wrapper (HF checkpoint)
  fixtures.py        synthetic 2-basin generator + ground-truth
  run_prototype.py   CLI: mesh + click → Selection + metrics
  tests/test_*.py    the metrics above
  REPORT.md          IoU table per phase → gate verdict
```

## Deps (isolated extra `[proto]`)

`numpy scipy scikit-image` (P0) · `torch sam2 huggingface_hub` (P1). SAM2-tiny on CPU
is slow but demo-fine.

## Risks (each caught by a test)

- Tiny coastal elements < 1 pixel → resolution floor (`test_res_sweep`).
- Scalar-field raster is out-of-distribution vs SAM2's natural-image training → may
  segment garbage. **This is the core risk; `test_iou_vs_gt` answers it.**
- Click → pixel off-by-one (`test_jitter_robust` + roundtrip).

## Start order

P0 first (cheap, CPU, no model) — synthetic fixture + `raster.py` + `project.py` +
`test_roundtrip_recall`. Eyeball numbers, then green-light P1.
