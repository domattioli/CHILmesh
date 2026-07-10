# Quickstart: Regenerate & Verify the README Hero GIF

How to regenerate `docs/gallery/readme_pipeline_annulus.gif` and verify every success
criterion after implementing feature 002.

## Prerequisites

```bash
cd /home/user/CHILmesh
source .venv/bin/activate            # venv: chilmesh 1.3.1, numpy, scipy, matplotlib, Pillow
```

If the venv is missing, build it per the repo's dev setup (`pip install -e ".[dev]"`).

## Regenerate the GIF (~2 min)

```bash
MPLBACKEND=Agg python scripts/generate_hero_animation.py
```

Always set `MPLBACKEND=Agg` (headless render). The script prints progress
(`Boundary rings…`, `Running ADMESH truss…`, `Rendering N frames…`, `Done: … (NNN KB)`).
The hero run enables `snapshot_retriangulate=True, snapshot_strict_interior=True` on
`distmesh2d_warmstart`; the ~200 extra Delaunay calls fit inside the ~2 min budget.

## Run the solver-kwarg tests

```bash
MPLBACKEND=Agg pytest tests/test_admesh_warmstart.py -v
```

Must pass, including the new asserts:
- `TestHistoryCapture` — default-path byte-identity (new kwargs at defaults → output
  identical to a no-hook run); existing `test_history_none_is_default_noop` still green.
- `TestSnapshotRetriangulate` — with `snapshot_retriangulate=True`, every captured `t` is a
  valid Delaunay of its own `p` (centroid `fd < 0`), zero inverted triangles, boundary rows
  pinned; and the returned `(p_out, t_out)` equals the hook-off return (physics unchanged).

All pre-existing tests in the file must pass **unmodified** (backward-compat gate).

## Success-criteria checklist (SC-001 … SC-007)

| SC | What it asserts | How to verify |
|---|---|---|
| **SC-001** | Zero truss-stage frames with a hole-spanning or grossly oversized/overlapping element. | Frame-by-frame audit of the rendered GIF with corrected metrics: `hole-span = edge-midpoint fd > +0.02`, `oversize = area > 8× snapshot-median`, `inverted = signed-area sign flip vs majority`. Probe 1 measured 0/200 frames degenerate under the shipped kwargs. Manual: step the GIF / re-run the probe audit script. |
| **SC-002** | Mean per-frame node displacement varies ≤ 3× between slowest/fastest truss thirds. | Scripted: recompute the thirds ratio over the play schedule (Probe 2 measured 1.261). A guard can assert `ratio <= 3`. |
| **SC-003** | Each stage boundary holds ≥ 1.5 s (≥ 15 frames @ 10 fps). | Inspect the schedule: all inter-stage holds are 18 f (1.8 s, per analyze-cycle-1 F-13); the peel k=0 hold is 18 f (FEM→Peel). Confirm in the rendered GIF's held runs. |
| **SC-004** | Peel reveals every layer once boundary-first, ≥ 0.8 s between reveals; histogram per-bin totals bit-identical across all peel frames while colors convert. | Step peel frames: k=1..4 held 8 f (0.8 s) each; assert `np.array_equal` of stacked per-bin totals across k=0..4 (Probe 3 = True). |
| **SC-005** | 100% of histogram-bearing frames carry both color-keyed lines + matching colored texts. | Sample frames from every stage: green dotted Median line + green text, red dotted Mean line + red text, neutral Min text present each time. |
| **SC-006** | GIF ≤ 4.2 MB and animates in the GitHub README. | `ls -l docs/gallery/readme_pipeline_annulus.gif` → ≤ 4.2 MB (post-process 256c+optimize projects ~2.3–3.0 MB). Confirm dims 864×396 (`python -c "from PIL import Image; print(Image.open('docs/gallery/readme_pipeline_annulus.gif').size)"`). Manual: view on the pushed `development` branch's README. |
| **SC-007** | Two runs produce structurally identical animations (same frame count + schedule). | Run the generator twice; compare printed frame count + schedule. Determinism holds under seed 11 (no new RNG; Delaunay + ADAPTIVE palette deterministic on fixed input). |

## Binary push rule (do NOT skip)

- The regenerated GIF is a **binary** file. Push it with the **git CLI only** — never via the
  GitHub MCP tools (`mcp__github__create_or_update_file` / `push_files`), which double-encode
  base64 and silently corrupt binaries (repo CLAUDE.md + DomI #85).
- Push to **`development` only** (repo single-branch policy). Do not create a feature branch.

```bash
git checkout development
git add docs/gallery/readme_pipeline_annulus.gif \
        scripts/generate_hero_animation.py \
        src/chilmesh/_vendor_admesh_truss.py \
        tests/test_admesh_warmstart.py
git commit -m "feat: refine README hero GIF (degeneracy, pacing, peel reveal, metrics)"
git push origin development
```

If commit signing fails in the cloud env, use `git -c commit.gpgsign=false commit` (repo
CLAUDE.md workaround). Confirm the pushed binary is intact (non-empty, GIF magic bytes) via
the raw URL rather than an MCP read.
