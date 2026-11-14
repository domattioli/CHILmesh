# Spec 008: v1.7.0 release prep

Branch: `release/v1.7.0` (off `development` @ 2cc4be9). Status: draft.

## Scope

1. #309: fort.14 open-boundary IBTYPEE survives read, write, re-read.
2. #306: `angle_based_smoother` never accepts a move that folds any quad corner.
3. #307: record a decision on NOAA validation data storage (no code).
4. Cut v1.7.0 from `development` after the operator settles the licence route.

Out of scope: `.chil` work (#201; next version), #301 harness (branch `feat/301-competitor-benchmarks`), Codespace-blocked #276 and #310.

## Requirements

- FR-1 (#309): `_read_fort14_open_segments` keeps the second field of each NVDLL line as `ibtype` (int), or `None` if absent.
- FR-2 (#309): the fort.14 writer emits `NVDLL IBTYPEE` when `ibtype` is not None, and `NVDLL` alone otherwise. Flow boundaries unchanged.
- FR-3 (#309): a round-trip test with a non-default open IBTYPEE compares equal.
- FR-4 (#306): `_local_min_quality_fast` scores a quad -1 when any corner i fails `cross(p[i+1]-p[i], p[i-1]-p[i]) > SJ_MIN*|u|*|w|`, SJ_MIN = 1e-6.
- FR-5 (#306): a test with the dart quad from #306 (p0=(1,0.5), p1=(0,0), p2=(2,0), p3=(1,2)) shows the check rejects it; a jittered 7x7 grid smoothing run creates no new inverted quads.
- FR-6: triangle behaviour and the C++ backend fallback stay unchanged; full pytest passes.
- FR-7: CHANGELOG entry for 1.7.0 lists #308, #304, #309, #306, and the Python 3.10 drop.

## Clarifications
- Q #306 default → on by default (bug fix). CHANGELOG notes quad smoothing output may differ.
- Q #307 → Valence, later. Close #307; file Valence issue for non-mesh assets; harness keeps direct fetch meanwhile.
- Q #309 absent IBTYPEE → keep absent (write `NVDLL` alone). General opt-in field-report/lint → new issue, not 1.7.0.
- Q #301 → not in 1.7.0. Whether #301 + `.chil` become 2.0.0 or two 1.x releases is decided later.

## Addendum: progress line (undocumented feature, operator 2026-11-14)

- FR-8: new private module `chilmesh/_progress.py`. One stderr line, redrawn in place with `\r`, cleared on finish, so stdout stays clean.
- FR-9: while running, the line uses real verbs plus elapsed time and the op's own counter, e.g. `✢ Peeling layers… (12s · layer 41)`. The glyph animates on a background thread tick. Show `%` only when the total is known.
- FR-10: on finish, the line becomes `✓ Chillin'. Peeled 87 layers in 2m 46s`. "Chillin'" appears only in the done state.
- FR-11: auto mode by default: print only when stderr is a TTY and the op has run longer than 1 s. `CHILMESH_PROGRESS=0|1|auto`, plus a private kwarg `_progress`. No public kwarg; public API unchanged.
- FR-12: wired ops: layerize, smoothing (angle-based, size-weighted), fort.14 read/write, quality analysis.
- FR-13: stdlib renderer by default. Use tqdm when it is installed; tqdm stays an optional dependency, never required.
- FR-14: undocumented, so no README entry. CHANGELOG gets one line. Tests cover: silent when not a TTY, silent under 1 s, env var forces on and off, the done line contains "Chillin'", and the tqdm and no-tqdm paths.
