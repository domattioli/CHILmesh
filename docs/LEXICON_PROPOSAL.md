# Layer-lexicon decision (#187) — RATIFIED 2026-07-10

Status: **DECIDED + swept.** Operator ratified names 2026-07-10; rename shipped same day. No consumers yet → clean rename, no compat aliases.

Caveman-compressed per operator (2026-07-10). Prior full-prose proposal in git history.

---

## Problem (was)

One concept, 3 drifting nouns (layer/peel/front) + 1 overloaded verb. `_skeletonize()` = faithful port of MATLAB `meshLayers` (inward ring-peel) but name implied medial-axis skeletonization → different construct. Reader can't tell rings from spine.

## Decision

Math frame: one inward ring = discrete level set of graph-distance-to-boundary. Iterative construction = **onion peeling** (Chazelle 1985 convex/onion layers). Fronts-collide locus = medial axis/straight skeleton = distinct op.

| Axis | Pick | Why |
|---|---|---|
| Noun | `layer` | stored ring = level set = onion/convex layer. Already dominant (290 uses). |
| Verb | `peel` | term of art = "onion peeling". Public entry `peel_layers()`. |
| `skeletonize` | reserved → medial-axis only (#223) | layers = level sets; skeleton = front-collision locus. Distinct → distinct names. |

`front` rejected as stored-ring noun → names the transient propagating interface, not the recorded ring (+ collides w/ advancing-front meshing in QuADMESH/MADMESHR). `onionize` = correct image but informal → term "onion peeling" kept in prose, not a code symbol.

## Rename map (shipped)

| Was | Now | Kind |
|---|---|---|
| `CHILmesh._layerize()` | `CHILmesh._peel()` | private |
| — | `CHILmesh.peel_layers()` | new public verb → returns `.layers` |
| `_skeletonize = _layerize` alias | **deleted** | no consumers → no shim |
| `mutations.reskeletonize_local()` | `mutations.repeel_local()` | public |
| `_skeletonize()` call sites (mutations.py) | `_peel()` | — |
| `layers` / `elements_in_layer` / `plot_layer` | unchanged | already on noun |
| (future) medial-axis op | `skeletonize()` | name freed, impl = #223 |

Test files renamed: `test_skeletonize_termination` → `test_peel_termination`, `test_reskeletonize_local_core` → `test_repeel_local_core`, `test_skeletonize_boundary_seed` → `test_peel_boundary_seed`.

## Invariant

`skeletonize` appears NOWHERE in code/tests post-sweep (grep-guarded). Reserved for #223 medial-axis. Concept distinction documented in `docs/CONCEPTS.md`.
