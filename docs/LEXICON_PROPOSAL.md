# Layer-lexicon consolidation proposal (#187)

Status: **DECISION REQUIRED** (operator picks names; rename deferred until then).
Scope: naming consistency for the inward layer-peeling concept across the codebase.
Constraint: public API is stable (constitution Principle VII) — any rename ships
with a deprecation shim, never a hard break.

---

## Problem

Three overlapping terms name one concept, and one verb is overloaded onto a
*different* concept:

- The inward ring-peeling operation is implemented as `_skeletonize()`, but it is
  a faithful port of the MATLAB **`meshLayers`** function — not medial-axis
  skeletonization. "Skeletonize" properly names medial-axis extraction, which is
  a related-but-distinct algorithm.
- Comments/internal prose drift between **layer**, **peel**, and **front** for the
  same rings. The opening post floats **"onionize"** as the verb.

Result: a reader cannot tell from the name whether `skeletonize` produces onion
rings or a medial skeleton. #187 asks: fix the noun, fix the verb, then make the
code consistent top-to-bottom — while keeping `skeletonize` / medial-axis as a
separate, honestly-named thing.

## Current usage inventory (src/, 2026-06-08)

| Term | Count | Where it lives |
|---|---:|---|
| `layer` | 290 | dominant noun. `self.layers` dict, `elements_in_layer()`, `plot_layer()`, `get_layer()`, layer indexing throughout |
| `skeletoniz*` | 48 | `_skeletonize()`, `reskeletonize_local()`, `skeletonize_diff()` — the **peel** op, mislabeled |
| `peel` | 7 | informal, comments + docstrings only |
| `front` | 1 | one docstring: "outer frontier" |
| `medial` | 1 | `layer_paths.py`: "medial axis branches" |
| `onion` | 0 | not used anywhere (only the issue title) |

Public API surface touched by a rename:
- `CHILmesh.layers` (dict; keys `OE/IE/OV/IV/bEdgeIDs`)
- `CHILmesh.elements_in_layer(layer_idx)`
- `chilplotting.plot_layer(mesh, layers=…)`
- `mutations.reskeletonize_local(...)`, `mutations.skeletonize_diff(prev_layers)`
- (private, free to rename now) `CHILmesh._skeletonize()`, `_snapshot_layers()`

MATLAB origin anchor: `_skeletonize` docstring — *"Faithful Python port of the
MATLAB `meshLayers` function."*

## Recommendation

1. **Noun = `layers`.** Already dominant (290) and matches the `self.layers`
   data structure + MATLAB `meshLayers`. Retire `peel`/`front` as concept-nouns;
   keep them only as informal verbs in prose where helpful.
2. **Verb = `layerize` (peel into layers).** Reserve **`skeletonize` exclusively
   for medial-axis extraction.** "onionize" is evocative but informal for a
   public API; "layerize" reads as the inverse of the `layers` noun and matches
   `meshLayers`. (Operator may overrule in favor of `onionize` — see decision.)
3. **Separate the concepts explicitly.** The ring-peeling op (`meshLayers` port)
   and any true medial-axis/skeleton op become distinctly named and distinctly
   documented, so neither name implies the other.

Proposed renames (deprecation-shim, not hard break):

| Today | Proposed | Notes |
|---|---|---|
| `CHILmesh._skeletonize()` | `CHILmesh._layerize()` | private — rename now, no shim needed |
| `reskeletonize_local()` | `relayerize_local()` | public → keep old name as deprecated alias |
| `skeletonize_diff()` | `layers_diff()` | public → deprecated alias |
| `layers` / `elements_in_layer` / `plot_layer` | unchanged | already on the chosen noun |
| (future) medial-axis op | `skeletonize()` | name freed up for its honest meaning |

## Migration approach (when names chosen)

- Private symbols (`_skeletonize`, `_snapshot_layers`) rename immediately.
- Public symbols get a thin alias + `DeprecationWarning` for one minor cycle
  (Principle VII); remove alias at the next minor bump.
- Single sweep PR; one `type: refactor` commit per file; tests updated in lockstep.
- `plot_layer` / `layers` / `elements_in_layer` stay put (no churn for the win).

## Decision required from operator

- [ ] Confirm noun = `layers` (vs `peel` / `front`).
- [ ] Confirm verb = `layerize` **or** override to `onionize` / other.
- [ ] Confirm `skeletonize` is reserved for the future medial-axis op only.

Once the two names are fixed, the consistency sweep is mechanical and ships as a
single deprecation-safe refactor PR. This doc is the decision artifact; no code
renamed until the boxes above are checked.
