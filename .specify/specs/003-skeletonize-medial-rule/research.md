# Research: the precise peel rule for `skeletonize()`

**Issue**: [#223](https://github.com/domattioli/CHILmesh/issues/223) — Implement a real
`skeletonize()`, distinct from the layer peel (`_peel`, formerly `_layerize`/`_skeletonize`).
**Status**: research note (advance-only; no implementation until the operator ratifies a rule).
**Author of the decision**: the operator (this is the CHILmesh author's specialized algorithm; a
routine session may *specify* the options but must not *pick* the semantics).

## Why this note exists

#223 is blocked on a single decision — **which peel rule `skeletonize()` implements** — and has
been unanswered for 5+ weeks (flagged in the #261 zoom-out as a decision-limbo bottleneck). The
issue offers three candidate rules; the author leans toward option (1) but explicitly asked for
confirmation "before implementation so we don't ship a third wrong version."

This note does three things so the decision becomes cheap:

1. Surfaces a **contradiction inside `docs/CONCEPTS.md`** about what `skeletonize()` even *is* — the
   real reason the decision is hard.
2. Specifies option (1) concretely against the *current* code (`_peel`, post-#187), so it is
   ratify-or-redirect rather than re-derive.
3. States a validation plan and the naming consequences of each choice.

It does **not** add code, and it does not flip any label on its own.

## The blocking contradiction: `CONCEPTS.md` defines `skeletonize()` two incompatible ways

`docs/CONCEPTS.md` is the repo's canonical explainer for the four constructs (distance field →
medial axis → skeleton → layers). It disagrees with itself on `skeletonize()`:

- **Line 19** (the reserved-name table row) says: *"`skeletonize()` — reserved, **medial-axis only**
  (#223, unimplemented)."* → this is **option (1)**, a *geometric* medial-axis extraction.
- **Lines 90–105** (`## The two mesh operations: _peel vs skeletonize`) describe `skeletonize()` as
  *"remove only removable (connectivity-preserving) boundary elements ... homotopy-equivalent
  thinning ... the medial spine remains"* → this is **option (2)**, a *topological* skeleton
  (Zhang–Suen-style thinning on the face complex).

The same document's own `## Skeleton` vs `## Medial axis` sections (lines 40–69) are explicit that
these are **different constructs produced by different procedures**: the medial axis is the exact
distance ridge (may branch to every corner; can be disconnected/spurious), while the skeleton is a
topology-preserving thinning (connected by construction, 1-wide). So the two `skeletonize()`
definitions in CONCEPTS.md are not two phrasings of one idea — they name two different outputs.

**This is why the three prototypes in #223 each failed one way or another:** they were chasing a
target the docs never pinned down. The decision #223 needs is not "which of three algorithms is
least buggy" — it is **"which construct does the reserved name `skeletonize()` denote?"** Everything
else follows.

## Current code the implementation would build on (post-#187)

`_peel()` (`src/chilmesh/CHILmesh.py:1128`) is a vectorized inward ring removal. Per pass `iL` it
records, in `self.layers`:

| key | contents |
|---|---|
| `OV[iL]` | outer-ring vertices bounding the layer's outer side |
| `OE[iL]` | elements adjacent to the layer's boundary edges |
| `IE[iL]` | active elements adjacent to any edge touching an `OV[iL]` vertex |
| `IV[iL]` | vertices of `OE ∪ IE` not in `OV[iL]` |
| `bEdgeIDs[iL]` | boundary edges defining this layer's outer frontier |

`n_layers = max_e ℓ(e) + 1` is the graph-distance-to-boundary in elements (a quantized,
mesh-discrete analogue of the distance field). The public verb is `peel_layers()`; the medial-axis
name `skeletonize` is reserved and **not even a stub** today (the `_skeletonize` compat alias
referenced in the #223 body was dropped by the #187 rename). So there is no back-compat surface to
preserve — `skeletonize()` is a clean-slate addition.

The important reuse fact: `_peel` already computes, for every element, the pass at which the
advancing front reached it. That *is* the grassfire arrival time. Option (1) below is essentially
"detect where fronts arriving from different boundaries meet," and all the arrival-time data it needs
is the layer decomposition `_peel` already produces.

## Option (1) — front-collision medial extraction (grassfire), specified

**Construct produced:** the medial axis (option (1) in #223; matches CONCEPTS.md **line 19**).

**Intuition:** the medial axis is the grassfire wavefront-collision locus (CONCEPTS.md line 48) —
the points where fires lit simultaneously along the whole boundary and burning inward at unit speed
meet head-on. In the mesh-discrete setting the "fire" is exactly `_peel`'s advancing front, and the
"arrival time" of element `e` is its layer index `ℓ(e)`.

**Rule.** An element is *medial* when the peel front reaches it from **two distinct boundary
directions at (nearly) the same time** — i.e. it is a local maximum of the graph-distance-to-boundary
field, or it sits on a ridge between two receding fronts. Concretely, a candidate definition to
ratify or amend:

1. Compute `ℓ(e)` for every element via `_peel` (already done; it is the layer index).
2. Build the element dual graph (elements are nodes; share-an-edge is an adjacency —
   `Edge2Elem` gives this directly).
3. Mark element `e` **medial** if it has no dual-neighbor with strictly greater `ℓ` — i.e. `e` is a
   local maximum / ridge element of the discrete distance field. (Equivalently: the front cannot
   advance past `e` into a deeper element, so two fronts collided at `e`.)
4. Record, per medial element, the peel pass `ℓ(e)` that produced it (the #223 output contract:
   "skeleton elements + the peel order that produced them").

**Why this is the right shape for a *medial* target:** it is connected-by-front-collision on
simply-connected regions, branches toward genuine corners (as the true medial axis does), and reuses
the layer-collision data the issue itself points at (`IV`/`IE`). It is `O(n)` on top of the existing
`O(n)` peel.

**Known caveats to decide in the plan, not now:**
- The medial axis is *allowed* to branch and to be locally >1 element wide at coarse resolution — it
  is not required to be a clean 1-wide curve (that is the *skeleton*, option 2). Acceptance tests
  must grade it as a medial axis (ridge coverage), not as a thinned skeleton (1-width).
- Multiply-connected domains (holes/islands): fronts also advance outward from inner boundaries, so
  the ridge forms between outer and inner fronts. `_peel` already seeds from all boundary edges
  (including inner rings), so `ℓ(e)` is the min distance to *any* boundary — correct for this.
- Local-max ties on flat plateaus (uniform `ℓ` bands) can thicken the ridge; a tie-break
  (e.g. keep the plateau element nearest the centroid of its `ℓ`-band component) is a plan-level
  refinement, not a semantic change.

## Option (2) — topology-preserving thinning (skeleton), for contrast

**Construct produced:** the topological skeleton (option (2) in #223; matches CONCEPTS.md
**lines 90–105**).

**Rule.** Iteratively remove *simple* boundary elements (deletion preserves homotopy type) while
*retaining endpoints*, à la Zhang–Suen but defined on a 2-D tri/quad face complex. Terminates at a
1-wide connected spine.

This is a materially larger piece of work: it needs precise definitions of "simple element" and
"endpoint element" for tri/quad faces, plus thinning-mask timing to avoid the two failure modes the
#223 prototypes already hit (local-thinning → 711 disconnected components; global
connectivity-preserving → collapse to 1 element). It is the *only* option that yields the
connectivity-preserving, homotopy-equivalent object the CONCEPTS.md peel-rule table describes.

## Option (3) — innermost-layer derivation, for completeness

`skeleton = elements in the deepest k layers`. The #223 prototype showed this is clean and connected
(34/2276 elems, 1 component on the L-shape) — but it is literally `_peel`'s *output sliced*, not a
distinct peel with its own rule. It cannot honor the issue's framing ("a distinct skeletonization
peel") and adds no construct the layer index doesn't already give. Include only as the trivial
fallback if the operator decides `skeletonize()` should not be a separate algorithm at all.

## Recommendation (routine session — advisory only)

Ratify **option (1)** and, in the same pass, **fix the CONCEPTS.md contradiction to match**:

- It aligns with CONCEPTS.md **line 19**'s existing reserved-name definition ("medial-axis only").
- It reuses the front-collision data the issue itself identifies, is `O(n)` on top of `_peel`, and is
  connected-by-construction on simply-connected regions.
- It is far less code than option (2) and — unlike option (3) — is a genuine distinct rule.

If instead the operator wants the **homotopy-preserving, 1-wide skeleton** described in CONCEPTS.md
lines 90–105, that is **option (2)**, a different construct, and this note's specification does not
apply — a separate thinning spec (simple-element + endpoint definitions) would be needed. Either way,
**CONCEPTS.md must be made self-consistent** as part of #223: the reserved-name row and the peel-rule
table currently promise two different objects under one name.

## Validation plan (applies once a rule is ratified)

- **Fixtures:** annulus, donut (has a hole), block_o, structured, plus the L-shape used in the #223
  prototype table — cover simply- and multiply-connected, tri and quad.
- **For option (1) (medial):** assert ridge elements are local maxima of `ℓ`; assert every interior
  element is within one dual-hop of a medial element (ridge coverage); assert the medial set is
  connected per connected component of the domain (allowing branches). Do **not** assert 1-width.
- **For option (2) (skeleton):** assert the result is 1-element-wide, connected, and
  homotopy-equivalent (same component + hole count as the input via Euler characteristic).
- **Both:** `skeletonize()` must not mutate `self.layers` / `_peel` output; a `peel_layers()` call
  before and after `skeletonize()` returns identical layers.

## What the operator needs to decide (one line)

Pick the construct the reserved name `skeletonize()` denotes: **(1) medial axis** (front-collision,
spec'd above — recommended), **(2) topological skeleton** (thinning; needs its own spec), or
**(3) innermost-layer slice** (not a distinct algorithm). The CONCEPTS.md contradiction is fixed to
match whichever is chosen.

## References

- #223 (issue), #221 (the `_layerize` rename), #187 (rename to `_peel`, dropped `_skeletonize` alias)
- `docs/CONCEPTS.md` — the four-construct explainer (and the contradiction this note flags)
- `src/chilmesh/CHILmesh.py:1128` (`_peel`), `:1251` (`peel_layers`)
- Blum (1967) grassfire/MAT; Zhang & Suen (1984) thinning
