# Distance, medial axis, skeleton, layers — four related, distinct constructs

These four terms are routinely conflated (CHILmesh itself shipped a layerization
routine misnamed `_skeletonize` until #221). They form a chain — each derives
from the one before — but they are **not** interchangeable. This note fixes the
vocabulary used throughout CHILmesh and explains the differences at a high,
algorithmic, and mathematical level.

> See also: the four-panel figure `docs/gallery/mesh_concepts.png`
> (`scripts/illustrate_mesh_concepts.py`) and the ownership decision in
> [#221](https://github.com/domattioli/CHILmesh/issues/221).

## Quick reference

| Construct | What it is | Output type | Common synonyms | CHILmesh |
|---|---|---|---|---|
| **Distance field** | scalar: distance from each interior point to the boundary | continuous field over the domain | distance transform, EDT, wall distance, `d(x)`; *signed* variant = SDF | not yet (roadmap; lives in ADMESH for generation) |
| **Medial axis** | the ridge of the distance field: centers of maximal inscribed disks | continuous 1-D curve/graph (the domain's "spine") | MAT (medial axis transform), symmetric axis, grassfire locus | not directly |
| **Skeleton** | a 1-wide, connectivity-preserving thinning of the shape | discrete 1-px / 1-element curve | thinning, morphological skeleton, centerline | `skeletonize()` (thinning peel) |
| **Layers** | the mesh's elements grouped into concentric bands peeled boundary-inward | partition of elements (integer layer index) | onion peeling, concentric layering, advancing-front layers | `_layerize()` → `.layers`, `.n_layers` |

The chain: **distance is a field → its ridge is the medial axis → the skeleton
is a thinned discrete approximation of that ridge → layers are concentric bands
of which the innermost approximate the medial region.**

## Distance field

- **High level.** A scalar `d(x)` defined everywhere inside the domain: how far
  is `x` from the nearest boundary point. Signed variant (SDF) is negative
  outside, positive inside (or vice versa).
- **Math.** `d(x) = min_{b ∈ ∂Ω} ‖x − b‖`. The SDF additionally carries the sign
  of inside/outside.
- **Algorithm.** Exact Euclidean distance transform (EDT), or the Fast Marching
  / Eikonal solution of `‖∇d‖ = 1`. `O(n)`–`O(n log n)` on a grid.
- **Note.** In this ecosystem a *domain* distance/SDF (pre-mesh, drives element
  size during generation) belongs to **ADMESH**, not CHILmesh (#221).

## Medial axis

- **High level.** The shape's central spine — the set of points that are
  equidistant from two or more nearest boundary points. Equivalently, the
  centers of all maximal inscribed disks.
- **Math.** `MA(Ω) = { x ∈ Ω : |argmin_{b ∈ ∂Ω} ‖x − b‖| ≥ 2 }`. It coincides
  with the **ridge set of the distance field** (where `∇d` is discontinuous).
  The *medial axis transform* (MAT) pairs each axis point with its inscribed
  radius `d(x)` — enough to reconstruct the shape exactly.
- **Algorithm.** Distance-ridge detection, Voronoi-diagram of the boundary
  sampling, or the grassfire (wavefront-collision) model.
- **Why it differs from the skeleton.** The medial axis is a *geometric* object
  defined by the continuous shape; it can be sensitive to tiny boundary
  perturbations (spurious branches) and need not be a clean 1-pixel curve.

## Skeleton

- **High level.** A *thinned* version of the shape: repeatedly shave the
  boundary inward, but never remove anything whose removal would break the
  shape's connectivity, until nothing more can be shaved. What remains is a
  1-wide connected curve.
- **Math.** A subset `S ⊆ Ω` that is **homotopy-equivalent** to `Ω` (same
  connectivity / holes) and is locally thin (width 1). Produced by removing
  *simple points* — points whose deletion preserves topology — until the result
  is idempotent.
- **Algorithm.** Iterative morphological thinning (e.g. Zhang–Suen): each pass
  deletes simple boundary points; repeat until a pass changes nothing.
- **Skeleton vs. medial axis.** Both are 1-D "central" descriptors and are often
  used interchangeably, but they are produced by different procedures and can
  differ: the medial axis is the exact distance ridge (may branch to every
  corner); the skeleton is whatever a topology-preserving thinning leaves
  behind (cleaner, connected by construction). See panels 2 vs 3 of the figure.

## Layers (layerization)

- **High level.** Take the actual mesh and peel it like an onion: the outermost
  ring of boundary-touching elements is layer 0, the next ring in is layer 1,
  and so on until no elements remain. Every element gets a layer index.
- **Math.** Define the layer index `ℓ(e)` of element `e` as the iteration at
  which it is peeled, where each iteration removes **all** elements currently
  adjacent to the (receding) boundary. `n_layers = max_e ℓ(e) + 1`. This is a
  *quantized, mesh-discrete analogue of the distance field* — `ℓ(e)` tracks
  graph-distance-to-boundary in elements, not Euclidean distance.
- **Algorithm.** CHILmesh `_layerize()` (formerly mis-named `_skeletonize`):
  vectorized boundary-ring removal, `O(n)` in element count, storing
  `OE`/`IE`/`OV`/`IV` (outer/inner elements and vertices) per layer.
- **Why it is not skeletonization.** Layerization removes **every** boundary
  element each pass with no topology test, so the active set shrinks to empty
  and the output is a full *partition into bands*. Skeletonization removes
  **only** elements whose deletion preserves connectivity, so the active set
  shrinks to a thin *spine* and stops there.

## The two mesh operations: `_layerize` vs `skeletonize`

Both are iterative inward peels on the same mesh; the **peel rule** is the only
difference.

| | `_layerize()` | `skeletonize()` |
|---|---|---|
| Peel rule | remove **all** current boundary elements each pass | remove only **removable** (connectivity-preserving) boundary elements |
| Topology | not preserved (peels through to nothing) | preserved (homotopy-equivalent thinning) |
| Terminates when | no elements remain (empty) | no element can be removed without disconnecting → the medial spine remains |
| Output | layer index per element; concentric bands | skeleton elements + the peel order that produced them |
| Discrete analogue of | the distance field (banded) | the medial axis (thinned) |

So `skeletonize()` *mimics* the layerize machinery — peel, record, repeat until
you cannot peel — but with the **skeletonization peel** (drop only simple,
topology-safe elements) rather than the **layerize peel** (drop the whole ring).

## Synonyms seen in the wild

- Distance field ≈ distance transform ≈ EDT ≈ wall distance ≈ `d(x)`; signed → SDF.
- Medial axis ≈ MAT ≈ symmetric axis ≈ (loosely) "topological skeleton".
- Skeleton ≈ thinning ≈ morphological skeleton ≈ centerline.
- Layers ≈ onion peeling ≈ concentric layering ≈ advancing-front layers ≈
  (the historical CHILmesh misnomer) "skeletonization".

When in doubt: **what is the input and what is the output?** A *domain* in,
a *field* out → distance. A field's ridge → medial axis. A topology-preserving
1-wide thinning → skeleton. A *mesh* in, an *integer band per element* out →
layers.
