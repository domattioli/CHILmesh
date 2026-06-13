# mesh-segmenter Context

Glossary seed for the proposed **mesh-segmenter** package — interactive, composable
sub-region *selection* over a 2D mesh ("SAM2 for meshes"). It selects; it never
generates. Lives in CHILmesh `docs/proposals/` until the sibling repo is created,
then lifts in whole. Decisions recorded in [`ADR-0001-architecture.md`](./ADR-0001-architecture.md).

## Language

**Selection**:
The canonical object — an immutable *mask over elements* bound to a parent mesh
(`int64` ids / `bool[n_elements]`). The analogue of a SAM2 mask. Set-algebra
(`|`, `&`, `~`) returns new Selections; everything else is a derived export.
_Avoid_: Region, mask (in user-facing API), subset.

**Element**:
A mesh cell (triangle or quad) — the unit a Selection is over. Canonical entity:
all mechanisms return element masks; nodal signals reduce onto elements.
_Avoid_: face, cell, triangle (when quads are also in play).

**Mechanism**:
A function that produces or refines a Selection — the analogue of a SAM2 prompt.
v1: `grow`, `by_distance`, `by_click`, `by_threshold`, `by_polygon`,
`by_skeleton_layer`.
_Avoid_: selector, filter, prompt (reserve "prompt" for the SAM2 analogy in prose).

**grow**:
Ring expansion / morphological *dilation* — expand a seed set outward by `n_rings`
of dual-graph adjacency. This is what the original issue called `by_layer`.
_Avoid_: by_layer, ring (the issue's collided name — see Flagged ambiguities), expand.

**Reduction rule**:
How a per-node signal (bathymetry, curvature) collapses onto the canonical element
mask. Default conservative *"all vertices in range"*; opt-in `any` / `mean`.
_Avoid_: aggregation, projection.

**Neutral artifact**:
What a Selection emits at the re-mesh seam, carrying no generator dependency:
`element_ids`, `boundary` (raw polygon rings, holes supported), and an optional
`submesh`. The umbrella — not the segmenter — wraps rings into an `admesh.Domain`
or feeds the submesh to quadmesh.
_Avoid_: output, result.

**Engine**:
The topology provider. mesh-segmenter depends on **chilmesh** only — adjacency
(`Edge2Elem`, `Vert2Elem`), `elements_in_layer`, `submesh`, spatial indices — built
from raw `(nodes, elements)`. Never depends on a *generator* (admesh / quadmesh).
_Avoid_: backend (reserve for chilmesh's C++/Rust compute backends).

**Umbrella**:
The proposed future "one-stop" mesh package that depends on chilmesh, admesh,
quadmesh, valence **and** mesh-segmenter, and wires `selection → re-mesh`. It owns
the selection→generator handoff; the segmenter stays a leaf.
_Avoid_: orchestrator, pipeline, one-stop (informal only).

## Flagged ambiguities

- **`Layer` is reserved.** CHILmesh `CONTEXT.md` binds **Layer** = a medial-axis
  *skeletonization peel* (OE/IE/OV/IV, global, inward). The issue's
  `by_layer(ring=…, n_layers=N)` meant *ring expansion outward from a seed* —
  a different operation. Resolution: that operation is **`grow`** (dilation);
  CHILmesh's true peels are exposed separately as **`by_skeleton_layer(idx)`**.
  Never let "layer" name the ring-expansion mechanism.

- **`Mesh` is overloaded** (inherited from CHILmesh/ADMESH-Domains). Here a
  Selection's *parent mesh* is a runtime topology object (a `CHILmesh`). The thin
  `admesh.Mesh` wire dataclass is an *input* that gets built into a `CHILmesh` for
  adjacency. Say "parent mesh" for the runtime object.

- **`node` vs `Element`.** fort.14 / ADCIRC say "node"; a Selection is over
  *elements*. Nodal fields exist (bathymetry per node) but never form the mask
  directly — they reduce. Keep the I/O word ("node") out of the selection API.

## Example dialogue

> **Dev:** "User clicks in the Gulf of Mexico — do we return the nodes or the
> triangles inside?"
> **Domain expert:** "Elements. A Selection is always an element mask. The click is
> a *mechanism* — `by_click` flood-fills connected elements until a predicate stops
> it."
> **Dev:** "But bathymetry is per node. How does 'shallower than 2 m' become elements?"
> **Domain expert:** "Through the *reduction rule*. Default: an element is in only if
> *all* its vertices are under 2 m. So `by_threshold` reads the nodal field, reduces,
> and hands back an element Selection — same type as every other mechanism."
> **Dev:** "Then the user wants to re-mesh that. We call `admesh.triangulate`?"
> **Domain expert:** "Not from in here. The Selection emits a *neutral artifact* — the
> boundary rings. The umbrella turns rings into an `admesh.Domain` and re-meshes. The
> segmenter never imports a generator; that's how it stays a chilmesh-only leaf."
> **Dev:** "And expanding three rings off the coastline — that's `by_layer`?"
> **Domain expert:** "Call it `grow`. 'Layer' is CHILmesh's skeletonization peel —
> different thing. `grow(seed, n_rings=3)` is dilation. The peels are
> `by_skeleton_layer`."
