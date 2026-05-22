# CHILmesh Context

Computational engine for 2D triangular, quadrilateral, and mixed-element
meshes. Loads, analyzes, and writes mesh files in formats sourced from
ADMESH-Domains and elsewhere. CHILmesh is a consumer of the ADMESH-Domains
glossary; it does not define `Domain`, `Mesh`, `Identity`, or related
terms itself.

## Upstream Terms

The following terms are defined canonically in ADMESH-Domains and must be
used as defined there:

- `Domain`
- `Mesh`
- `Mesh Boundary` vs `Domain Boundary`
- `Identity`
- `Canonical form`
- `Containment`
- `Relatedness`
- `Mesh Record`
- `Curator`

See: https://github.com/domattioli/ADMESH-Domains/blob/main/CONTEXT.md

CHILmesh's I/O implementations must respect the canonicalization rules
defined upstream. Any drift between CHILmesh's serialization output and
the upstream Canonical-form contract is a bug, not a design choice.

## Language

**CHILmesh** (the class):
The in-memory computational mesh object: vertices, elements, precomputed
adjacency structures (vert↔vert, vert↔elem, elem↔elem, edge↔elem), and
optional Layer skeletonization. Distinct from a `Mesh Record` (ADMESH-
Domains catalog entry) — CHILmesh is the runtime representation, the
record is the catalog row that points at the file the runtime loads.
_Avoid_: "the mesh" in code paths where ambiguity with `Mesh Record` is
possible — say `CHILmesh` or `Mesh Record` explicitly.

**Vertex**:
A 2D point in a CHILmesh, indexed and carrying optional bathymetry. We
call it `Vertex` everywhere except inside fort.14 I/O code paths, where
the ADCIRC term `node` is used to match the format's own vocabulary.
_Avoid_: node (only allowed in fort.14 I/O modules).

**Element**:
A connectivity group of 3 (triangle) or 4 (quad) vertices.
_Avoid_: cell (FVM term), face (overloaded with 3D faces).

**Element Type**:
`triangle` | `quadrilateral` | `mixed`. A property of a CHILmesh instance,
derivable from element connectivity.
_Avoid_: shape, kind (`kind` is the ADMESH-Domains `Mesh Record` field
meaning `"mesh"` vs `"boundary"` — different concept).

**Edge**:
An ordered pair of vertices shared between elements. Stored in the
internal half-edge structure `EdgeMap`.

**EdgeMap**:
Internal half-edge adjacency structure used by topology queries.

**Layer**:
A nested boundary ring produced by Skeletonization. Layer 0 = outermost
elements; layer N = innermost (medial axis approximation).
_Avoid_: ring, shell.

**Skeletonization**:
Recursive boundary-peeling that yields Layers. Triggered by
`_skeletonize()`. Gated by the `compute_layers` flag at construction.

**compute_layers**:
Boolean flag (default `True`) on readers and constructors. `False` skips
Skeletonization at load time — fast init for boundary-only meshes or very
large meshes. Trades startup time for Layer-query availability.

**Reader** / **Writer**:
- `read_from_fort14(path)` — ADCIRC `.fort.14`. Triangles + quads.
- `read_from_2dm(path)` — SMS `.2dm`. Triangles + quads.
- `write_fort14(mesh, path)` — ADCIRC `.fort.14` writer.
- `.chil` reader / writer — TBD; design in progress in ADMESH-Domains.

**from_admesh_domain(record)**:
Bridge entry point: takes an ADMESH-Domains `Mesh Record`, routes by its
`type` field to the right Reader. See `specs/002-admesh-domain-loader/`.
_Avoid_: load, fetch (overloaded with `Mesh Record.load()` on the
ADMESH-Domains side, which downloads bytes — different operation).

## Flagged ambiguities

**`Mesh` vs `CHILmesh` vs `Mesh Record`** — when the distinction matters,
use the qualified term. `CHILmesh` for the runtime object, `Mesh Record`
for the catalog entry, plain `Mesh` only when speaking abstractly (and
upstream-defined per ADMESH-Domains).

**`kind`** — overloaded across repos. Inside CHILmesh, `kind` is unused;
do not introduce. Inside ADMESH-Domains, `Mesh Record.kind` is
`"mesh"` vs `"boundary"`. Inside fort.14 boundary segments, the ADCIRC
`IBTYPE` integer codes serve a similar role but are NOT called `kind`.

**`node`** — fort.14 vocabulary. Use only in fort.14 I/O code paths.
Elsewhere: `Vertex`.

## Example dialogue

> Dev: "Why does the reader take `compute_layers=False`?"
>
> Reviewer: "Big meshes — WNAT-class — pay a Skeletonization cost on load
> that boundary-only consumers don't need. False skips it; layer queries
> raise until the user explicitly opts in. Identity is unaffected — it's
> derived from the `Mesh Record` and the file bytes, not the runtime
> CHILmesh's Layers."
