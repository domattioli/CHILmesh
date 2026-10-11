# `.chil` Format: Codebase Integration Investigation (#201)

CHILmesh already holds about 90% of what a `.chil` writer needs for the geometry payload. The remaining 10% is metadata (CRS, domain boundary, hash) that #154's adversarial review assigns to the registry. A `.chil` adapter fits into the existing `CHILmesh.save`/`CHILmesh.load` suffix dispatch.

This document maps where the proposed `.chil` mesh container (designed in #154, verdict **RESHAPE**) would slot into CHILmesh's I/O surface. It also lists what is built and what the object model lacks.

> **Scope (per #201):** brainstorm only. Authoring or saving `.chil` files into the repo beyond examples is out of scope. This document maps the integration surface and the gaps. It does not implement a reader or writer.

**Contents**

1. Current I/O surface
2. Object model the I/O must serialize
3. Gaps for `.chil` v1
4. Where `.chil` would integrate
5. Recommendation
6. Frontmatter and split payload
7. Borrowing from `admesh/admesh`

---

## 1. Current I/O surface

| Format | Read | Write | Entry points |
|---|---|---|---|
| ADCIRC `fort.14` (`.14`/`.fort14`) | ✅ | ✅ | `CHILmesh.read_from_fort14`, `write_to_fort14`, module `write_fort14` |
| SMS `.2dm` | ✅ | ✅ (private `_write_2dm`) | `CHILmesh.read_from_2dm`, `_write_2dm` |
| Gmsh `.msh` (v2.2 + v4.1) | ✅ | ✅ | `gmsh_io.read_msh`/`write_msh`, `CHILmesh.read_from_msh`/`write_to_msh` |
| Valence registry record | ✅ | - | `CHILmesh.from_admesh_domain(record)` |
| `.chil` | ❌ | ❌ | **none; does not exist anywhere in `src/`** |

### Unified dispatch

`CHILmesh.save(filename)` and `CHILmesh.load(filename)` branch on file suffix (`CHILmesh.py:2346`, `:2354`). A `.chil` adapter is a third suffix branch. It is not a new public API shape.

The adapter pattern has precedent: `gmsh_io.py` is a standalone module that the class delegates to.

### Principle VI doc-drift

Constitution Principle VI says 2dm "write support future". `_write_2dm` already ships (private). Either promote it to a public `write_to_2dm` or correct the constitution text. This is flagged separately and is not part of this investigation.

---

## 2. Object model the I/O must serialize

A runtime `CHILmesh` carries these fields (`CHILmesh.py:274–288`):

- `points`: `(N,3)` float (z padded to 0 for 2D)
- `connectivity_list`: `(M,3|4)` int, mixed tri/quad (quads padded)
- `boundary_segments`: `list[dict]` with ADCIRC IBTYPE / kind metadata
- `boundary_condition`, `seed_ibtypes`, `seed_boundary_kinds`
- `adjacencies`: derived (Elem2Vert, Edge2Vert, Edge2Elem, …); **not serialized, recomputed on load**
- `layers`: skeletonization output (OE/IE/OV/IV/bEdgeIDs); **derived, not serialized**

The minimalism rule in #154 says to store nothing derivable. The existing formats already follow it: only points, connectivity and boundary records persist. Adjacencies and layers recompute.

A `.chil` writer therefore needs the same three things the fort.14 writer already emits, plus the new metadata in §3.

---

## 3. Gaps for `.chil` v1

`.chil` (#154) requires fields that CHILmesh does **not** currently model:

| `.chil` requirement (#154) | Present in CHILmesh? | Gap |
|---|---|---|
| Native CRS declaration (Q8) | ❌ | No CRS field; Principle V makes coords **opaque** → direct conflict (#154 D9 HIGH) |
| Domain Boundary distinct from Mesh Boundary (Q11) | ❌ | Only `boundary_segments` (mesh edges) exist; no domain-outline concept |
| `fort.13` nodal attributes for lossless mesh round-trip (Q11) | ❌ | No fort.13 read/write anywhere → "lossless fort.14+fort.13" unmet |
| Deterministic content hash / `content_uid` (Q15) | ❌ | No hashing in CHILmesh (lives in Valence schema, not here) |
| Multi-ring boundary (holes, Bermuda) (Q10/Q11) | partial | boundary_segments can list multiple rings, but no canonical winding/start-vertex normalization |
| Quantization to fixed grid (Q9) | ❌ | No coordinate quantization step |

Two of the three things a `.chil` writer needs (points, connectivity+boundary) already come from the fort.14 path. The **new** work is metadata the object does not hold (CRS, domain boundary, hash). It is not mesh geometry.

---

## 4. Where `.chil` would integrate

Integration is **localized**, once the RESHAPE blockers clear (the answer to #201's question).

1. Add the module `src/chilmesh/chil_io.py`, mirroring `gmsh_io.py`. It provides `read_chil()` and `write_chil()`, which return or consume a `CHILmesh`.
2. The container per #154 Q12 is a zip of `manifest.json` + `vertices.npy` + `triangles.npy`/`quads.npy` + `boundaries.toml` + `origin.toml`.
3. Add two suffix branches in `CHILmesh.save`/`load` (`CHILmesh.py:2346`/`2354`). `.chil` calls `chil_io.write_chil`/`read_chil`. Other call sites do not change.
4. Make `from_admesh_domain` the primary producer. Per #154 Q14, `.chil` is the registry's canonical store. The registry-bridge classmethod is where a CHILmesh first meets a `.chil` payload, so wire reading there too.
5. Reuse the existing export adapters for output: fort.14 (lossless), 2dm (lossy, drops typed boundaries) and msh. They become the "write OUT of `.chil`" adapters that #154 Q14 describes, with no new code.

### Blockers (from #154 RESHAPE, CHILmesh side)

- Principle V (coord-agnostic) conflicts with WGS84-canonical Identity. The constitution needs reconciliation before a CRS field lands.
- Principle VI (format pluralism) conflicts with `.chil` as *canonical*, because promoting one format is what VI forbids. CHILmesh treats `.chil` as **one more readable/writable adapter**, not a privileged one. Registry-canonical is a Valence decision, not a CHILmesh one.
- The fort.13 nodal-attribute round-trip is a prerequisite for `kind="mesh"` losslessness. It is entirely unbuilt.

---

## 5. Recommendation

1. **Add a `.chil` export adapter first.** `write_chil` serializes only what CHILmesh already holds (points, connectivity, boundary_segments) under `kind="mesh"`. CRS is recorded as `native = "unknown"`/passthrough, with no transform, which respects Principle V. There is no Identity hash and no WGS84 transform. The adapter is purely additive and needs no constitution amendment.
2. **Defer the Identity/CRS/Domain-Boundary half to Valence.** The registry, hashing and curation already live there. #154 D6/D9/D15 put the heavy deps (shapely/pyproj/numpy) on the registry side. CHILmesh's minimal base install should not absorb them.
3. **Build the fort.13 reader/writer.** It is the highest-value independent prerequisite. It is useful on its own (ADCIRC nodal attributes) regardless of `.chil`.

Net: CHILmesh's I/O is **already 90% of the way** to being a `.chil` producer for the geometry payload. The missing 10% is the part that #154's adversarial review places upstream in the registry.

---

## 6. Frontmatter and split payload

The `.chil` manifest can act as frontmatter, so LLMs never read the bulk arrays.

> Operator ask (#201, 2026-06-14): *"use a context-management technique like
> frontmatter so we don't have LLM models reading and writing the mesh data
> (which can be huge) themselves, let the .py files do it."*

### Problem

A WNAT-scale mesh is ~10⁶–10⁷ nodes. A single `cat mesh.14` exceeds an LLM context window and costs tokens. It also risks corrupting binary data passed through a tool call. An agent should never load array bytes into context. Only Python should touch them.

### Pattern

Frontmatter is text and LLM-owned. The payload is binary and Python-owned.

| Layer | Content | Who reads/writes | Size |
|---|---|---|---|
| **Frontmatter** | `manifest.json`: counts (NP/NE), bbox, per-field stats (min/max/mean of depth, quality), layer summary (`n_layers`, per-layer counts), `origin`, sha256 of each payload, CRS | LLM + Python | ~1 KB |
| **Payload** | points / connectivity / nodal-attribute arrays (`.npy`/`.npz`, or referenced `fort.14`+`fort.13` by path+sha) | **Python only** | MB–GB |

The Fable review (comment 4680372228 §3) already recommends a "stdlib-queryable `manifest.json` as first archive entry". The new framing here is the *agent contract*. An agent inspecting a mesh reads **only** the frontmatter. Any operation on the arrays goes through a Python entry point.

### Already half-built

A `fort.14` header line (`NE NP`) is proto-frontmatter. `CHILmesh` already computes everything a manifest needs (`elem_quality`, bbox, `layers`, boundary counts).

One piece is missing: a single `summarize(mesh) -> dict` / `chilmesh summary <file>` that emits the frontmatter **without** the agent opening the raw file.

### Steps (additive, no constitution touch)

1. Add `chilmesh.summary(path_or_mesh) -> dict` and a `python -m chilmesh summary <file>` CLI. It parses header/metadata only (lazy). It never loads the full array block for counts that live in the header.
2. Emit that dict as the `.chil` `manifest.json` frontmatter when the `write_chil` export adapter (§5.1) lands.
3. Add an agent guardrail: a local development hook (not part of the package). It refuses a raw `cat`/`head`/`tail` of `*.14`/`*.grd`/`*.2dm`/`*.13`/`*.msh`/`*.npy`/`*.npz`/`fort.NNN` files over 64 KB. It points to `chilmesh summary <file>` instead.

Net: the format already points this way (manifest-first). Formalizing `summarize()` and the read-only-frontmatter agent contract is the cheap win. It is useful immediately for mesh inspection in routine sessions, independent of `.chil`.

---

## 7. Borrowing from `admesh/admesh`

CHILmesh can borrow one idea from `admesh/admesh` and no source code.

> Operator ask (#201, 2026-06-14): *"can we borrow any smart data structures,
> I/O, anything from https://github.com/admesh/admesh"*

`admesh/admesh` (Anthony Martin's ADMesh; **unrelated** to `domattioli/ADMESH`) is a C library for **repairing STL triangle-soup** for 3-D printing.

### License blocker

admesh/admesh is **GPL-2.0**. CHILmesh is **PolyForm Noncommercial 1.0.0**. The licenses are incompatible, so GPL code **cannot** be vendored or copied into CHILmesh. Ideas and algorithms are not copyrightable, so only those can be borrowed. **Never copy source.**

### Overlap

| admesh/admesh | CHILmesh equivalent | Verdict |
|---|---|---|
| `stl_hash_edge {key[6]; facet_number; which_edge; *next}`; hash table matching shared edges from coordinate keys, O(1) | **`EdgeMap`** (Phase-1, hash O(1) edge lookup) | **Already independently adopted.** Same idea; nothing to borrow. |
| `stl_neighbors {neighbor[3]; which_vertex_not[3]}`; per-facet compact face-adjacency | `Edge2Elem` + `Elem2Edge` | Equivalent info; admesh's encoding is tri-only/denser but CHILmesh is mixed tri/quad → not portable. |
| **Mesh repair**: unconnected-edge stitching (tolerance match) + hole-filling by facet insertion + degenerate-facet removal + normal repair | **None**; CHILmesh assumes watertight input (fort.14/ADMESH output is watertight by construction) | **The only additive idea.** Low value *today* (inputs are clean), but the *algorithm* (tolerance edge-stitch → hole-fill) is the recipe if CHILmesh ever ingests dirty soup (CAD/STL import, untrusted `.2dm`). |
| I/O: ASCII/binary **STL**, **OFF**, DXF, VRML export | fort.14, `.2dm`, gmsh `.msh`, fort.13 | Different domain; STL is **3-D single-precision surface** (`float x,y,z`); CHILmesh is **2-D double-precision** planar hydro mesh. STL is lossy (single-precision, no boundary types, no depth). `OFF` is the only cheap, plausibly-useful export bridge; STL/DXF/VRML are out of scope. |

### Mismatch

admesh is `float`-precision, 3-D and **tri-only**. It repairs surfaces and has no layers, skeletonization, quads or smoothing. CHILmesh is double-precision, 2-D and mixed-element, with topology, skeleton and FEM smoothing. The two solve different problems.

### Recommendation

Borrow **one idea** and copy no source. File a *future* `request: research` note for a tolerance-based **mesh-heal pass** (unconnected-edge stitch + hole-fill). Do this **iff** a dirty-input ingest path appears (STL/CAD or untrusted `.2dm`).

CHILmesh already has the hash-edge connectivity that admesh is known for (`EdgeMap`). The borrowed idea needs no code, no dependency and no license entanglement.

---

_Investigation for #201. References: #154 (`.chil` design + RESHAPE verdict),
constitution Principles V/VI/VII, `CHILmesh.py` I/O surface, `gmsh_io.py`
adapter precedent, `admesh/admesh` (GPL-2.0, STL repair)._
