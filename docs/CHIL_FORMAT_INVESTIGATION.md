# `.chil` Format — Codebase Integration Investigation (#201)

Investigation of where the proposed `.chil` mesh container (designed in #154,
verdict **RESHAPE**) would slot into CHILmesh's existing I/O surface, what is
already built, and what the object model is still missing.

> **Scope (per #201):** brainstorm only. Authoring/saving `.chil` files into the
> repo beyond examples is explicitly **out of scope**. This doc maps the
> integration surface and the gaps; it does not implement a reader/writer.

---

## 1. Current I/O surface (what's already built)

| Format | Read | Write | Entry points |
|---|---|---|---|
| ADCIRC `fort.14` (`.14`/`.fort14`) | ✅ | ✅ | `CHILmesh.read_from_fort14`, `write_to_fort14`, module `write_fort14` |
| SMS `.2dm` | ✅ | ✅ (private `_write_2dm`) | `CHILmesh.read_from_2dm`, `_write_2dm` |
| Gmsh `.msh` (v2.2 + v4.1) | ✅ | ✅ | `gmsh_io.read_msh`/`write_msh`, `CHILmesh.read_from_msh`/`write_to_msh` |
| Valence registry record | ✅ | — | `CHILmesh.from_admesh_domain(record)` |
| `.chil` | ❌ | ❌ | **none — does not exist anywhere in `src/`** |

**Unified dispatch already exists.** `CHILmesh.save(filename)` and
`CHILmesh.load(filename)` branch on file suffix (`CHILmesh.py:2346`,
`:2354`). This is the natural seam: a `.chil` adapter is a third suffix branch,
not a new public API shape. The adapter pattern is precedent —
`gmsh_io.py` is a standalone module the class delegates to.

**Note on Principle VI doc-drift:** constitution Principle VI says 2dm
"write support future", but `_write_2dm` already ships (private). Either
promote it to public `write_to_2dm` or correct the constitution text — flagged
separately, not part of this investigation.

---

## 2. Object model the I/O must serialize

A runtime `CHILmesh` carries (`CHILmesh.py:274–288`):

- `points` — `(N,3)` float (z padded to 0 for 2D)
- `connectivity_list` — `(M,3|4)` int, mixed tri/quad (quads padded)
- `boundary_segments` — `list[dict]` with ADCIRC IBTYPE / kind metadata
- `boundary_condition`, `seed_ibtypes`, `seed_boundary_kinds`
- `adjacencies` — derived (Elem2Vert, Edge2Vert, Edge2Elem, …) — **not serialized; recomputed on load**
- `layers` — skeletonization output (OE/IE/OV/IV/bEdgeIDs) — **derived, not serialized**

The minimalism rule in #154 (store nothing derivable) matches what the existing
formats already do: only points + connectivity + boundary records persist;
adjacencies and layers recompute. So a `.chil` writer's required payload is the
same three things the fort.14 writer already emits, plus the new metadata below.

---

## 3. Gaps — what the object model lacks for `.chil` v1

`.chil` (#154) requires fields CHILmesh does **not** currently model:

| `.chil` requirement (#154) | Present in CHILmesh? | Gap |
|---|---|---|
| Native CRS declaration (Q8) | ❌ | No CRS field; Principle V makes coords **opaque** → direct conflict (#154 D9 HIGH) |
| Domain Boundary distinct from Mesh Boundary (Q11) | ❌ | Only `boundary_segments` (mesh edges) exist; no domain-outline concept |
| `fort.13` nodal attributes for lossless mesh round-trip (Q11) | ❌ | No fort.13 read/write anywhere → "lossless fort.14+fort.13" unmet |
| Deterministic content hash / `content_uid` (Q15) | ❌ | No hashing in CHILmesh (lives in Valence schema, not here) |
| Multi-ring boundary (holes, Bermuda) (Q10/Q11) | partial | boundary_segments can list multiple rings, but no canonical winding/start-vertex normalization |
| Quantization to fixed grid (Q9) | ❌ | No coordinate quantization step |

**Conclusion:** of the three things a `.chil` writer needs, two
(points, connectivity+boundary) are already produced by the fort.14 path; the
**new** work is metadata the object does not hold (CRS, domain boundary, hash),
not mesh geometry.

---

## 4. Where `.chil` would integrate (the answer to #201's question)

If/when the RESHAPE blockers clear, integration is **localized**:

1. **New module `src/chilmesh/chil_io.py`** — mirror `gmsh_io.py`: `read_chil()`
   / `write_chil()` returning/consuming a `CHILmesh`. Container per #154 Q12 is a
   zip of `manifest.toml` + `vertices.npy` + `triangles.npy`/`quads.npy` +
   `boundaries.toml` + `provenance.toml`.
2. **Two suffix branches** in `CHILmesh.save`/`load` (`CHILmesh.py:2346`/`2354`)
   — `.chil` → `chil_io.write_chil`/`read_chil`. Zero churn to other call sites.
3. **`from_admesh_domain` becomes the primary producer.** Per #154 Q14 `.chil`
   is the registry's canonical store; the registry-bridge classmethod is where a
   CHILmesh first meets a `.chil` payload, so wire reading there too.
4. **Export adapters out of `.chil`** — fort.14 (lossless), 2dm (lossy, drops
   typed boundaries), msh — all already exist; they become the "write OUT of
   `.chil`" adapters #154 Q14 describes with no new code.

**What blocks shipping it (from #154 RESHAPE, CHILmesh-side):**
- Principle V (coord-agnostic) vs WGS84-canonical Identity — needs constitution
  reconciliation before a CRS field lands.
- Principle VI (format pluralism) vs `.chil` as *canonical* — promoting one
  format is exactly what VI forbids; CHILmesh should treat `.chil` as **one more
  readable/writable adapter**, not privileged. Registry-canonical is an
  Valence decision, not a CHILmesh one.
- fort.13 nodal-attribute round-trip is a prerequisite for `kind="mesh"`
  losslessness and is entirely unbuilt.

---

## 5. Recommendation

1. **Smallest first step that is constitution-safe:** add a `.chil` **export
   adapter** (`write_chil`) that serializes only what CHILmesh already holds
   (points, connectivity, boundary_segments) under `kind="mesh"`, with CRS
   recorded as `native = "unknown"`/passthrough (no transform — respects
   Principle V). No Identity hash, no WGS84 transform. This is a pure additive
   adapter, no constitution amendment required.
2. **Defer the Identity/CRS/Domain-Boundary half** to Valence where the
   registry, hashing, and curation already live (#154 D6/D9/D15 put the heavy
   deps — shapely/pyproj/numpy — on the registry side, which CHILmesh's minimal
   base install should not absorb).
3. **fort.13 reader/writer** is the highest-value independent prerequisite and
   is useful on its own (ADCIRC nodal attributes) regardless of `.chil`.

Net: CHILmesh's I/O is **already 90% of the way** to being a `.chil` producer
for the geometry payload; the missing 10% is deliberately the part #154's
adversarial review says belongs upstream in the registry, not here.

---

## 6. Frontmatter / split-payload so LLMs never read the bulk arrays

> Operator ask (#201, 2026-06-14): *"use a context-management technique like
> frontmatter so we don't have LLM models reading and writing the mesh data
> (which can be huge) themselves — let the .py files do it."*

**Problem.** A WNAT-scale mesh is ~10⁶–10⁷ nodes. A single `cat mesh.14` blows
an LLM context window, costs tokens, and risks the MCP base64-corruption hazard
already logged in CLAUDE.md. An agent should never load array bytes into
context; only Python should touch them.

**Pattern — frontmatter (text, LLM-owned) + payload (binary, Python-owned).**

| Layer | Content | Who reads/writes | Size |
|---|---|---|---|
| **Frontmatter** | `manifest.toml`: counts (NP/NE), bbox, per-field stats (min/max/mean of depth, quality), layer summary (`n_layers`, per-layer counts), `provenance`/`lineage`, sha256 of each payload, CRS | LLM + Python | ~1 KB |
| **Payload** | points / connectivity / nodal-attribute arrays (`.npy`/`.npz`, or referenced `fort.14`+`fort.13` by path+sha) | **Python only** | MB–GB |

This is the same **"stdlib-queryable `manifest.toml` as first archive entry"**
already recommended in the Fable review (comment 4680372228 §3). The new framing
is the *agent contract*, not just the format: an agent inspecting a mesh reads
**only** the frontmatter; any operation on the arrays goes through a Python
entry point.

**Already half-built.** A `fort.14` header line (`NE NP`) is proto-frontmatter,
and `CHILmesh` already computes everything a manifest needs (`elem_quality`,
bbox, `layers`, boundary counts). The missing piece is a single
`summarize(mesh) -> dict` / `chilmesh summary <file>` that emits the
frontmatter **without** the agent ever opening the raw file.

**Concrete steps (additive, no constitution touch):**
1. `chilmesh.summary(path_or_mesh) -> dict` + `python -m chilmesh summary <file>`
   CLI — parses header/metadata only (lazy; never loads the full array block for
   counts that live in the header).
2. Emit that dict as the `.chil` `manifest.toml` frontmatter when the
   `write_chil` export adapter (§5.1) lands.
3. **Agent guardrail** — ✅ shipped as the `scripts/hooks/mesh_read_guard.sh`
   PreToolUse hook (wired into `.claude/settings.json` for the `Bash` and `Read`
   matchers, same class as `branch_guard`/`secret_path_guard`). Refuses a raw
   `Read`/`cat`/`head`/`tail`/… of `*.14`/`*.grd`/`*.2dm`/`*.13`/`*.msh`/`*.npy`/
   `*.npz`/`fort.NNN` over `CHILMESH_MESH_READ_MAX_KB` (default 64 KB) and reroutes
   to `chilmesh summary <file>`; bypass with `CHILMESH_MESH_READ_GUARD_BYPASS=1`
   (logged to `~/.claude/hook-bypass.log`). Mirrors the `mcp-binary-push`
   refuse-and-reroute pattern in DomI. Smoke test:
   `scripts/hooks/tests/mesh_read_guard_smoke.sh` (11 scenarios).

**Net:** the format already points this way (manifest-first); formalizing
`summarize()` + the read-only-frontmatter agent contract is the cheap win and is
useful immediately (mesh inspection in routine sessions) independent of `.chil`.

---

## 7. Can we borrow from `admesh/admesh` (the STL library)?

> Operator ask (#201, 2026-06-14): *"can we borrow any smart data structures,
> I/O, anything from https://github.com/admesh/admesh"*

`admesh/admesh` (Anthony Martin's ADMesh — **unrelated** to `domattioli/ADMESH`)
is a C library for **repairing STL triangle-soup** for 3-D printing.

**License blocker first.** admesh/admesh is **GPL-2.0**; CHILmesh is **PolyForm
Noncommercial 1.0.0**. GPL code **cannot** be vendored/copied into CHILmesh —
incompatible licenses. So this is *borrow ideas/algorithms* (not copyrightable)
only — **never copy source**.

**What it has, and the overlap:**

| admesh/admesh | CHILmesh equivalent | Verdict |
|---|---|---|
| `stl_hash_edge {key[6]; facet_number; which_edge; *next}` — hash table matching shared edges from coordinate keys, O(1) | **`EdgeMap`** (Phase-1, hash O(1) edge lookup) | **Already independently adopted.** Same idea; nothing to borrow. |
| `stl_neighbors {neighbor[3]; which_vertex_not[3]}` — per-facet compact face-adjacency | `Edge2Elem` + `Elem2Edge` | Equivalent info; admesh's encoding is tri-only/denser but CHILmesh is mixed tri/quad → not portable. |
| **Mesh repair**: unconnected-edge stitching (tolerance match) + hole-filling by facet insertion + degenerate-facet removal + normal repair | **None** — CHILmesh assumes watertight input (fort.14/ADMESH output is watertight by construction) | **The only genuinely additive idea.** Low value *today* (inputs are clean), but the *algorithm* (tolerance edge-stitch → hole-fill) is the recipe if CHILmesh ever ingests dirty soup (CAD/STL import, untrusted `.2dm`). |
| I/O: ASCII/binary **STL**, **OFF**, DXF, VRML export | fort.14, `.2dm`, gmsh `.msh`, fort.13 | Different domain — STL is **3-D single-precision surface** (`float x,y,z`); CHILmesh is **2-D double-precision** planar hydro mesh. STL is lossy (single-precision, no boundary types, no depth). `OFF` is the only cheap, plausibly-useful export bridge; STL/DXF/VRML are out of scope. |

**Mismatch summary:** admesh is `float`-precision, 3-D, **tri-only** surface
repair with no layers/skeletonization/quad/smoothing. CHILmesh is double-precision
2-D mixed-element with topology + skeleton + FEM smoothing. Different problem.

**Recommendation:** borrow **one idea, not code** — file a *future* `request:
research` note for a tolerance-based **mesh-heal pass** (unconnected-edge stitch +
hole-fill) **iff** a dirty-input ingest path appears (STL/CAD or untrusted
`.2dm`). The hash-edge connectivity admesh is famous for, CHILmesh already has
(`EdgeMap`). No code, no dependency, no license entanglement.

---

_Investigation for #201. References: #154 (`.chil` design + RESHAPE verdict),
constitution Principles V/VI/VII, `CHILmesh.py` I/O surface, `gmsh_io.py`
adapter precedent, `admesh/admesh` (GPL-2.0, STL repair), CLAUDE.md
token-hygiene + `mcp-binary-push` reroute precedent._
