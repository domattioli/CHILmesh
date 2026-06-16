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

_Investigation for #201. References: #154 (`.chil` design + RESHAPE verdict),
constitution Principles V/VI/VII, `CHILmesh.py` I/O surface, `gmsh_io.py`
adapter precedent._
