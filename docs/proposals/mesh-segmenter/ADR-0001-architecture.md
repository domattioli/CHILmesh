# ADR-0001 — mesh-segmenter: standalone package, chilmesh engine, neutral mask

| Field | Value |
|---|---|
| Status | Proposed |
| Date | 2026-06-13 |
| Origin | [CHILmesh #153](https://github.com/domattioli/CHILmesh/issues/153) (filed as the `admesh-segmenter` proposal; ADMESH-side as #9) |
| Method | grill-with-docs design session |
| Related | ADMESH [`docs/adr/ADR-001-chilmesh-boundary.md`](https://github.com/domattioli/ADMESH/blob/main/docs/adr/ADR-001-chilmesh-boundary.md) |

## Context

The proposal is a composable sub-region **selection** API over a 2D mesh ("SAM2 for
meshes"): pick elements by layer / distance / click / threshold / polygon, combine
with set-algebra, usually to identify a sub-region to **re-mesh**. The issue body
assumed it would live in ADMESH and consume `admesh.Mesh`; a later comment proposed
"a package within admesh"; another raised a SAM2-style *learned* segmentation.

Two facts from the codebase contradict the admesh framing:

1. **`admesh.Mesh` is adjacency-free** (`frozen`, `slots`; `nodes / elements /
   boundaries / bathymetry / quality / title`). Every adjacency-walk mechanism the
   proposal needs (ring expansion, flood-fill, connected components) requires
   topology ADMESH deliberately does not keep.
2. **CHILmesh already ships that topology and ~60–70% of the proposed API** —
   `Edge2Elem`/`Vert2Elem`, `elements_in_layer`, `submesh` (renumbering),
   `_skeletonize`, `build_spatial_indices`, fort.14 I/O — all built from raw
   `(nodes, elements)`.

ADMESH's ADR-001 had already classified segmentation as *consumer-side* and
reaffirmed "heavy consumer-side functionality earns its own package." This ADR
finishes that thread by naming the package, its engine, and its dependency rules.

## Decision

**A standalone `mesh-segmenter` package whose only topology dependency is
`chilmesh`. It selects; it never generates.**

1. **Placement — standalone repo `mesh-segmenter`.** Not an ADMESH submodule
   (would bloat the generator import surface and contradict ADR-001) and not
   in-tree CHILmesh (keeps shapely / future ML deps out of CHILmesh core).

2. **Engine — chilmesh only.** Accepts `admesh.Mesh` / fort.14 / raw
   `(nodes, elements)` and builds a `CHILmesh` internally for adjacency, layers,
   spatial index, and `submesh`. Because CHILmesh constructs from raw arrays, the
   segmenter is mesh-library-agnostic — justifying the generic name.

3. **No generator dependency.** mesh-segmenter MUST NOT import `admesh` or
   `quadmesh`. A future "one-stop" **umbrella** package (depends on chilmesh,
   admesh, quadmesh, valence, *and* mesh-segmenter) owns the `selection → re-mesh`
   wiring. This keeps the segmenter a leaf and prevents a dependency cycle the
   moment the umbrella composes them.

4. **Core object — `Selection`, a canonical immutable element mask** bound to a
   parent mesh. Set-algebra `|` `&` `~` returns new Selections. Nodal signals
   (bathymetry, curvature) reduce onto elements via a documented rule (default
   conservative "all vertices in range"; opt-in `any` / `mean`).

5. **Output — neutral artifact at the re-mesh seam.** A Selection emits
   `element_ids`, `boundary` (raw polygon rings, holes supported), and an optional
   `submesh`; never an `admesh.Domain` directly. The umbrella wraps rings → Domain
   for re-gen, or feeds the submesh → quadmesh tri2quad.

6. **Scope — two phases.**
   - **v1 (deterministic):** `Selection` + set-algebra; `grow` (ring dilation from a
     seed), `by_distance` (shapely), `by_click` (flood-fill w/ predicate),
     `by_threshold` (scalar + nodal reduction), `by_polygon`, `by_skeleton_layer`.
   - **v2 (research spike):** SAM2-*inspired* learned mask-from-size-field
     ("click in the gulf → the gulf, respecting bathymetry"). Transfer-learning,
     **not** train-from-scratch; uncertain, explicitly fenced out of v1. Plugs into
     the same `Selection`. De-risking bridge: v1
     `by_click(criterion=size_field_gradient)` is already region-growing bounded by
     the size field — the deterministic MVP of the SAM2 idea.

## Consequences

**Enables**
- A buildable v1 with no ML risk, immediately useful for the re-mesh pipeline.
- Clean layering: umbrella → {generators, segmenter}; segmenter → chilmesh.
- v2 ML work slots in behind a stable `Selection` contract — no API churn.

**Costs / forecloses**
- A naming correction: the issue's `by_layer` becomes **`grow`**; "Layer" stays
  reserved for CHILmesh's skeletonization peel (`by_skeleton_layer`). See
  [`CONTEXT.md`](./CONTEXT.md).
- The umbrella must own selection→generator glue; the segmenter cannot offer a
  one-call `.re_mesh()` convenience without breaking the leaf rule.
- ADMESH ADR-001 should gain a back-reference: the "segmenter sibling" it
  anticipated is engined by **chilmesh**, not admesh.

## Follow-up

- [ ] Operator: create the `mesh-segmenter` repo; lift this folder in as `docs/adr/`
      + `CONTEXT.md`.
- [ ] Confirm package/import name (`mesh_segmenter`? `meshseg`?).
- [ ] Note in ADMESH ADR-001 that the anticipated sibling is chilmesh-engined.
- [ ] v2: separate research spec before any training/transfer-learning work.
