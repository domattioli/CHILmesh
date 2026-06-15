# Handoff — mesh_segmenter sandbox session

For a fresh session that boots **inside the new `mesh_segmenter` repo** to prototype.
Everything here was designed in two grill-with-docs sessions on CHILmesh #153; this
folder (`docs/proposals/mesh-segmenter/`) is the seed — **lift it whole into the new
repo as `docs/`** (ADR → `docs/adr/`, the rest alongside).

## 0. Mission (one line)

`mesh_segmenter` = interactive, composable sub-region **selection** over a 2D mesh
("SAM2 for meshes"). **It selects; it never generates.**

## 1. Locked decisions (full rationale: `ADR-0001-architecture.md`)

- **Standalone repo**, provisional import name `mesh_segmenter` (a cooler name wanted
  before first publish — not urgent).
- **Engine = `chilmesh` ONLY.** Accept `admesh.Mesh` / fort.14 / raw `(nodes,
  elements)` → build a `CHILmesh` internally for adjacency, layers, spatial index,
  `submesh`. Because CHILmesh builds from raw arrays, the segmenter is
  mesh-library-agnostic.
- **Never import a generator** (`admesh` / `quadmesh`). A future "one-stop umbrella"
  (chilmesh + admesh + quadmesh + valence + mesh_segmenter) owns the
  `selection → re-mesh` wiring. Keeps the segmenter a leaf, avoids a dependency cycle.
- **Core object = `Selection`** — an immutable **element mask** bound to a parent
  mesh. Set-algebra `|` `&` `~` returns new Selections.
- **Exports are lazy/derived** — `element_ids`, `boundary` (raw polygon rings), an
  optional `submesh`. Never emit an `admesh.Domain`; the umbrella wraps rings.

## 2. v1 mechanism contracts (grill round 2)

- **Adjacency** = per-call `connectivity="edge"|"vertex"`, default `"edge"` (no
  corner-bleed at pinch points). Components + `Selection.boundary` perimeters are
  **edge-only** regardless of how the mask was grown.
- **`Selection.boundary`** = per-component list of `(outer_ring, holes[])`; never
  auto-cleaned. `Selection.components()` yields per-component sub-Selections. Pinch
  rings flagged (warning), never silently merged — a mask must not self-edit.
- **`by_click` criterion** = edge-crossing predicate `fn(from_elem, to_elem) -> bool`
  (`True` = stop). Named shortcuts (`"connected"`, `("field_jump", field, delta)`)
  wrap it. v2 learned model is just another crossing predicate.
- **Fields** = plain numpy arrays (`n_nodes` | `n_elements`, auto-detect, nodal
  reduces via "all verts in range" default). Sources: mesh-attached, chilmesh-computed
  (quality / edge-length / area), umbrella-supplied admesh size-field components.
  "All information" reaches the segmenter as an array — never an admesh import.

## 3. v1 surface to build

`Selection` + set-algebra; mechanisms `grow` (ring dilation), `by_distance` (shapely),
`by_click` (flood-fill + predicate), `by_threshold` (scalar + reduction), `by_polygon`,
`by_skeleton_layer` (exposes CHILmesh peels). Naming note: the issue's `by_layer` is
**`grow`**; "Layer" stays CHILmesh's skeletonization peel.

## 4. FIRST TASK this session — the prototype, not v1

Before committing to v2, run the **rasterize → SAM2 bootstrap prototype** in
`PROTOTYPE-rasterize-sam2.md`. It decides (no training) whether SAM2 beats the cheap
deterministic baseline. Build order:

1. **P0 (CPU, no model)** — `prototype/raster.py` + `project.py` + synthetic 2-basin
   `fixtures.py` + `test_roundtrip_recall`. Watershed stub = baseline IoU.
2. Eyeball numbers, then **P1** (real SAM2 via `huggingface_hub`), then **P2** gate.

Gate: `SAM2 IoU > baseline IoU` AND IoU ≥ 0.70 on ≥ 2 cases AND jitter IoU ≥ 0.80. Miss
any → drop SAM2, ship deterministic `by_click(field_gradient)`.

Keep prototype under `prototype/` with its own `[proto]` extra — throwaway, not the
shipped API, never in chilmesh.

## 5. Repo setup (cold start)

```bash
# in the new mesh_segmenter repo root
python -m venv .venv && . .venv/bin/activate
pip install -e ../CHILmesh          # engine, editable from sibling checkout
pip install -e ".[dev,proto]"       # numpy scipy scikit-image (+ torch sam2 hf for P1)
# package skeleton: mesh_segmenter/{__init__,selection,mechanisms/}.py
# prototype lives in prototype/ (see PROTOTYPE-rasterize-sam2.md file map)
```

CHILmesh entry points the engine gives you for free: `CHILmesh(connectivity=elements,
points=nodes, build_spatial_indices=True)`, `.elements_in_layer(i)`, `.submesh(ids)`,
`Edge2Elem` (−1 = boundary), `Vert2Elem`, fort.14 read/write, `from_admesh_domain`.

## 6. Hard truth to keep in view

The pipeline's real bottleneck is **NOT segmentation** — it's **conforming re-stitch**
(re-inserting a re-meshed sub-region into the parent mesh with matched boundary nodes,
no T-junctions). That lives in the **umbrella**, not here. The segmenter only owes a
clean boundary ring. Don't over-invest in segmentation polish until the re-stitch path
is proven viable elsewhere.

## 7. Open questions

- Cooler final package name (provisional `mesh_segmenter`).
- Real demand beyond the operator — coastal "select shallow elements" plausible but
  unvalidated. v1 is cheap enough that this doesn't gate building it.
- Whether the umbrella repo exists yet — segmenter's standalone value is limited until
  it does.

## 8. Refs

- CHILmesh #153 — issue + both grill comments (round 1 architecture, round 2
  contracts).
- `ADR-0001-architecture.md`, `CONTEXT.md`, `PROTOTYPE-rasterize-sam2.md` (this
  folder).
- ADMESH `docs/adr/ADR-001-chilmesh-boundary.md` — the consumer-side / sibling-package
  precedent (segmenter is chilmesh-engined; no back-ref added per operator).
- DomI #268 — skill-load recurrence (why a session may run a DomI skill via SKILL.md
  emulation instead of the registered Skill).

## 9. Conventions in the new repo

- Branch: work on `development`, draft PR `development → main` (mirror CHILmesh policy).
- Coding dispatch: code → Haiku subagent; main session plans/reviews/integrates.
- Caveman mode active for orchestrator/technical exchange.
