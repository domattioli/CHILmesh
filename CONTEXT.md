# CHILmesh Context

**Purpose:** Vocabulary for CHILmesh runtime terms and cross-repo coordination with ADMESH-Domains (canonical domain/mesh registry).

## CHILmesh-local terms

- **CHILmesh** — the runtime Python class `src/chilmesh/CHILmesh.py`; represents an instantiated mesh object
- **Vertex** — a mesh point (x, y); used in this codebase, NOT "node" (fort.14 I/O files use "node" but internal code uses "vertex")
- **Element** — a triangular or quadrilateral cell; stored in `Elem2Vert` adjacency
- **ElementType** — `"tri"`, `"quad"`, or `"mixed"`; stored in `CHILmesh.kind` (distinct from ADMESH-Domains `kind` — see Flagged Ambiguities)
- **Edge** — a pair of vertices (v_i, v_j) with v_i < v_j (canonical ordering); stored in `Edge2Vert`
- **EdgeMap** — hash-map from sorted vertex pair → edge ID; core O(1) adjacency structure from Phase 1 optimization
- **Layer** — one "peel" in the medial-axis skeletonization; each layer stores OE (outer edges), IE (inner edges), OV (outer vertices), IV (inner vertices)
- **Skeletonization** — the iterative layer-peeling algorithm yielding `CHILmesh.layers`; controlled by `compute_layers=True` at init
- **compute_layers** — `CHILmesh.__init__` kwarg; `True` = full init with skeletonization (~5× slower), `False` = fast adjacency-only init
- **Reader / Writer** — fort.14 (`read_fort14`, `write_fort14`) and 2dm I/O; fort.14 is primary, 2dm read-only
- **from_admesh_domain** — bridge constructor: `CHILmesh.from_admesh_domain(domain)` builds a CHILmesh object from an ADMESH-Domains `Domain`/`Mesh` record

## Upstream terms (ADMESH-Domains canonical)

ADMESH-Domains owns the canonical definitions of: **Domain**, **Mesh** (registry record), **Mesh Record**, **Domain Boundary**, **Mesh Boundary**, **Identity**, **Relatedness**, **Containment**, **content_uid**, **kind** (registry sense), **Curator**, **Submitter**. Cite the upstream: `domattioli/ADMESH-Domains → CONTEXT.md`.

## Flagged ambiguities

- **`Mesh` vs `CHILmesh` vs `Mesh Record`** — "Mesh" is overloaded: ADMESH-Domains uses it for the catalog entry; CHILmesh uses it colloquially for the runtime object. Prefer `CHILmesh` (class name) in this codebase; `Mesh Record` for the catalog entry.
- **`kind` overload** — CHILmesh uses `kind ∈ {"tri", "quad", "mixed"}` for element topology; ADMESH-Domains uses `kind ∈ {"mesh", "boundary"}` for record type. Different meanings; do not conflate.
- **`node` vs `Vertex`** — fort.14 format and ADCIRC community say "node"; CHILmesh internal code uses "vertex". Keep the distinction: I/O layer can emit "node" in docs; compute layer always says "vertex".

## Pending: .chil format design (CHILmesh#154)

A `.chil` file format was designed in a Socratic session and an adversarial review returned verdict "RESHAPE". The design affects CHILmesh in three ways:

- **Constitution V (CRS Agnosticism)**: `.chil` Identity forcing WGS84 conflicts with CHILmesh's coordinate-opaque stance. Resolution pending.
- **Constitution VI (Format Pluralism)**: `.chil` as the canonical registry format could privilege one format over fort.14/2dm. Resolution: `.chil` should be one of many readable formats but the format the registry *writes*; CHILmesh remains pluralist at the library layer.
- **Constitution VII (API Stability)**: `CHILmesh.kind` rename implications if `.chil` uses a different `kind` vocabulary. Requires a deprecation cycle.

Until #154 is resolved and constitutions amended, `.chil` I/O is NOT in scope for CHILmesh.

## Cross-references

- CHILmesh#154 (`.chil` format design + adversarial review)
- CHILmesh#149 (prior CONTEXT.md draft, closed)
- ADMESH-Domains#79 (companion CONTEXT.md + ADR-0001 draft)
- ADMESH-Domains#80 (Domain lineage graph)
- `.specify/memory/constitution.md` (CHILmesh governance contract)
