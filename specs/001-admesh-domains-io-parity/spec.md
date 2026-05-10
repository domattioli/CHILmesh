# Feature Specification: ADMESH-Domains I/O Parity

**Feature Branch**: `planning-optimize_modernize`
**Created**: 2026-04-26
**Status**: Draft
**Input**: CHILmesh must load any mesh from ADMESH-Domains catalog correctly, save it back without data loss, and expose metadata fields compatible with ADMESH-Domains schema — so CHILmesh can serve as computational engine behind catalog.

---

## User Scenarios & Testing *(mandatory)*

### User Story 1 — Load Real Coastal Mesh from Catalog (Priority: P1)

Researcher uses ADMESH-Domains to find coastal simulation mesh, downloads it, opens it in CHILmesh to analyse element quality and layer structure. Today this fails silently or crashes because most real-world meshes contain quad or mixed elements and CHILmesh only accepts triangles.

**Why this priority**: Blocking bug — catalog exists to deliver meshes to tools like CHILmesh; if loading fails, entire integration story breaks.

**Independent Test**: Load one quad-element and one mixed-element mesh from catalog. CHILmesh initialises without error, reports correct node and element counts, and analysis methods (`elem_quality`, `interior_angles`) return valid results.

**Acceptance Scenarios**:

1. **Given** valid ADCIRC `.fort.14` file containing quad elements, **When** researcher opens it in CHILmesh, **Then** mesh loads successfully with correct number of nodes, elements, and element type reported as "Quadrilateral" or "Mixed-Element".

2. **Given** valid SMS `.2dm` file from catalog, **When** researcher opens it in CHILmesh, **Then** mesh loads successfully and all geometry is preserved.

3. **Given** mesh loaded and analysed, **When** researcher saves it back to ADCIRC format, **Then** saved file round-trips identically (same node and element counts, same geometry).

---

### User Story 2 — Bulk Metadata Queries on Catalog (Priority: P2)

Developer queries ADMESH-Domains catalog to find all meshes with fewer than 50,000 nodes in specific region. Catalog's `node_count` and `bounding_box` fields need to be accurate. Developer also wants to load several meshes quickly to verify metadata — without triggering slow full initialisation.

**Why this priority**: Catalog usefulness depends on accurate, consistently populated metadata. Fast loading needed for bulk workflows.

**Independent Test**: Load ten catalog meshes with layers disabled. Each loads in under 2 seconds. Call `admesh_metadata()` on each; returned `node_count` and `element_count` match file's actual contents.

**Acceptance Scenarios**:

1. **Given** researcher who needs only metadata, **When** they load mesh with fast option, **Then** mesh initialises in under 2 seconds (vs. ~30s today for large meshes) and `node_count`, `element_count`, `element_type`, `bounding_box` all available.

2. **Given** mesh loaded via catalog record, **When** `admesh_metadata()` called, **Then** returned dict contains `node_count`, `element_count`, `element_type`, `bounding_box`, and each value matches file.

3. **Given** catalog `Mesh` record object, **When** developer calls `CHILmesh.from_admesh_domain(record)`, **Then** CHILmesh returned correctly without developer needing to know which file format or reader to use.

---

### User Story 3 — Contribute a New Domain with Validated Metadata (Priority: P3)

Hydrologist contributes new mesh to ADMESH-Domains. Before submitting, wants to verify metadata fields (`node_count`, `element_type`, `bounding_box`) match what CHILmesh actually computes from file.

**Why this priority**: Prevents catalog data drift. Lower priority because catalog can function without it, but quality degrades over time.

**Independent Test**: Load mesh and call `admesh_metadata()`. Compare each field against values contributor entered in manifest. Any mismatch surfaces clearly.

**Acceptance Scenarios**:

1. **Given** mesh file contributor is preparing to submit, **When** they load it in CHILmesh and call `admesh_metadata()`, **Then** returned dict contains all four schema-compatible fields with correct values they can copy into manifest.

2. **Given** mesh where contributor entered wrong `node_count`, **When** `admesh_metadata()` compared to manifest entry, **Then** mismatch immediately visible.

---

### Edge Cases

- File does not exist on disk → clear error with guidance to download from catalog first.
- File contains unsupported element type (e.g., 5-node) → raises descriptive error naming unsupported element type and element ID.
- Mesh loaded without layers, then `get_layer()` called → raises descriptive error directing user to enable layers on reload.
- `from_admesh_domain()` called with `Mesh` record whose `type` field is unrecognised → falls back to ADCIRC reader and logs warning.

---

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST load ADCIRC `.fort.14` meshes containing triangular, quad, and mixed elements without raising error.

- **FR-002**: System MUST save any loaded mesh (including quad and mixed-element) back to ADCIRC `.fort.14` format with correct element node counts, preserving round-trip fidelity.

- **FR-003**: System MUST load SMS `.2dm` format meshes containing triangular and/or quad elements.

- **FR-004**: System MUST expose metadata method returning at minimum: node count, element count, element type classification, and axis-aligned bounding box — all consistent with ADMESH-Domains Mesh schema fields.

- **FR-005**: System MUST support initialisation without computing skeletonisation layers, completing in under 2 seconds for meshes up to ~15,000 elements.

- **FR-006**: System MUST provide single entry point accepting ADMESH-Domains `Mesh` record (or equivalent dict) and returning correctly initialised mesh object, routing to right reader based on record's `type` field.

- **FR-007**: System MUST raise descriptive, actionable error when mesh file is not found, when unsupported element type encountered, or when layer-dependent methods called on mesh initialised without layers.

### Key Entities

- **Mesh record**: ADMESH-Domains catalog entry with at minimum filename, format type, and optionally metadata fields (node_count, element_count, element_type, bounding_box, kind).
- **Mesh object**: In-memory representation of loaded mesh, with geometry, topology, and optionally skeletonisation layers.
- **Metadata dict**: Format-agnostic dict of computed mesh properties compatible with ADMESH-Domains schema field names.

---

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: All mesh files in ADMESH-Domains catalog marked as `test_case: true` load without error.

- **SC-002**: Round-trip fidelity: mesh loaded then saved produces file with identical node count and element count to original.

- **SC-003**: Fast loading completes in under 2 seconds for mesh of up to ~15,000 elements when layers not requested (currently ~30 seconds).

- **SC-004**: `admesh_metadata()` values match file contents to 100% accuracy for node count and element count across all test fixtures.

- **SC-005**: Zero regressions — all tests that pass today continue to pass after this feature delivered.

---

## Assumptions

- ADMESH-Domains `Mesh` records are duck-typed; no hard import of `admesh-domains` package required inside CHILmesh.
- ADCIRC_GRD files use same line format as standard ADCIRC `.fort.14` and can be read by same reader.
- Element types beyond triangles (3 nodes) and quads (4 nodes) are out of scope; 5-node or higher elements raise error.
- Coordinate values in mesh files treated as unitless; CHILmesh does not interpret whether they represent geographic degrees or projected metres.
- Boundary meshes (`kind: "boundary"`) do not require skeletonisation layers; computing layers on boundary mesh issues warning but does not error.
- Existing 4 triangular test fixtures (annulus, donut, block_o, structured) remain primary regression suite; new minimal quad fixture added for FR-001/FR-002 coverage.
