# Feature Specification: ADMESH-Domains I/O Parity

**Feature Branch**: `planning-optimize_modernize`
**Created**: 2026-04-26
**Status**: Draft
**Input**: CHILmesh must load any mesh from the ADMESH-Domains catalog correctly,
save it back without data loss, and expose metadata fields compatible with the
ADMESH-Domains schema — so that CHILmesh can serve as the computational engine
behind the catalog.

---

## User Scenarios & Testing *(mandatory)*

### User Story 1 — Load Real Coastal Mesh from Catalog (Priority: P1)

A researcher uses ADMESH-Domains to find a coastal simulation mesh, downloads
it, and opens it in CHILmesh to analyse element quality and layer structure.
Today this fails silently or crashes because most real-world meshes contain
quad or mixed elements and CHILmesh only accepts triangles.

**Why this priority**: Blocking bug — the catalog exists to deliver meshes to
tools like CHILmesh; if loading fails, the entire integration story breaks.

**Independent Test**: Load one quad-element and one mixed-element mesh from the
catalog. CHILmesh initialises without error, reports correct node and element
counts, and analysis methods (`elem_quality`, `interior_angles`) return valid
results.

**Acceptance Scenarios**:

1. **Given** a valid ADCIRC `.fort.14` file containing quad elements,
   **When** the researcher opens it in CHILmesh,
   **Then** the mesh loads successfully with the correct number of nodes,
   elements, and element type reported as "Quadrilateral" or "Mixed-Element".

2. **Given** a valid SMS `.2dm` file from the catalog,
   **When** the researcher opens it in CHILmesh,
   **Then** the mesh loads successfully and all geometry is preserved.

3. **Given** a mesh that was loaded and analysed,
   **When** the researcher saves it back to ADCIRC format,
   **Then** the saved file round-trips identically (same node and element
   counts, same geometry).

---

### User Story 2 — Bulk Metadata Queries on Catalog (Priority: P2)

A developer queries the ADMESH-Domains catalog to find all meshes with fewer
than 50,000 nodes in a specific region. The catalog's `node_count` and
`bounding_box` fields need to be accurate. The developer also wants to load
several meshes quickly just to verify metadata — without triggering a slow
full initialisation on each one.

**Why this priority**: Catalog usefulness depends on accurate, consistently
populated metadata. Fast loading is needed for bulk workflows.

**Independent Test**: Load ten catalog meshes with layers disabled. Each loads
in under 2 seconds. Call `admesh_metadata()` on each; the returned
`node_count` and `element_count` match the file's actual contents.

**Acceptance Scenarios**:

1. **Given** a researcher who needs only metadata,
   **When** they load a mesh with the fast option,
   **Then** the mesh initialises in under 2 seconds (vs. ~30s today for
   large meshes) and `node_count`, `element_count`, `element_type`, and
   `bounding_box` are all available.

2. **Given** a mesh loaded via the catalog record,
   **When** `admesh_metadata()` is called,
   **Then** the returned dict contains `node_count`, `element_count`,
   `element_type`, and `bounding_box`, and each value matches the file.

3. **Given** a catalog `Mesh` record object,
   **When** the developer calls `CHILmesh.from_admesh_domain(record)`,
   **Then** CHILmesh is returned correctly without the developer needing to
   know which file format or reader to use.

---

### User Story 3 — Contribute a New Domain with Validated Metadata (Priority: P3)

A hydrologist contributes a new mesh to ADMESH-Domains. Before submitting,
they want to verify that the metadata fields they plan to enter (`node_count`,
`element_type`, `bounding_box`) match what CHILmesh actually computes from
the file.

**Why this priority**: Prevents catalog data drift. Lower priority because the
catalog can function without it, but quality degrades over time without it.

**Independent Test**: Load a mesh and call `admesh_metadata()`. Compare each
field against the values the contributor entered in the manifest. Any mismatch
is surfaced clearly.

**Acceptance Scenarios**:

1. **Given** a mesh file the contributor is preparing to submit,
   **When** they load it in CHILmesh and call `admesh_metadata()`,
   **Then** the returned dict contains all four schema-compatible fields with
   correct values they can copy into the manifest.

2. **Given** a mesh where the contributor entered a wrong `node_count`,
   **When** `admesh_metadata()` is compared to the manifest entry,
   **Then** the mismatch is immediately visible (correct value from CHILmesh
   vs. incorrect value in manifest).

---

### Edge Cases

- What happens when the mesh file does not exist on disk?
  → Clear error with guidance to download from catalog first.
- What happens when the file contains an unsupported element type (e.g., 5-node)?
  → Raises a descriptive error naming the unsupported element type and the element ID.
- What happens when a mesh is loaded without layers, then `get_layer()` is called?
  → Raises a descriptive error directing the user to enable layers on reload.
- What happens when `from_admesh_domain()` is called with a `Mesh` record
  whose `type` field is unrecognised?
  → Falls back to ADCIRC reader and logs a warning.

---

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The system MUST load ADCIRC `.fort.14` meshes containing
  triangular, quad, and mixed elements without raising an error.

- **FR-002**: The system MUST save any loaded mesh (including quad and
  mixed-element) back to ADCIRC `.fort.14` format with correct element node
  counts, preserving round-trip fidelity.

- **FR-003**: The system MUST load SMS `.2dm` format meshes containing
  triangular and/or quad elements.

- **FR-004**: The system MUST expose a metadata method returning at minimum:
  node count, element count, element type classification, and axis-aligned
  bounding box — all consistent with the ADMESH-Domains Mesh schema fields.

- **FR-005**: The system MUST support initialisation without computing
  skeletonisation layers, completing in under 2 seconds for meshes up to
  ~15,000 elements.

- **FR-006**: The system MUST provide a single entry point that accepts an
  ADMESH-Domains `Mesh` record (or equivalent dict) and returns a correctly
  initialised mesh object, routing to the right reader based on the record's
  `type` field.

- **FR-007**: The system MUST raise a descriptive, actionable error when a
  mesh file is not found, when an unsupported element type is encountered,
  or when layer-dependent methods are called on a mesh initialised without
  layers.

### Key Entities

- **Mesh record**: An ADMESH-Domains catalog entry with at minimum a filename,
  a format type, and optionally metadata fields (node_count, element_count,
  element_type, bounding_box, kind).
- **Mesh object**: The in-memory representation of a loaded mesh, with
  geometry, topology, and optionally skeletonisation layers.
- **Metadata dict**: A format-agnostic dict of computed mesh properties
  compatible with ADMESH-Domains schema field names.

---

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: All mesh files in the ADMESH-Domains catalog marked as
  `test_case: true` load without error.

- **SC-002**: Round-trip fidelity: a mesh loaded then saved produces a file
  with identical node count and element count to the original.

- **SC-003**: Fast loading completes in under 2 seconds for a mesh of up to
  ~15,000 elements when layers are not requested (currently ~30 seconds).

- **SC-004**: `admesh_metadata()` values match file contents to 100% accuracy
  for node count and element count across all test fixtures.

- **SC-005**: Zero regressions — all tests that pass today continue to pass
  after this feature is delivered.

---

## Assumptions

- ADMESH-Domains `Mesh` records are duck-typed; no hard import of the
  `admesh-domains` package is required inside CHILmesh.
- ADCIRC_GRD files use the same line format as standard ADCIRC `.fort.14`
  and can be read by the same reader.
- Element types beyond triangles (3 nodes) and quads (4 nodes) are out of
  scope; 5-node or higher elements raise an error.
- Coordinate values in mesh files are treated as unitless; CHILmesh does not
  interpret whether they represent geographic degrees or projected metres.
- Boundary meshes (`kind: "boundary"`) do not require skeletonisation layers;
  computing layers on a boundary mesh issues a warning but does not error.
- The existing 4 triangular test fixtures (annulus, donut, block_o, structured)
  remain the primary regression suite; a new minimal quad fixture will be
  added for FR-001/FR-002 coverage.
