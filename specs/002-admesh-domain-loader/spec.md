# Feature Specification: ADMESH-Domains Loader

**Feature Branch**: `claude/peaceful-goldberg-CCD9l`  
**Created**: 2026-04-26  
**Status**: Draft  
**Input**: Create `from_admesh_domain()` classmethod to canonicalize CHILmesh loading from ADMESH-Domains Mesh records.

## User Scenarios & Testing

### User Story 1 - Load CHILmesh from ADMESH-Domains Mesh instance (Priority: P1)

A developer has an `admesh_domains.Mesh` object (a structured record with filename, type, and kind fields) and wants to load it as a CHILmesh with skeletonization enabled by default.

**Why this priority**: This is the primary use case for ADMESH-Domains integration; developers will call `from_admesh_domain(mesh_record)` as the canonical entry point.

**Independent Test**: Can load an ADCIRC-format Mesh instance, get a populated CHILmesh with adjacencies and layers computed.

**Acceptance Scenarios**:

1. **Given** an `admesh_domains.Mesh` instance with `filename="path/to/mesh.fort.14"`, `type="ADCIRC"`, `kind="mesh"`, **When** `CHILmesh.from_admesh_domain(mesh_record)` is called, **Then** a CHILmesh is returned with `_is_boundary=False`, adjacencies built, and skeletonization layers computed.
2. **Given** an `admesh_domains.Mesh` with `kind="boundary"`, **When** called with `from_admesh_domain(mesh_record)`, **Then** layers are not computed and `_is_boundary=True` is set.
3. **Given** an `admesh_domains.Mesh` with SMS 2DM format, **When** called, **Then** the correct reader (2DM, not FORT14) is invoked and mesh is loaded correctly.

---

### User Story 2 - Load CHILmesh from plain dict (Priority: P1)

A developer may have mesh metadata in a plain Python dict (e.g., from JSON config or API response) instead of an admesh_domains.Mesh instance. The method should accept both transparently.

**Why this priority**: Maximizes interoperability and allows use without admesh-domains package installed.

**Independent Test**: Can load a mesh from a `{"filename": "...", "type": "ADCIRC"}` dict and get identical result to instance variant.

**Acceptance Scenarios**:

1. **Given** a plain dict `{"filename": "path/to/mesh.fort.14", "type": "ADCIRC", "kind": "mesh"}`, **When** `from_admesh_domain(mesh_dict)` is called, **Then** a CHILmesh is returned with same behavior as instance variant.
2. **Given** a minimal dict with only `"filename"` key (no type/kind), **When** called, **Then** defaults are applied (`type="ADCIRC"`, `kind="mesh"`) and mesh loads successfully.

---

### User Story 3 - Skip layer computation when requested (Priority: P2)

A developer may want to load the mesh topology without the expensive skeletonization step (e.g., for preview, inspection, or batch processing).

**Why this priority**: Supports use cases where skeletonization is not needed; reduces latency for boundary-only access.

**Independent Test**: Can load a mesh with `compute_layers=False`, adjacencies exist, but layers dict is empty or undefined.

**Acceptance Scenarios**:

1. **Given** a mesh record and `compute_layers=False`, **When** `from_admesh_domain(mesh_record, compute_layers=False)` is called, **Then** adjacencies are built but skeletonization is skipped.
2. **Given** `kind="boundary"` (which auto-skips layers regardless of flag), **When** called with `compute_layers=True`, **Then** layers are still skipped (boundary logic overrides flag).

---

### Edge Cases

- What happens if the file referenced in `filename` does not exist? → Raise `FileNotFoundError` with clear message ("Mesh file not found: {path}")
- What if `type` field is missing or unrecognized? → Default to `"ADCIRC"` and attempt fort14 reader
- What if `admesh-domains` package is not installed? → Method should work via duck typing; no import of admesh_domains at module level
- What if `kind` field is missing? → Default to `"mesh"` (not boundary)
- What if mesh is already skeletonized (has layers)? → No re-computation; return as-is

## Requirements

### Functional Requirements

- **FR-001**: Method MUST accept an `admesh_domains.Mesh` instance (duck-typed via `getattr`)
- **FR-002**: Method MUST accept a plain `dict` with at minimum a `"filename"` key
- **FR-003**: Method MUST route to `read_from_2dm()` if `type="SMS_2DM"`, else to `read_from_fort14()`
- **FR-004**: Method MUST accept `compute_layers: bool` parameter (default `True`)
- **FR-005**: Method MUST skip layer computation if `compute_layers=False`
- **FR-006**: Method MUST skip layer computation if `kind="boundary"` (boundary logic overrides flag)
- **FR-007**: Method MUST set internal `_is_boundary` flag to `True` if `kind="boundary"`, else `False`
- **FR-008**: Method MUST raise `FileNotFoundError` with descriptive message if file does not exist
- **FR-009**: Method MUST work without importing `admesh-domains` at module level (no hard dependency)
- **FR-010**: Method MUST default `type` to `"ADCIRC"` if not provided
- **FR-011**: Method MUST default `kind` to `"mesh"` if not provided

### Key Entities

- **Mesh Record**: A data object (admesh_domains.Mesh instance or dict) with fields:
  - `filename` (required): Path to mesh file
  - `type` (optional): Format identifier (`"ADCIRC"` or `"SMS_2DM"`), defaults to `"ADCIRC"`
  - `kind` (optional): Mesh class (`"mesh"` or `"boundary"`), defaults to `"mesh"`

## Success Criteria

### Measurable Outcomes

- **SC-001**: Can load ADCIRC-format meshes from both Mesh instances and dicts without code duplication in caller
- **SC-002**: Can load SMS 2DM-format meshes correctly via `from_admesh_domain()`
- **SC-003**: Boundary meshes skip layer computation and complete in <50% of time vs. full mesh (verified via performance test)
- **SC-004**: `FileNotFoundError` is raised with filename in message for missing files (tested in error scenario)
- **SC-005**: Method functions identically whether `admesh-domains` is installed or not (tested both ways)
- **SC-006**: All existing tests continue to pass (zero regression)

## Assumptions

- **API Stability**: ADMESH-Domains `Mesh` record structure with `filename`, `type`, `kind` fields is stable (or will change with notice)
- **File Paths**: Mesh filenames are absolute or relative to current working directory; no special path resolution needed
- **Reader Availability**: `read_from_2dm()` and `read_from_fort14()` classmethods already exist and work correctly
- **Boundary Flag**: The `_is_boundary` internal flag is the correct mechanism for marking boundary meshes (no validation of boundary-specific logic needed in this method)
- **Lazy vs. Eager**: Skeletonization is computed eagerly (during `from_admesh_domain()` call), not lazily deferred
- **No Material IDs**: Material boundary condition metadata is ignored (as per Issue #44 notes)
