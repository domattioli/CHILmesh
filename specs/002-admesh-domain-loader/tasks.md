# Task Checklist: ADMESH-Domains Loader

**Status:** Ready to Start  
**Last Updated:** 2026-05-14  

---

## Phase 1: Design & Contracts (Days 1-2)

### T1.1 Review ADMESH-Domains Structure (2h)
- **Description:** Examine ADMESH-Domains `Mesh` record schema (fields, types, optional/required)
- **Acceptance Criteria:** 
  - [ ] Identified all fields used by CHILmesh integration
  - [ ] Confirmed `filename`, `type`, `kind` are the key attributes
  - [ ] Verified no breaking schema changes in current version
- **Test Scenario:** (N/A — design task)
- **Owner:** TBD
- **Status:** Not started

### T1.2 Finalize Duck-Typing Approach (1h)
- **Description:** Decide which attributes to check for Mesh instance vs. dict duck typing
- **Acceptance Criteria:**
  - [ ] Method uses `getattr(record, 'filename', ...)` for duck typing
  - [ ] Fallback to dict-style access `record.get('filename', ...)`
  - [ ] No type checking with `isinstance(record, admesh_domains.Mesh)`
- **Test Scenario:** (N/A — design task)
- **Owner:** TBD
- **Status:** Not started

### T1.3 Document Type Hints & Signatures (1h)
- **Description:** Write method signature with full type hints
- **Acceptance Criteria:**
  - [ ] Signature: `from_admesh_domain(mesh_record: dict | object, compute_layers: bool = True) -> CHILmesh`
  - [ ] Parameter docstrings explain dict vs. Mesh instance usage
  - [ ] Return type documented
- **Test Scenario:** (N/A — design task)
- **Owner:** TBD
- **Status:** Not started

---

## Phase 2: Implementation (Days 3-5)

### T2.1 Implement `from_admesh_domain()` Classmethod (4h)
- **Description:** Add classmethod to CHILmesh class with full functionality
- **Acceptance Criteria:**
  - [ ] Method accepts both Mesh instance and dict (duck typed)
  - [ ] Parameter validation: `filename` present and valid
  - [ ] Default `type` to `"ADCIRC"` if missing (FR-010)
  - [ ] Default `kind` to `"mesh"` if missing (FR-011)
  - [ ] Raise `FileNotFoundError` with filename in message if file missing (FR-008)
  - [ ] Route to `read_from_2dm()` if `type="SMS_2DM"` (FR-003)
  - [ ] Route to `read_from_fort14()` otherwise (FR-003)
  - [ ] Load mesh successfully (FR-001, FR-002)
- **Test Scenario:** 
  - `from_admesh_domain({"filename": "test_annulus.fort.14"})` returns CHILmesh with correct elements/vertices
  - `from_admesh_domain(mock_mesh_instance)` returns identical result
- **Owner:** TBD
- **Status:** Not started

### T2.2 Implement Boundary & Layer Control Logic (2h)
- **Description:** Handle `_is_boundary` flag and conditional skeletonization
- **Acceptance Criteria:**
  - [ ] Set `_is_boundary=True` if `kind="boundary"` (FR-007)
  - [ ] Set `_is_boundary=False` otherwise (FR-007)
  - [ ] Compute layers if `compute_layers=True` AND `kind != "boundary"` (FR-005, FR-006)
  - [ ] Skip layers if `compute_layers=False` (FR-005)
  - [ ] Skip layers if `kind="boundary"` (FR-006) — boundary logic overrides flag
  - [ ] Docstring explains priority: boundary always skips, else check flag
- **Test Scenario:**
  - `from_admesh_domain({"filename": "...", "kind": "boundary"})` → `mesh._is_boundary=True`, no layers computed
  - `from_admesh_domain({"filename": "..."}, compute_layers=False)` → layers empty/undefined
  - `from_admesh_domain({"filename": "...", "kind": "boundary"}, compute_layers=True)` → layers still skipped
- **Owner:** TBD
- **Status:** Not started

### T2.3 Add Duck-Typing Error Handling (1h)
- **Description:** Handle missing/invalid attributes gracefully
- **Acceptance Criteria:**
  - [ ] If `filename` missing from dict, raise clear error
  - [ ] If `filename` missing from object, raise clear error
  - [ ] If `filename` is empty string, raise error
  - [ ] If `type` unrecognized (not ADCIRC or SMS_2DM), log warning and default to ADCIRC
- **Test Scenario:**
  - `from_admesh_domain({})` → raises ValueError with "filename required"
  - `from_admesh_domain({"filename": "", "type": "ADCIRC"})` → raises error
  - `from_admesh_domain({"filename": "...", "type": "UNKNOWN"})` → warns and defaults
- **Owner:** TBD
- **Status:** Not started

### T2.4 Add Comprehensive Docstring (1h)
- **Description:** Document method with examples and edge cases
- **Acceptance Criteria:**
  - [ ] Docstring includes: purpose, parameters, return type
  - [ ] Example 1: Load from dict
  - [ ] Example 2: Load from Mesh instance
  - [ ] Example 3: Load boundary mesh with layer skipping
  - [ ] Example 4: Load without layer computation
  - [ ] Notes section: explains duck typing, file path resolution, admesh-domains optional
- **Test Scenario:** (N/A — documentation task)
- **Owner:** TBD
- **Status:** Not started

---

## Phase 3: Testing & Validation (Days 6-7)

### T3.1 Write Unit Tests - ADCIRC Format (2h)
- **Description:** Test loading ADCIRC-format meshes (fort.14)
- **Acceptance Criteria:**
  - [ ] Test case: Load from dict with ADCIRC type → mesh loads correctly (SC-001)
  - [ ] Test case: Load from Mesh instance → same result as dict (SC-001)
  - [ ] Test case: Default type to ADCIRC if missing → works (FR-010)
  - [ ] Test case: Adjacencies and elements match original mesh (FR-001, FR-002)
- **Test Scenario:** Use annulus fixture
  - `m1 = from_admesh_domain({"filename": "annulus.fort.14"})`
  - `m2 = CHILmesh.read_from_fort14("annulus.fort.14")`
  - Assert `m1.elements == m2.elements`, `m1.points == m2.points`
- **Owner:** TBD
- **Status:** Not started
- **File:** `tests/test_admesh_domains_loader.py`

### T3.2 Write Unit Tests - SMS 2DM Format (1h) [If 2DM support exists]
- **Description:** Test loading SMS 2DM-format meshes
- **Acceptance Criteria:**
  - [ ] Test case: Load from dict with SMS_2DM type → mesh loads correctly (SC-002)
  - [ ] Test case: Correct reader invoked (2DM, not fort14) (FR-003)
- **Test Scenario:** If structured fixture is SMS 2DM format
  - `m = from_admesh_domain({"filename": "structured.2dm", "type": "SMS_2DM"})`
  - Assert valid mesh returned
- **Owner:** TBD
- **Status:** Not started
- **File:** `tests/test_admesh_domains_loader.py`

### T3.3 Write Unit Tests - Boundary Mesh Handling (1.5h)
- **Description:** Test boundary mesh logic and layer skipping
- **Acceptance Criteria:**
  - [ ] Test case: Load boundary mesh → `_is_boundary=True` (FR-007)
  - [ ] Test case: Boundary mesh skips layers (SC-003) even with `compute_layers=True` (FR-006)
  - [ ] Test case: Regular mesh with `compute_layers=False` → no layers (FR-005)
  - [ ] Performance: Boundary mesh load < 50% time of full mesh (SC-003)
- **Test Scenario:**
  - `m = from_admesh_domain({"filename": "annulus.fort.14", "kind": "boundary"})`
  - Assert `m._is_boundary == True`, `len(m.layers) == 0`
  - Measure load time, verify <50% vs full mesh
- **Owner:** TBD
- **Status:** Not started
- **File:** `tests/test_admesh_domains_loader.py`

### T3.4 Write Unit Tests - Error Cases (1h)
- **Description:** Test error handling for invalid inputs
- **Acceptance Criteria:**
  - [ ] Test case: Missing file → FileNotFoundError with filename (FR-008, SC-004)
  - [ ] Test case: Missing filename field → ValueError with clear message
  - [ ] Test case: Empty filename → ValueError
  - [ ] Test case: Invalid type field → logs warning, defaults to ADCIRC
- **Test Scenario:**
  - `from_admesh_domain({"filename": "/nonexistent/mesh.fort.14"})` → raises FileNotFoundError("Mesh file not found: /nonexistent/mesh.fort.14")
  - `from_admesh_domain({})` → raises ValueError with "filename"
- **Owner:** TBD
- **Status:** Not started
- **File:** `tests/test_admesh_domains_loader.py`

### T3.5 Test Duck Typing - Mesh Instance Compatibility (1h)
- **Description:** Test with real/mock admesh_domains.Mesh objects
- **Acceptance Criteria:**
  - [ ] Method works with mock Mesh object (attributes: filename, type, kind)
  - [ ] Method works without admesh_domains installed (SC-005)
  - [ ] Duck typing requires no isinstance() checks
- **Test Scenario:**
  - Create mock object with `filename`, `type`, `kind` attributes
  - Pass to method, verify loads correctly
  - Verify no import of admesh_domains in method body
- **Owner:** TBD
- **Status:** Not started
- **File:** `tests/test_admesh_domains_loader.py`

### T3.6 Run Full Regression Test Suite (1h)
- **Description:** Ensure no regressions in existing tests
- **Acceptance Criteria:**
  - [ ] `pytest tests/ -v` passes with no failures (SC-006)
  - [ ] No new test warnings
  - [ ] Coverage for new code ≥85%
- **Test Scenario:**
  - `pytest tests/ -v` from repo root
  - `pytest tests/test_admesh_domains_loader.py --cov=src/chilmesh`
- **Owner:** TBD
- **Status:** Not started
- **File:** All tests

---

## Phase 4: Documentation & Review (Day 8)

### T4.1 Update README or Examples (0.5h)
- **Description:** Add example of loading from ADMESH-Domains
- **Acceptance Criteria:**
  - [ ] README or examples/ file shows `from_admesh_domain()` usage
  - [ ] Example shows both dict and Mesh instance variants
  - [ ] Example shows boundary mesh handling
- **Owner:** TBD
- **Status:** Not started

### T4.2 Code Review & Feedback Integration (0.5h)
- **Description:** Address peer review comments
- **Acceptance Criteria:**
  - [ ] All review comments resolved or explicitly noted as deferred
  - [ ] Code follows CHILmesh style (per CLAUDE.md, constitution.md)
  - [ ] Type hints correct and verified with mypy (if available)
- **Owner:** TBD
- **Status:** Not started

### T4.3 Update CHANGELOG (0.5h)
- **Description:** Document feature in CHANGELOG.md
- **Acceptance Criteria:**
  - [ ] Entry under 0.3.0 (or next version)
  - [ ] Describes `from_admesh_domain()` new feature
  - [ ] Links to spec or issue
- **Owner:** TBD
- **Status:** Not started

---

## Summary

| Phase | Task Count | Est. Hours | Owner | Status |
|-------|-----------|-----------|-------|--------|
| 1: Design | 3 | 4h | TBD | Not started |
| 2: Implementation | 4 | 8h | TBD | Not started |
| 3: Testing | 6 | 7h | TBD | Not started |
| 4: Docs | 3 | 1.5h | TBD | Not started |
| **TOTAL** | **16** | **20.5h** | — | — |

*Estimate updated from plan.md (23h total); testing effort refined based on fixture availability.*

---

## Progress Tracking

- [ ] All Phase 1 tasks complete
- [ ] All Phase 2 tasks complete
- [ ] All Phase 3 tasks complete
- [ ] All Phase 4 tasks complete
- [ ] Merged to main; released in 0.3.0
