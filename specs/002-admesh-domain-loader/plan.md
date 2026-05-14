# Implementation Plan: ADMESH-Domains Loader

**Feature:** Create `from_admesh_domain()` classmethod  
**Target Version:** 0.3.0  
**Owner:** CHILmesh Team  
**Created:** 2026-04-26  

---

## High-Level Approach

Implement a classmethod `CHILmesh.from_admesh_domain(mesh_record, compute_layers=True)` that:
1. Accepts both `admesh_domains.Mesh` instances and plain dicts (duck typing)
2. Routes to appropriate file reader (`read_from_fort14()` for ADCIRC, `read_from_2dm()` for SMS)
3. Sets internal state (`_is_boundary` flag) based on mesh `kind`
4. Conditionally computes skeletonization layers (skips for boundaries)
5. Works without hard dependency on `admesh-domains` package

**Key Principle:** Bridge layer between ADMESH-Domains and CHILmesh, enabling seamless integration for downstream projects.

---

## Phase Breakdown

### Phase 1: Design & Contracts (Days 1-2)

**Goal:** Define method signature, API contracts, and integration points.

**Tasks:**
- Review ADMESH-Domains `Mesh` record structure (fields, types, constraints)
- Finalize duck-typing approach (which attributes to check)
- Document type hints and signatures
- Update type stubs if needed

**Deliverables:**
- Finalized method signature
- Integration contract documentation
- Type hints ready for implementation

**Success:** Design review passed; no ambiguity in implementation requirements

---

### Phase 2: Implementation (Days 3-5)

**Goal:** Implement classmethod with full functionality for both Mesh instances and dicts.

**Tasks:**
- Implement `from_admesh_domain()` classmethod in `CHILmesh.py`
- Add parameter validation and default handling
- Route file reader based on `type` field
- Set `_is_boundary` flag correctly
- Handle `compute_layers` parameter and boundary override logic
- Add error handling for missing files and invalid records

**Deliverables:**
- Working classmethod with all FR features
- Docstring with examples
- Type hints for all parameters/returns

**Success:** All unit tests pass (see tasks.md)

---

### Phase 3: Testing & Validation (Days 6-7)

**Goal:** Comprehensive testing across fixtures and edge cases.

**Tasks:**
- Write parametrized tests for both ADCIRC and SMS 2DM formats (if available)
- Test Mesh instance variant vs. dict variant for equivalence
- Test boundary mesh handling (layer skipping)
- Test `compute_layers=False` parameter
- Test error cases (missing file, bad format, missing fields)
- Test that method works without `admesh-domains` installed
- Run full test suite for regression

**Deliverables:**
- ≥10 test cases covering all user stories and edge cases
- Performance baseline: boundary mesh load time <50% of full mesh (SC-003)

**Success:** All tests pass; SC-001 through SC-006 verified

---

### Phase 4: Documentation & Review (Days 8)

**Goal:** Final documentation and peer review.

**Tasks:**
- Add method to README examples or docstring
- Write integration guide for ADMESH-Domains users
- Code review and feedback incorporation
- Update CHANGELOG

**Deliverables:**
- Updated docstrings
- Integration example in tests or examples/
- CHANGELOG entry

**Success:** Peer review approved; documentation complete

---

## Milestones

| Milestone | Target | Gate |
|-----------|--------|------|
| **Design Complete** | Day 2 EOD | Integration contract approved |
| **Implementation Complete** | Day 5 EOD | All unit tests passing |
| **Testing Complete** | Day 7 EOD | All functional tests passing, no regressions |
| **Shipped** | Day 8 EOD | Merged to main; released in 0.3.0 |

---

## Known Blockers & Dependencies

### No Hard Blockers

- ✅ `read_from_fort14()` and `read_from_2dm()` already exist
- ✅ `_is_boundary` flag mechanism already established

### Soft Dependencies

- Performance baseline in SC-003 assumes `_skeletonize()` is main cost center (true per v0.2.0 benchmarks)

---

## Risk Assessment

### Low Risk

- ✅ New method, no changes to existing APIs
- ✅ Duck typing (no hard dependency on admesh-domains) limits integration surface
- ✅ Defaults preserve backward compatibility
- ✅ Boundary flag already used internally

### Medium Risk

- 🟡 **ADMESH-Domains schema change:** If Mesh record structure changes, method breaks. Mitigation: document expectations in docstring; add version check if needed.
- 🟡 **File path resolution:** Assumes absolute/relative paths work. Mitigation: test with various path formats.

### Mitigation Summary

- Add integration tests with real ADMESH-Domains objects (if available)
- Document assumptions clearly in docstring
- Add `kind="boundary"` override logic to prevent accidental layer computation

---

## Success Metrics

| Metric | Target | How to Measure |
|--------|--------|-----------------|
| **Feature Complete** | Day 8 | All FR-001 through FR-011 working |
| **Test Coverage** | ≥85% | Run pytest --cov on test_admesh_domains_loader.py |
| **Backward Compat** | 100% | Full test suite passes without modification |
| **User Satisfaction** | TBD | ADMESH-Domains integration users confirm seamless loading |

---

## Effort Estimate

| Phase | Estimate | Confidence |
|-------|----------|------------|
| Design | 4h | High (straightforward contract) |
| Implementation | 8h | High (no algorithm complexity) |
| Testing | 8h | Medium (edge cases may appear) |
| Documentation | 3h | High (straightforward) |
| **Total** | **23h** | **High** |

---

## Next Steps

1. Get design approval (contract review)
2. Assign implementation owner
3. Create tasks.md from this plan
4. Begin Phase 1 tasks
