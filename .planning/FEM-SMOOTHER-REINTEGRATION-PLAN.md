# FEM Smoother Re-Integration Plan

**Issue:** #105 (Re-integrate FEM smoother work with proper documentation)  
**Context:** PR #104 merged, then reverted during documentation compression (PR #97)  
**Status:** Blocked on understanding prior state + finding preserved code

---

## Current Understanding

### What Happened
1. **PR #104:** FEM smoother quad-stiffness symmetry fix + mixed-element support + regression tests
   - Fixed: asymmetric diagonal decomposition causing element collapse
   - Added: `_quad_stiffness_assembly()` symmetric algorithm
   - Tests: 7 new regression tests in `TestMixedElementFEMSmoother`
   - Results: min quality 0.187 → 0.274 (+46.5%)

2. **PR #97:** Documentation compression (caveman-mode)
   - Reverted: FEM smoother changes (to match compressed docs baseline)
   - Effect: Lost quad/mixed-element FEM fixes

3. **Now:** Code exists on branch, docs compressed; need clean re-merge

### Key Changes to Restore
```
src/chilmesh/CHILmesh.py:
  - _quad_stiffness_assembly() → symmetric dual-diagonal averaging
  - _mixed_stiffness_assembly() → handles tri + quad cohesively
  - boundary_edges() → lazy-compute for compute_layers=False
  - redistribute_outer_perim() → h_min 0.05→0.18 (eliminate degenerates)

tests/test_smoothing.py:
  - TestMixedElementFEMSmoother (7 new tests)
  - test_quad_stiffness_is_symmetric
  - test_detect_element_types_mixed
  - test_fem_smoother_mixed_no_element_collapse
  - test_fem_smoother_mixed_interior_moves_toward_equilibrium
  (+ 3 more)

scripts/generate_mixed_truss_demo.py:
  - Updated visualization (wireframe, skeletonization, quality colormaps)
  - Correct plot_layer() invocation with autoscale_view()

docs/README:
  - FEM smoother section (restored + compressed)
  - Mixed-element capability documentation
  - Performance/quality results
```

---

## Integration Blockers

### 1. Code Source Uncertainty
**Blocker:** Where is the preserved FEM smoother code?
- PR #104 branch deleted?
- Code on `daily-maintenance` branch?
- Need to locate + verify correctness

**Resolution:** Requires:
- Git search: `git log --all --grep="FEM\|quad.*stiffness"` → find commit SHAs
- Branch inspection: `git show <branch>:src/chilmesh/CHILmesh.py` → verify code present

### 2. Merge Conflict Risk
**Risk:** Documentation was compressed post-PR #104; merging old code into compressed docs = conflict
- Issue #97 removed/simplified docstrings
- FEM smoother tests added multiline docstrings
- Mixed-element code is verbose

**Mitigation:**
- Rebase PR #104 onto current main (not just cherry-pick)
- Resolve conflicts: keep current docstring style (compressed/caveman-mode)
- Update test docstrings to match current convention

### 3. Test Regression Check
**Blocker:** Current test suite (439 passing); adding 7 FEM tests = 446
- Need to verify: tests still pass + new tests pass
- Parametrization: all 4 fixtures (annulus, donut, block_o, structured)

### 4. Documentation Sync
**Blocker:** README, CHANGELOG, docstrings must align
- Current: Compressed (caveman-mode)
- PR #104 added: Verbose multiline docs

**Resolution:**
- Extract documentation from PR #104
- Re-write in caveman style (fragments, no fluff)
- Cross-ref with Constitution Principle VIII (docs = contract)

---

## Integration Steps (3-4 hours)

### Phase 1: Locate Code (30 min)
```bash
# Find FEM smoother commits
git log --all --oneline --grep="FEM\|quad.*stiffness\|PR #104" | head -20

# Inspect branch
git show origin/daily-maintenance:src/chilmesh/CHILmesh.py | grep "_quad_stiffness_assembly" | head -5

# Or inspect PR:
git log --all --grep="104" --format="%h %s" | head -5
git show <SHA>:src/chilmesh/CHILmesh.py
```

### Phase 2: Extract Code (1 hour)
- Copy `_quad_stiffness_assembly()`, `_mixed_stiffness_assembly()`, related methods
- Identify test cases from PR #104
- List documentation changes needed

### Phase 3: Merge onto Current Main (1 hour)
```bash
# Option A: Cherry-pick commits
git cherry-pick <FEM-SHA>..<END-SHA>

# Option B: Rebase preserved branch
git rebase main origin/daily-maintenance -- tests/ src/

# Resolve conflicts (docstrings, test structure)
```

### Phase 4: Verify + Document (1-1.5 hours)
- Run full test suite: `pytest -v`
- Verify: all 439 original + 7 new FEM tests pass
- Update CHANGELOG: summarize FEM smoother changes
- Update README: restore FEM section (caveman-style)
- Ensure Constitution Principle VIII (docs = contract) satisfied

---

## Success Criteria

- [ ] Code restored: `_quad_stiffness_assembly()`, `_mixed_stiffness_assembly()` in CHILmesh.py
- [ ] Tests pass: all 446 (439 original + 7 FEM new)
- [ ] No test regressions on block_o (critical large mesh)
- [ ] Documentation updated: README + docstrings in caveman-mode
- [ ] Changelog entry: documents FEM smoother quad support, mixed-element improvements, quality gains
- [ ] Constitution Principle VIII: docs match implementation

---

## Known Issues from PR #104

**Visualization bugs fixed by PR #104:**
- Left panel: wireframe only (white elements, black edges)
- Middle panel: skeletonization layers (requires `autoscale_view()` after `add_collection()`)
- Right panel: element quality (cool colormap)

**These fixes must be restored in `scripts/generate_mixed_truss_demo.py`**

---

## Recommendation

**Priority:** High (unblocks MADMESHR adaptive workflows)  
**Effort:** 3-4 hours (code extraction, merge, testing, docs)  
**Risk:** Medium (merge conflicts likely; needs careful conflict resolution)  
**Timeline:** Next session after reviewing PR #104 history

**Precondition:** Locate preserved code (git search in Phase 1)

---

## Related

- Issue #105: FEM smoother re-integration (this plan)
- Issue #95, #100: Original bugs (quad collapse, mixed-element quality)
- PR #104: Original fix (now reverted)
- PR #97: Documentation compression (reason for revert)
- Issue #103: Code structure audit (relates to monolithic CHILmesh.py)
- MADMESHR: Primary downstream consumer (adaptive refinement)

---

**Next Step:** Execute Phase 1 (locate code) to unblock phases 2-4.
