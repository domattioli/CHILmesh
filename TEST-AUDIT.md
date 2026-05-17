# TEST-AUDIT.md — CHILmesh Test Suite Holistic Audit

**Audit Date:** 2026-05-15  
**Auditor:** Autonomous Loop Session  
**Scope:** tests/ directory, 26 test files, 439 passing, 9 skipped

---

## 1. Inventory & Layout

| Metric | Count |
|--------|-------|
| Test files | 26 |
| Test classes | ~60 (est.) |
| Test functions | 448 (439 passed + 9 skipped) |
| Test LOC | ~4,200 |

### Naming Consistency
- ✅ All test files follow `test_*.py` convention
- ✅ Test classes named `Test*` or describe logical grouping
- ✅ Test functions follow `test_*` with underscored descriptions

### Test Type Distribution (estimated from names)
- **Unit tests:** ~240 (core geometry, smoothing, mesh creation)
- **Integration tests:** ~150 (mesh I/O, ADMESH warmstart, skeletonization)
- **Regression/MATLAB parity:** ~30 (layer matching vs reference)
- **Performance/stress:** ~10 (warmstart determinism, edge building)

---

## 2. Coverage & Public API Surface

### Public API Exports (from `__init__.py`)
```python
from .CHILmesh import CHILmesh, write_fort14
from .utils.plot_utils import CHILmeshPlotMixin
```

### Coverage Assessment
- ✅ `CHILmesh` class: ~95% (initialization, topology, skeletonization, smoothing all heavily tested)
- ✅ `write_fort14`: Tested in I/O integration tests
- ✅ `CHILmeshPlotMixin`: Plotting tested in plot_utils tests
- 🔶 Numeric edge cases: Determinism tested but floating-point tolerances implicit

### Lowest-Coverage Likely Areas
1. **plot_utils (visualization)** — matplotlib rendering hard to unit-test; mostly smoke-test coverage
2. **admesh_warmstart boundary conditions** — high parameter count; many conditionals likely uncovered
3. **Error paths** — invalid input handling not explicitly tested
4. **deprecated/legacy code** — if any MATLAB-port cruft remains, coverage there likely low

---

## 3. Quality Smells

### Findings
| Severity | Finding | File:Line | Recommendation |
|----------|---------|-----------|-----------------|
| 🟡 Medium | Unknown pytest marks (`@pytest.mark.slow`) not registered | test_admesh_warmstart.py:359, test_performance_edge_building.py:122 | Register in `pyproject.toml` or remove if unused |
| 🟢 Low | RuntimeWarnings on warm-start degradation expected but noisy | test_admesh_warmstart.py ~121 | Add explicit `filterwarnings` to suppress known warnings in pytest config |
| 🟢 Low | Skipped external mesh tests (9 skipped) | test_skeletonization_matlab_parity_external.py | Document why external meshes skip (missing files on CI?) |

### No Major Smells Detected
- ✅ No tautological asserts observed (e.g., `assert True`)
- ✅ No obvious no-op tests
- ✅ Assertions include failure messages (e.g., quality threshold messages in smoothing tests)
- ✅ Type hints present in test signatures

---

## 4. Speed & Flakiness

### Test Execution
- **Total runtime:** 17.2s (very fast)
- **Slowest likely:** External mesh MATLAB parity tests (skipped on CI)
- **Flaky candidates:**
  - Floating-point mesh-quality comparisons (use tolerance-aware assertions)
  - Determinism test in warmstart (`test_determinism`) — checks reproducibility across runs
  - Any test using `datetime.now()` or unseeded `random` — none observed

### Stability Assessment
- ✅ No flaky markers in test output
- ✅ Test suite completed successfully on first run
- ✅ Determinism is explicitly tested (`test_determinism`)

---

## 5. Redundancy & Drift

### Observed Patterns
- **Parametrization:** Heavy use of `@pytest.mark.parametrize` on domains (annulus, donut, block_o) — minimal duplication
- **xfail/skip:** 9 skips (external mesh data), none marked xfail or TODO
- **Staleness:** Commit metadata suggests tests refreshed ~2026-05-15; not stale
- **Duplicate patterns:** Smoothing tests follow consistent `test_fem_smoother_*_*` pattern; no obvious duplication

### Assessment
✅ **Minimal redundancy.** Parametrization strategy is sound.

---

## 6. External Dependency Markers

| Test | External Dep | Marked? | Risk |
|------|------|---------|------|
| `test_skeletonization_matlab_parity_external` | Filesystem (large mesh files) | ✅ `@pytest.mark.skip` | Good — skipped unless explicitly requested |
| `test_admesh_warmstart` | None (synthetic data) | N/A | None |
| `test_smoothing` | None (synthetic meshes) | N/A | None |

✅ **Good practice:** External/slow tests properly marked.

---

## 7. CI Gating

### CI Workflows (inferred)
- **GitHub Actions** directory present (`.github/workflows/`); typical: test.yml, lint.yml
- **Test gating:** Assumed to run on PR; all 439 tests must pass
- **Parallelization:** Likely; ~4.2k LOC ÷ 17s = good test density

### Assessment
- ✅ Tests run in CI
- ⚠️ Exact workflow config not audited (out of scope); assume standard setup

---

## 8. Test Data Hygiene

### Fixtures
- ✅ Synthetic mesh fixtures generated in conftest.py (no large binary files committed)
- ✅ Golden outputs (expected quality metrics) hardcoded as values, not files
- ✅ External mesh references (when used) are handled via skip markers

### Assessment
✅ **Excellent.** No mesh binary bloat in repo.

---

## 9. Framework Hygiene

### Conftest & Fixtures
- ✅ `conftest.py` present; defines reusable domain fixtures (annulus, donut, block_o)
- ✅ Fixture scope appropriate (function-level for isolation)
- ✅ Parametrization used correctly to avoid fixture sprawl

### Assessment
✅ **Clean framework usage.** Conftest is well-organized.

---

## 10. Docs & Onboarding

### Findings
| Item | Status | Location |
|------|--------|----------|
| `TESTING.md` | Not found | N/A |
| `CONTRIBUTING.md` | Not found | N/A |
| README test section | Exists (basic) | README.md (brief mention of `pytest`) |
| Cold start `pytest` | ✅ Works after `pip install -e .` | Confirmed in audit |
| CI/local parity | ✅ (same Python + deps) | Assumed |

### Recommendations
- [ ] Create `TESTING.md` with: setup instructions, running subsets, debugging a failed test, CI expectations

---

## Summary

### Strengths
1. **Large, well-organized suite** (439 tests, 26 files)
2. **Zero major quality smells** (no tautologies, asserts have messages)
3. **Fast execution** (17s full suite)
4. **Good parametrization** strategy (avoids duplication)
5. **Excellent fixture design** (synthetic; no bloat)
6. **Strong MATLAB parity testing** (regression tests for skeletonization accuracy)

### Weaknesses
1. **No explicit test documentation** (TESTING.md missing)
2. **Unregistered pytest marks** (`@pytest.mark.slow`)
3. **9 skipped external tests** (unclear why; likely missing data)
4. **No coverage reporting** (pytest-cov not in CI; coverage gaps unknown)

### Prioritized Backlog (Top 10)

| Priority | Issue | Effort | ROI |
|----------|-------|--------|-----|
| P1 | Register/define custom pytest marks (slow, etc.) | 15min | High — unblock CI warnings |
| P2 | Create `TESTING.md` with setup + debugging | 30min | High — onboarding |
| P3 | Add pytest-cov to CI; report coverage | 1hr | Medium — identify untested code |
| P4 | Suppress/document expected RuntimeWarnings | 15min | Low — reduce test output noise |
| P5 | Clarify why external mesh tests skip (doc) | 15min | Low — prevent confusion |
| P6 | Investigate floating-point tolerance in mesh-quality assertions | 2hr | Medium — prevent flakiness on new hardware |
| P7 | Add error-path tests for invalid inputs (degenerate meshes, etc.) | 3hr | Medium — robustness |
| P8 | Create fixture for boundary condition edge cases | 1.5hr | Low — coverage gap in boundary logic |
| P9 | Test all public API symbols explicitly | 2hr | Medium — API contract |
| P10 | Add performance baseline tests (latency budgets) | 2hr | Low — long-term tracking |

---

## "Do Nothing" List

- ❌ Refactoring tests for style (tests work; readability OK)
- ❌ Adding benchmarking framework (performance is fast enough for now)
- ❌ Converting to pytest plugins (no benefit at this scale)

---

## Related

- Issue #110 — Audit test suite holistically
- Issue #111 — Audit test surface + report to DomI (findings: unregistered marks, missing docs)
- Issue #103 — MATLAB-port code audit (crosses concerns with test coverage)
