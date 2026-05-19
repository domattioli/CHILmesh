# Testing CHILmesh

Guide for running, writing, and maintaining tests.

## Quick Start

```bash
# Install development dependencies
pip install -e ".[dev]"

# Run full test suite
pytest -v

# Run fast tests only (exclude slow fixtures)
pytest -m "not slow" -v

# Run specific test
pytest tests/test_smoothing.py::TestTriangleSmoother::test_fem_smoother_triangle_preserves_boundary -v
```

## Test Suite Overview

**26 test files | 439 passing | 9 skipped | Runtime: 17.2s**

### Fixtures (Parametrized)

All tests run against 4 built-in meshes via `@pytest.mark.parametrize('mesh', [...])`:

| Fixture | Type | Size | Load Time | Use Case |
|---------|------|------|-----------|----------|
| **annulus** | Triangle | Small | ~0.1s | Quick smoke tests |
| **donut** | Triangle | Medium | ~0.5s | Geometry validation |
| **block_o** | Triangle | Large (100k elems) | ~30s | Performance baseline |
| **structured** | Quad | Medium | ~1s | Mixed-element tests |

### Test Categories

- **Unit tests** (~240): Core geometry, adjacency, skeletonization
- **Integration tests** (~150): Mesh I/O, ADMESH warmstart, layer validation
- **Regression tests** (~30): MATLAB parity, bug fixes
- **Performance tests** (~10): Determinism, edge-building speed

## Running Tests

### All tests
```bash
pytest -v
```

### Exclude slow tests
```bash
pytest -m "not slow" -v        # Skip block_o fixtures
pytest -k "not block_o" -v     # Alternative syntax
```

### Specific test class/function
```bash
pytest tests/test_smoothing.py::TestTriangleSmoother -v
pytest tests/test_smoothing.py::TestTriangleSmoother::test_fem_smoother_triangle_preserves_boundary -v
```

### With coverage
```bash
pip install pytest-cov
pytest --cov=src/chilmesh --cov-report=html tests/
# Open htmlcov/index.html
```

## Writing Tests

### TDD Workflow

1. **Write test first** — describe desired behavior
2. **Run test** — verify it fails (Red)
3. **Implement code** — make test pass (Green)
4. **Refactor** — improve clarity while tests pass

### Test Template

```python
import pytest
from chilmesh import CHILmesh

@pytest.mark.parametrize("mesh", ["annulus", "donut"])
def test_my_feature(mesh):
    """
    Describe what is being tested.
    
    Args:
        mesh: Built-in fixture name (annulus, donut, block_o, structured)
    """
    m = CHILmesh.from_fixture(mesh)
    
    # Arrange: setup
    result = m.some_operation()
    
    # Act: execute
    assert result is not None
    assert len(result) == expected_length, "descriptive assertion message"

@pytest.mark.slow  # Mark slow tests
def test_expensive_operation():
    m = CHILmesh.from_fixture("block_o")
    # expensive test...
```

### Assertions

- **Always include failure messages:**
  ```python
  assert quality > 0.3, f"quality {quality} below threshold"
  ```

- **Use parametrize for fixture variants:**
  ```python
  @pytest.mark.parametrize("mesh", ["annulus", "donut", "block_o"])
  def test_layers_consistent(mesh):
      m = CHILmesh.from_fixture(mesh)
      assert m.n_layers > 0
  ```

- **Test MATLAB parity explicitly:**
  ```python
  def test_layer_count_matches_reference():
      m = CHILmesh.from_fixture("block_o")
      assert m.n_layers == BLOCK_O_LAYER_COUNT_FROM_MATLAB
  ```

## Debugging Tests

### Print debug info
```bash
pytest -v -s tests/test_smoothing.py  # -s = capture stdout
```

### Drop into debugger on failure
```bash
pip install pytest-pdb
pytest --pdb tests/test_smoothing.py  # Drops to (Pdb) on failure
```

### Run single test with verbose output
```bash
pytest -vv -s tests/test_smoothing.py::TestTriangleSmoother::test_fem_smoother_triangle_preserves_boundary
```

## Known Issues & Skipped Tests

**9 tests skipped:** External mesh MATLAB parity tests. Skipped because:
- Large mesh files not bundled (stored separately)
- Used for regression validation, not CI gating
- Run manually for verification: `pytest tests/test_skeletonization_matlab_parity_external.py -v`

## CI & Release Gates

- **Test gate:** All tests must pass before PR merge
- **Coverage gate:** Public API ≥90% coverage (smoke test on matplotlib)
- **Performance gate:** No regression in block_o initialization time

### CI Matrix Policy

| Event | Matrix | Rationale |
|---|---|---|
| **Pull request** | `ubuntu-latest × py3.11`, `pytest -n auto -m "not slow"` | Fast feedback (~2 min). Skips `block_o` fixture. |
| **Push to `main` / `release/**`** | `{ubuntu, macos, windows}-latest × py{3.10, 3.11, 3.12}`, `pytest -n auto -v` | Full portability gate; runs `block_o`. |

Windows is exercised on every merge to `main` (issue #121). PR cycles stay fast — Linux-only — and Windows-specific regressions are caught at the merge gate, not on every push. The `build-and-smoke` job runs ubuntu-only and uses a `$RUNNER_TEMP`-based venv path that is OS-portable (works under any shell).

Platform support claim (constitution §Release Checklist): "All tests pass on Python 3.10 / 3.11 / 3.12 × Ubuntu / macOS / Windows."

## Troubleshooting

### "ModuleNotFoundError: No module named 'chilmesh'"
```bash
pip install -e .
```

### "RuntimeWarning: warm-start truss did not improve quality"
Expected behavior. ADMESH warmstart may return input unchanged if optimization fails. Not a test failure.

### Slow test execution
- Skip block_o: `pytest -m "not slow" -v`
- On first run, pytest caches fixtures (~30s for block_o)
- Subsequent runs use cache

### Flaky tests
Report with: **mesh name**, **fixture type**, **error message**, **reproduction steps**

Example:
```
Flaky: test_determinism on block_o
Error: median quality difference > tolerance
Reproducible: Run 10 times with seed=42
```

## Related

- Issue #110: Test suite holistic audit (.planning/TEST-AUDIT.md)
- Issue #111: Test surface audit + upstream report
- `.specify/memory/constitution.md`: Test-first principle (Principle III)
