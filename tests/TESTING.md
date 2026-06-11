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

**42 test files | 982 tests collected**

- Fast PR mode (`pytest -n auto -m "not slow"`): **928 passed, 52 skipped (~27s)**
- Full suite adds 2 `slow`-marked `block_o` tests, run on push-to-`main` / `release/**`

### Fixtures (Parametrized)

Most tests run against the built-in meshes via `@pytest.mark.parametrize('mesh', [...])`.
`conftest.py` exposes 5 fixtures (`FIXTURE_NAMES`); the 4 triangle fixtures form
`TRI_FIXTURE_NAMES` for triangle-only parametrizations:

| Fixture | Type | Size | Load Time | Use Case |
|---------|------|------|-----------|----------|
| **annulus** | Triangle | Small | ~0.1s | Quick smoke tests |
| **donut** | Triangle | Medium | ~0.5s | Geometry validation |
| **block_o** | Triangle | Large (100k elems) | ~30s | Performance baseline (`slow`) |
| **structured** | Quad | Medium | ~1s | Mixed-element tests |
| **quad_2x2** | Quad | Tiny (4 elems) | <0.01s | Quad/mixed unit cases |

### Test Categories

- **Unit tests**: Core geometry, adjacency, skeletonization, smoothing
- **Integration tests**: Mesh I/O, ADMESH warmstart, layer validation, backend equivalence
- **Regression tests**: MATLAB parity, fixed-bug reproducers
- **Performance tests**: Determinism, edge-building speed

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

**52 tests skipped (fast PR mode).** Skips are environment- or geometry-conditional, not failures:
- **External MATLAB parity** (`test_skeletonization_matlab_parity_external.py`): require `admesh-domains` + large mesh files not bundled. Run manually: `pip install admesh-domains && pytest tests/test_skeletonization_matlab_parity_external.py -v`
- **C++ backend equivalence** (`test_backend_equivalence.py`): skipped when `chilmesh_cpp` extension is not built (`pip install -e ".[cpp]"` or build the extension to exercise these).
- **Geometry-conditional** (e.g. `test_spatial_indexing.py`): point-location cases that don't apply to holed fixtures (annulus/donut/block_o).

## CI & Release Gates

- **Test gate:** All tests must pass before PR merge
- **Coverage gate:** `--cov-fail-under=80` on push-to-`main` / `release/**`. PR runs skip coverage to keep feedback fast (#122, TEST-AUDIT F13). Floor was set after F1 raised `CHILmesh.py` line coverage from 73 % → 89 % and overall to 83 %.
- **Performance gate:** No regression in block_o initialization time

### CI Matrix Policy

| Event | Matrix | Rationale |
|---|---|---|
| **Pull request** | `ubuntu-latest × py3.11`, `pytest -n auto -m "not slow"` | Fast feedback (~2 min). Skips `block_o` fixture. |
| **Push to `main` / `release/**`** | `{ubuntu, macos, windows}-latest × py{3.10, 3.11, 3.12}`, `pytest -n auto -v` | Full portability gate; runs `block_o`. |

Windows is exercised on every merge to `main` (issue #121). PR cycles stay fast — Linux-only — and Windows-specific regressions are caught at the merge gate, not on every push. The `build-and-smoke` job runs ubuntu-only and uses a `$RUNNER_TEMP`-based venv path that is OS-portable (works under any shell).

#### xdist worker semantics

`pytest -n auto` spawns one worker process per CPU. Each worker reimports `tests/conftest.py`, so the session-scope mesh cache (`_MESH_CACHE`) is per-worker — no cross-worker race, but each worker pays the `block_o` warmup once on push-to-main runs (PR runs exclude `slow`-marked fixtures). PR-mode local run: `842 passed, 47 skipped, 26.4s` (verified 2026-05-23, issue #122).

#### Runner cache hygiene

- `actions/setup-python` uses `check-latest: true` so a cached interpreter is replaced when GitHub publishes a newer patch release for the requested minor — `pip cache` still hits because `cache-dependency-path: pyproject.toml`.
- `MPLCONFIGDIR` is set to a per-run `$RUNNER_TEMP` subdir on every job. matplotlib's first-run font-cache build dominates macOS wall-clock; routing it to tmp keeps the cost predictable and prevents `$HOME` writes on shared runners.

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
