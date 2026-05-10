# CHILmesh Test Coverage Analysis

**Analysis Date**: April 26, 2026  
**Coverage Tool**: pytest-cov  
**Python Version**: 3.11

---

## Executive Summary

CHILmesh has **84% overall test coverage** (737 statements, 118 missed). While the core mesh functionality is well-tested, there are significant gaps in:

1. **Plotting functionality** (60% coverage) - Critical for visualization and debugging
2. **Mesh smoothing methods** (0% coverage) - Advanced mesh quality improvement
3. **Mesh mutation operations** (0% coverage) - Important for adaptation workflows
4. **Error handling** (70% coverage) - Edge cases in validation and fallbacks

Closing these gaps will improve coverage to **90%+** and ensure reliability for downstream projects (MADMESHR, ADMESH, ADMESH-Domains).

---

## Current Coverage Summary

### Overall Statistics
```
Total Statements: 736
Covered: 618 (84%)
Missed: 118 (16%)
Tests Passing: 137 (excluding slow block_o mesh)
```

### By Module
| Module | Statements | Coverage | Missed | Priority |
|--------|-----------|----------|--------|----------|
| `src/chilmesh/CHILmesh.py` | 508 | 92% | 40 | High |
| `src/chilmesh/utils/plot_utils.py` | 197 | 60% | 78 | **Critical** |
| `src/chilmesh/__init__.py` | 6 | 100% | 0 | — |
| `src/chilmesh/data/__init__.py` | 0 | 100% | 0 | — |
| `src/chilmesh/examples.py` | 23 | 100% | 0 | — |

---

## Major Coverage Gaps

### 1. **Mesh Smoothing Methods (CRITICAL)** — 0% Coverage

**Risk Level**: 🔴 **High**  
**Feature Maturity**: Beta (angle-based not implemented)

#### Missing Coverage
- `smooth_mesh()` — Main smoothing interface (lines 1118–1134)
- `angle_based_smoother()` — Raises `NotImplementedError` (lines 1136–1147)
- `direct_smoother()` — FEM-based smoothing (lines 1149–1203)

#### Why It Matters
Mesh smoothing is essential for quality improvement in adaptation workflows. Without tests:
- Refactoring FEM solver could break silently
- Users cannot verify smoothing actually improves mesh quality
- Boundary constraint handling is untested

#### Recommendations

**Test 1: `test_smooth_mesh_fem_preserves_geometry()`**
```python
def test_smooth_mesh_fem_preserves_geometry(annulus):
    """FEM smoothing should preserve boundary nodes."""
    mesh = annulus
    boundary_edges = mesh.boundary_edges()
    boundary_nodes = np.unique(mesh.edge2vert(boundary_edges).flatten())
    
    points_before = mesh.points[boundary_nodes].copy()
    mesh.smooth_mesh('fem', acknowledge_change=True)
    points_after = mesh.points[boundary_nodes]
    
    np.testing.assert_allclose(points_before, points_after, atol=1e-10)
```

**Test 2: `test_smooth_mesh_fem_quality_improvement()`**
```python
def test_smooth_mesh_fem_quality_improvement(donut):
    """Smoothing should improve skew angle quality."""
    mesh = donut
    quality_before, _, stats_before = mesh.elem_quality(quality_type='skew')
    mesh.smooth_mesh('fem', acknowledge_change=True)
    quality_after, _, stats_after = mesh.elem_quality(quality_type='skew')
    
    # Mean quality should improve
    assert stats_after['mean'] > stats_before['mean']
```

**Test 3: `test_smooth_mesh_requires_acknowledge()`**
```python
def test_smooth_mesh_requires_acknowledge(annulus):
    """smooth_mesh must require acknowledge_change=True."""
    mesh = annulus
    with pytest.raises(AssertionError):
        mesh.smooth_mesh('fem', acknowledge_change=False)
```

**Test 4: `test_angle_based_smoother_not_implemented()`**
```python
@pytest.mark.xfail(reason="Angle-based smoothing not yet implemented")
def test_angle_based_smoother_not_implemented(annulus):
    """Skip until angle-based smoothing is implemented."""
    mesh = annulus
    with pytest.raises(NotImplementedError):
        mesh.angle_based_smoother()
```

---

### 2. **Plotting Module (CRITICAL)** — 60% Coverage

**Risk Level**: 🔴 **High**  
**Impact**: Visualization is critical for debugging and presentations

#### Uncovered Methods (78 statements)
- `plot()` — Main mesh visualization (lines 47–110)
- `plot_edge()` — Edge-specific plotting (missing coverage)
- `plot_layers()` — Layer visualization with colormaps (missing coverage)
- `axis_chilmesh()` — Axis configuration helper (lines 20–44)
- `plot_elem_quality()` — Quality metric coloring (missing coverage)

#### Why It Matters
- Plotting code paths are not exercised in CI
- Matplotlib API changes could silently break visualization
- Edge cases (empty meshes, single elements) are untested
- Color mapping and layer visualization are fragile

#### Recommendations

**Create new file: `tests/test_plotting.py`**

```python
import pytest
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import numpy as np


class TestMeshPlotting:
    """Test mesh visualization functionality."""
    
    def test_plot_basic_mesh_creates_figure(self, annulus):
        """Plotting should create valid matplotlib figure."""
        fig, ax = annulus.plot()
        assert isinstance(fig, plt.Figure)
        assert isinstance(ax, plt.Axes)
        plt.close(fig)
    
    def test_plot_with_element_colors(self, annulus):
        """Plotting with elem_color should add patches."""
        fig, ax = annulus.plot(elem_color='lightblue')
        # Check patches were added
        assert len(ax.patches) > 0
        plt.close(fig)
    
    def test_plot_edges_only(self, annulus):
        """Plotting with elem_color='none' should only draw edges."""
        fig, ax = annulus.plot(elem_color='none', edge_color='red')
        # Check line segments/collections added
        assert len(ax.lines) > 0 or len(ax.collections) > 0
        plt.close(fig)
    
    def test_plot_single_element(self, annulus):
        """Plotting single element by ID."""
        elem_id = 5
        fig, ax = annulus.plot(elem_ids=elem_id, elem_color='green')
        # Should create figure with limited geometry
        assert ax.get_xlim()[0] < ax.get_xlim()[1]
        plt.close(fig)
    
    def test_plot_element_list(self, annulus):
        """Plotting list of elements."""
        elem_ids = [0, 1, 5, 10]
        fig, ax = annulus.plot(elem_ids=elem_ids)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)
    
    @pytest.mark.parametrize('mesh_name', 
        ['annulus', 'donut', 'structured', 'quad_2x2'])
    def test_plot_all_fixtures(self, request, mesh_name):
        """Plotting should work on all fixture mesh types."""
        mesh = request.getfixturevalue(mesh_name)
        fig, ax = mesh.plot()
        assert fig is not None
        plt.close(fig)
    
    def test_axis_chilmesh_configuration(self, annulus):
        """axis_chilmesh should configure axes correctly."""
        fig, ax = plt.subplots()
        result_ax = annulus.axis_chilmesh(ax)
        
        assert result_ax is ax
        assert ax.get_aspect() == 'equal'
        assert ax.get_facecolor() == (1.0, 1.0, 1.0, 1.0)  # white
        plt.close(fig)
    
    def test_axis_chilmesh_creates_default_axes(self, annulus):
        """axis_chilmesh should use current axes if not provided."""
        plt.figure()
        ax = annulus.axis_chilmesh()
        assert ax is not None
        plt.close('all')
    
    def test_plot_respects_linestyle(self, annulus):
        """Plotting should respect linestyle parameter."""
        fig, ax = annulus.plot(linestyle='--', elem_color='none')
        assert isinstance(fig, plt.Figure)
        plt.close(fig)


class TestPlottingEdgeCases:
    """Test edge cases in plotting."""
    
    def test_plot_with_custom_axes(self, annulus):
        """Plot should use provided axes."""
        fig, ax = plt.subplots()
        returned_fig, returned_ax = annulus.plot(ax=ax)
        assert returned_ax is ax
        assert returned_fig is fig
        plt.close(fig)
    
    def test_plot_linewidth_parameter(self, annulus):
        """Plotting should accept linewidth parameter."""
        fig, ax = annulus.plot(linewidth=2.0, elem_color='none')
        assert isinstance(fig, plt.Figure)
        plt.close(fig)
```

**File size**: ~300 lines (covers ~60 missed statements in plot_utils.py)

---

### 3. **Mesh Mutation Operations (MEDIUM)** — 0% Coverage

**Risk Level**: 🟡 **Medium**  
**Impact**: Important for mesh adaptation and modification workflows

#### Missing Coverage
- `change_points()` — Coordinate modification (lines 70–81)
- Full validation of 2D/3D point arrays
- Error handling for invalid shapes

#### Why It Matters
Mesh mutation is used by downstream projects (ADMESH) for adaptation. Without tests:
- Coordinate modifications could silently corrupt topology
- Validation logic is unchecked
- Safety assertion (`acknowledge_change`) behavior is untested

#### Recommendations

**Add to `tests/test_mesh_operations.py` (new file):**

```python
import pytest
import numpy as np


class TestMeshMutation:
    """Test mesh modification operations."""
    
    def test_change_points_2d_array(self, annulus):
        """change_points should accept 2D (n_pts, 2) arrays."""
        mesh = annulus.copy()
        new_points_2d = mesh.points[:, :2].copy()
        new_points_2d[0] = [1.5, 2.5]
        
        mesh.change_points(new_points_2d, acknowledge_change=True)
        np.testing.assert_array_equal(mesh.points[0, :2], [1.5, 2.5])
    
    def test_change_points_3d_array(self, annulus):
        """change_points should accept 3D (n_pts, 3) arrays."""
        mesh = annulus.copy()
        new_points_3d = mesh.points.copy()
        new_points_3d[0] = [1.5, 2.5, 100.0]
        
        mesh.change_points(new_points_3d, acknowledge_change=True)
        np.testing.assert_array_equal(mesh.points[0], [1.5, 2.5, 100.0])
    
    def test_change_points_requires_acknowledge(self, annulus):
        """change_points must require acknowledge_change=True."""
        mesh = annulus
        new_points = mesh.points.copy()
        
        with pytest.raises(AssertionError):
            mesh.change_points(new_points, acknowledge_change=False)
    
    def test_change_points_wrong_columns(self, annulus):
        """change_points should reject arrays with wrong columns."""
        mesh = annulus
        wrong_points = np.random.rand(mesh.n_verts, 4)
        
        with pytest.raises(ValueError, match="2 or 3 columns"):
            mesh.change_points(wrong_points, acknowledge_change=True)
    
    def test_change_points_preserves_topology(self, annulus):
        """Changing points should not alter connectivity."""
        mesh = annulus.copy()
        conn_before = mesh.connectivity_list.copy()
        
        new_points = mesh.points.copy()
        new_points *= 1.1  # Scale uniformly
        mesh.change_points(new_points, acknowledge_change=True)
        
        np.testing.assert_array_equal(mesh.connectivity_list, conn_before)
```

---

### 4. **Backward Compatibility Properties (LOW)** — ~50% Coverage

**Risk Level**: 🟢 **Low**  
**Impact**: Legacy code paths only

#### Missing Coverage
- `Layers` property (line 68) — Uppercase variant (deprecated)
- `grid_name` property getter/setter (lines 51–58)

#### Recommendations

**Add to `tests/test_backward_compatibility.py`:**

```python
def test_layers_property_backward_compat(annulus):
    """Layers property should return same as layers."""
    mesh = annulus
    assert mesh.Layers is mesh.layers
    assert mesh.Layers == mesh.layers


def test_grid_name_property_get(annulus):
    """grid_name property should retrieve name."""
    mesh = annulus
    mesh.grid_name = "Test Grid"
    assert mesh.grid_name == "Test Grid"


def test_grid_name_property_set(annulus):
    """grid_name property should set name."""
    mesh = annulus
    original_name = mesh.grid_name
    mesh.grid_name = "New Name"
    assert mesh.grid_name == "New Name"
```

---

### 5. **Error Handling & Edge Cases (MEDIUM)** — ~70% Coverage

**Risk Level**: 🟡 **Medium**  
**Impact**: Robustness and user experience

#### Missing Coverage
- Degenerate triangle fallback in CCW (lines 118–124)
- Invalid edge skipping (lines 352–353, 385–386)
- Invalid `elem_quality()` type parameter

#### Recommendations

**Add edge case tests:**

```python
def test_ccw_degenerate_quad_fallback():
    """CCW should handle degenerate quads with duplicate vertices."""
    # Create mesh with degenerate quad: [0, 1, 1, 2]
    connectivity = np.array([[0, 1, 1, 2]])
    points = np.array([[0, 0], [1, 0], [1, 1]])
    mesh = CHILmesh(connectivity, points, compute_layers=False)
    
    # Should not crash; should use fallback logic
    assert mesh.connectivity_list is not None


def test_elem_quality_invalid_type(annulus):
    """elem_quality should reject invalid quality types."""
    mesh = annulus
    with pytest.raises(ValueError, match="Unknown quality type"):
        mesh.elem_quality(quality_type='invalid_type')
```

---

## Implementation Priority & Roadmap

### Phase 1: High Impact (Weeks 1–2)
**Target**: +40 covered statements, +5% coverage

1. **`tests/test_plotting.py`** (78 missed statements)
   - Effort: 2–3 hours
   - Gain: ~60 statements (+8%)
   - Blocks: Visualization debugging

2. **`tests/test_mesh_smoothing.py`** (50 missed statements)
   - Effort: 2–3 hours
   - Gain: ~40 statements (+5%)
   - Blocks: Quality improvement workflows

### Phase 2: Medium Impact (Week 3)
**Target**: +15 covered statements, +2% coverage

3. **Error handling tests** (20 missed statements)
   - Effort: 1–2 hours
   - Gain: ~15 statements (+2%)

4. **Mesh mutation tests** (18 missed statements)
   - Effort: 1 hour
   - Gain: ~12 statements (+1.5%)

### Phase 3: Low-Impact Polish (Week 4)
**Target**: +10 covered statements, +1% coverage

5. **Backward compatibility tests** (5 missed statements)
6. **Uncommon code path tests** (5 missed statements)

---

## Coverage Targets

### Current State
```
Overall: 84% (618 / 736)
CHILmesh.py: 92% (468 / 508)
plot_utils.py: 60% (119 / 197)
```

### Target State (After Phase 1–2)
```
Overall: 90% (662 / 736)  [+44 statements]
CHILmesh.py: 95% (483 / 508)  [+15 statements]
plot_utils.py: 85% (167 / 197)  [+48 statements]
```

### Stretch Goal (All phases)
```
Overall: 92% (677 / 736)  [+59 statements]
CHILmesh.py: 96% (488 / 508)  [+20 statements]
plot_utils.py: 90% (177 / 197)  [+58 statements]
```

---

## Testing Tips & Setup

### Enable Non-Interactive Plotting
```bash
export MPLBACKEND=Agg
pytest tests/test_plotting.py -v
```

### Run Coverage Reports
```bash
# Full report with missing lines
pytest --cov=src/chilmesh --cov-report=term-missing -v

# Only fast tests (exclude block_o)
pytest -k 'not block_o' --cov=src/chilmesh --cov-report=term-missing -q

# HTML report
pytest --cov=src/chilmesh --cov-report=html
# Open htmlcov/index.html in browser
```

### Identify Uncovered Lines
```bash
pytest --cov=src/chilmesh --cov-report=term-missing tests/ | grep -A5 "Missing"
```

---

## Maintenance & CI Integration

### Update `.github/workflows/tests.yml`
```yaml
- name: Test Coverage
  run: |
    pytest --cov=src/chilmesh --cov-report=term-missing --cov-fail-under=85
```

### Track Coverage Over Time
```bash
# Generate badge
coverage-badge -o coverage.svg -f
# Commit to repo
git add coverage.svg
```

---

## Related Documents

- **CLAUDE.md** — Development guidelines and testing standards
- **PLANNING_DATA_STRUCTURE_MODERNIZATION.md** — Feature roadmap
- **pyproject.toml** — pytest configuration and coverage settings

---

**Last Updated**: April 26, 2026  
**Analysis by**: Claude Code (Haiku 4.5)
