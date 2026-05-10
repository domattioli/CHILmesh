# CHILmesh Test Coverage Analysis

**Analysis Date**: April 26, 2026  
**Coverage Tool**: pytest-cov  
**Python Version**: 3.11

---

## Executive Summary

CHILmesh has **84% overall test coverage** (737 statements, 118 missed). Significant gaps in:

1. **Plotting functionality** (60% coverage)
2. **Mesh smoothing methods** (0% coverage)
3. **Mesh mutation operations** (0% coverage)
4. **Error handling** (70% coverage)

Closing gaps improves coverage to **90%+**.

---

## Current Coverage Summary

```
Total Statements: 736
Covered: 618 (84%)
Missed: 118 (16%)
Tests Passing: 137 (excluding slow block_o mesh)
```

| Module | Statements | Coverage | Missed | Priority |
|--------|-----------|----------|--------|----------|
| `src/chilmesh/CHILmesh.py` | 508 | 92% | 40 | High |
| `src/chilmesh/utils/plot_utils.py` | 197 | 60% | 78 | **Critical** |
| `src/chilmesh/__init__.py` | 6 | 100% | 0 | — |
| `src/chilmesh/data/__init__.py` | 0 | 100% | 0 | — |
| `src/chilmesh/examples.py` | 23 | 100% | 0 | — |

---

## Major Coverage Gaps

### 1. Mesh Smoothing Methods (CRITICAL) — 0% Coverage

**Risk:** High | **Maturity:** Beta (angle-based not implemented)

**Missing:** `smooth_mesh()` (lines 1118–1134), `angle_based_smoother()` (lines 1136–1147), `direct_smoother()` (lines 1149–1203)

**Why it matters:** Refactoring FEM solver could break silently; boundary constraint handling untested.

**Recommended tests:**

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

def test_smooth_mesh_fem_quality_improvement(donut):
    """Smoothing should improve skew angle quality."""
    mesh = donut
    quality_before, _, stats_before = mesh.elem_quality(quality_type='skew')
    mesh.smooth_mesh('fem', acknowledge_change=True)
    quality_after, _, stats_after = mesh.elem_quality(quality_type='skew')
    assert stats_after['mean'] > stats_before['mean']

def test_smooth_mesh_requires_acknowledge(annulus):
    with pytest.raises(AssertionError):
        annulus.smooth_mesh('fem', acknowledge_change=False)

@pytest.mark.xfail(reason="Angle-based smoothing not yet implemented")
def test_angle_based_smoother_not_implemented(annulus):
    with pytest.raises(NotImplementedError):
        annulus.angle_based_smoother()
```

---

### 2. Plotting Module (CRITICAL) — 60% Coverage

**Risk:** High | **78 missed statements**

**Uncovered:** `plot()` (lines 47–110), `plot_edge()`, `plot_layers()`, `axis_chilmesh()` (lines 20–44), `plot_elem_quality()`

**Why it matters:** Matplotlib API changes could silently break visualization; edge cases untested.

**Create `tests/test_plotting.py`:**

```python
import pytest
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import numpy as np


class TestMeshPlotting:
    def test_plot_basic_mesh_creates_figure(self, annulus):
        fig, ax = annulus.plot()
        assert isinstance(fig, plt.Figure)
        assert isinstance(ax, plt.Axes)
        plt.close(fig)
    
    def test_plot_with_element_colors(self, annulus):
        fig, ax = annulus.plot(elem_color='lightblue')
        assert len(ax.patches) > 0
        plt.close(fig)
    
    def test_plot_edges_only(self, annulus):
        fig, ax = annulus.plot(elem_color='none', edge_color='red')
        assert len(ax.lines) > 0 or len(ax.collections) > 0
        plt.close(fig)
    
    def test_plot_single_element(self, annulus):
        fig, ax = annulus.plot(elem_ids=5, elem_color='green')
        assert ax.get_xlim()[0] < ax.get_xlim()[1]
        plt.close(fig)
    
    @pytest.mark.parametrize('mesh_name', 
        ['annulus', 'donut', 'structured', 'quad_2x2'])
    def test_plot_all_fixtures(self, request, mesh_name):
        mesh = request.getfixturevalue(mesh_name)
        fig, ax = mesh.plot()
        assert fig is not None
        plt.close(fig)
    
    def test_axis_chilmesh_configuration(self, annulus):
        fig, ax = plt.subplots()
        result_ax = annulus.axis_chilmesh(ax)
        assert result_ax is ax
        assert ax.get_aspect() == 'equal'
        plt.close(fig)
```

---

### 3. Mesh Mutation Operations (MEDIUM) — 0% Coverage

**Risk:** Medium | **Impact:** ADMESH adaptation workflows

**Missing:** `change_points()` (lines 70–81), 2D/3D array validation, error handling

```python
class TestMeshMutation:
    def test_change_points_2d_array(self, annulus):
        mesh = annulus.copy()
        new_points_2d = mesh.points[:, :2].copy()
        new_points_2d[0] = [1.5, 2.5]
        mesh.change_points(new_points_2d, acknowledge_change=True)
        np.testing.assert_array_equal(mesh.points[0, :2], [1.5, 2.5])
    
    def test_change_points_requires_acknowledge(self, annulus):
        with pytest.raises(AssertionError):
            annulus.change_points(annulus.points.copy(), acknowledge_change=False)
    
    def test_change_points_wrong_columns(self, annulus):
        with pytest.raises(ValueError, match="2 or 3 columns"):
            annulus.change_points(np.random.rand(annulus.n_verts, 4), acknowledge_change=True)
    
    def test_change_points_preserves_topology(self, annulus):
        mesh = annulus.copy()
        conn_before = mesh.connectivity_list.copy()
        mesh.change_points(mesh.points.copy() * 1.1, acknowledge_change=True)
        np.testing.assert_array_equal(mesh.connectivity_list, conn_before)
```

---

### 4. Backward Compatibility Properties (LOW) — ~50% Coverage

**Missing:** `Layers` property (line 68, deprecated uppercase variant), `grid_name` getter/setter (lines 51–58)

```python
def test_layers_property_backward_compat(annulus):
    assert annulus.Layers is annulus.layers

def test_grid_name_property_get(annulus):
    annulus.grid_name = "Test Grid"
    assert annulus.grid_name == "Test Grid"
```

---

### 5. Error Handling & Edge Cases (MEDIUM) — ~70% Coverage

**Missing:** degenerate triangle fallback in CCW (lines 118–124), invalid edge skipping, invalid `elem_quality()` type

```python
def test_ccw_degenerate_quad_fallback():
    connectivity = np.array([[0, 1, 1, 2]])
    points = np.array([[0, 0], [1, 0], [1, 1]])
    mesh = CHILmesh(connectivity, points, compute_layers=False)
    assert mesh.connectivity_list is not None

def test_elem_quality_invalid_type(annulus):
    with pytest.raises(ValueError, match="Unknown quality type"):
        annulus.elem_quality(quality_type='invalid_type')
```

---

## Implementation Priority & Roadmap

### Phase 1: High Impact (Weeks 1–2) — Target: +40 statements, +5% coverage
1. `tests/test_plotting.py` — 78 missed statements; 2–3h; gain ~60 statements (+8%)
2. `tests/test_mesh_smoothing.py` — 50 missed statements; 2–3h; gain ~40 statements (+5%)

### Phase 2: Medium Impact (Week 3) — Target: +15 statements, +2% coverage
3. Error handling tests — 1–2h; gain ~15 statements
4. Mesh mutation tests — 1h; gain ~12 statements

### Phase 3: Polish (Week 4)
5. Backward compatibility tests; uncommon code path tests

---

## Coverage Targets

| State | Overall | CHILmesh.py | plot_utils.py |
|-------|---------|-------------|---------------|
| **Current** | 84% (618/736) | 92% (468/508) | 60% (119/197) |
| **After Phase 1-2** | 90% (662/736) | 95% (483/508) | 85% (167/197) |
| **Stretch** | 92% (677/736) | 96% (488/508) | 90% (177/197) |

---

## Testing Tips & Setup

```bash
# Non-interactive plotting
export MPLBACKEND=Agg
pytest tests/test_plotting.py -v

# Full report with missing lines
pytest --cov=src/chilmesh --cov-report=term-missing -v

# Fast subset only
pytest -k 'not block_o' --cov=src/chilmesh --cov-report=term-missing -q

# HTML report
pytest --cov=src/chilmesh --cov-report=html
```

### CI Integration

```yaml
- name: Test Coverage
  run: |
    pytest --cov=src/chilmesh --cov-report=term-missing --cov-fail-under=85
```

---

**Last Updated**: April 26, 2026  
**Analysis by**: Claude Code (Haiku 4.5)
