# CHILmesh Code Structure Audit

**Issue:** #103 (Code structure audit - overly faithful MATLAB port)  
**Date:** 2026-05-15  
**Finding:** Monolithic CHILmesh.py (80 KB) contains 48 methods + 1 function in single file.

---

## Current Structure

### File: `src/chilmesh/CHILmesh.py` (80 KB, 2,500+ LOC)

**Single class:** `CHILmesh` with 48 methods

**Distribution by concern:**
```
│ Initialization & Topology (10 methods)
├─ __init__, _initialize_mesh, _ensure_ccw_orientation
├─ _build_adjacencies, _validate_adjacencies, _identify_edges
├─ _build_elem2edge, _build_vert2edge, _build_vert2elem, _build_edge2elem

│ I/O Operations (5 methods)
├─ read_from_fort14, write_to_fort14, read_from_2dm
├─ from_admesh_domain, admesh_metadata

│ Topology Queries (8 methods)
├─ boundary_edges, boundary_node_indices, get_vertex_edges
├─ get_vertex_elements, get_layer, elements_in_layer
├─ edge2vert, elem2edge, edge2elem, etc.

│ Skeletonization (3 methods)
├─ _skeletonize, _mesh_layers, get_layer (partial)

│ Quality Metrics (2 methods)
├─ elem_quality, interior_angles

│ Smoothing & Optimization (15 methods)
├─ smooth_mesh, angle_based_smoother, direct_smoother
├─ _tri_stiffness_assembly, _quad_stiffness_assembly, _mixed_stiffness_assembly
├─ advancing_front_boundary_edges, add_advancing_front_element, pinch_points
├─ _ordered_vertex_ring, remove_boundary_loop, _detect_element_types

│ Properties & Utilities (5 methods)
├─ grid_name (property), Layers (property), copy
├─ signed_area, _elem_type, change_points
```

---

## Assessment

### Cohesion: GOOD
Methods are logically related (all manipulate same mesh state). Decomposing won't improve logic clarity.

### User API: PROBLEMATIC
**Current:** `from chilmesh import CHILmesh`  
**Problem:** User must import full CHILmesh to use any submodule. No way to access only smoothing, only I/O, etc.

**Example:** User wants to write a mesh to fort.14:
```python
# Current: Must import full CHILmesh
from chilmesh import CHILmesh
m = CHILmesh(...)
m.write_to_fort14("output.14")

# Desired: Could import I/O only
from chilmesh.io import write_fort14
write_fort14(mesh_data, "output.14")
```

### Pythonicity: MATLAB-LIKE
Method naming and organization follow MATLAB QuADMesh+ structure closely (faithful port). Idiomatic Python would separate concerns into modules.

---

## Refactoring Options

### Option A: Leave As-Is (Simplest)
- **Pros:** No breaking changes; fast execution
- **Cons:** Monolithic; hard to navigate; users must understand full API
- **Recommendation:** v1.0 safe choice; defer to Phase 5 after stabilization

### Option B: Extract Submodules (Medium Effort)
Decompose into:
```
chilmesh/
├─ __init__.py              (re-export public API)
├─ mesh.py                  (core CHILmesh class, 80 KB → 40 KB)
├─ io.py                    (read_from_fort14, write_to_fort14, read_from_2dm)
├─ skeletonization.py       (_skeletonize, _mesh_layers, get_layer)
├─ smoothing.py             (smooth_mesh, *_smoother, stiffness assemblies)
├─ quality.py               (elem_quality, interior_angles)
└─ utils.py                 (signed_area, _elem_type, copy, etc.)
```

**Pros:**
- Users can `from chilmesh.io import read_fort14` (specific imports)
- Code organization mirrors concern separation
- Easier to find methods
- Reduces cognitive load (small files)

**Cons:**
- Requires careful design of module interdependencies
- May expose internal implementation details
- Breaks some internal method privacy (_skeletonize would become public)

**Backward compatibility:** Can re-export all public methods from `__init__.py`:
```python
# chilmesh/__init__.py
from .mesh import CHILmesh
from .io import read_fort14, write_fort14
# ... rest of exports

# Old code still works:
from chilmesh import CHILmesh
# New code also works:
from chilmesh.io import write_fort14
```

### Option C: Mixin Pattern (Light-Touch)
Keep CHILmesh as main class, extract concerns as mixins:
```python
class IOmixin:
    def read_from_fort14(self): ...
    def write_to_fort14(self): ...

class SmootherMixin:
    def smooth_mesh(self): ...

class CHILmesh(IOMixin, SmootherMixin, TopologyMixin):
    # Core initialization + properties
```

**Pros:**
- Non-breaking (no import changes)
- Modest code organization improvement
- No module reorganization needed

**Cons:**
- Mixins still live in single file (doesn't reduce monolithism)
- Mixin composition less discoverable than modules
- Adds metaclass complexity

---

## Current Issues Caused by Structure

1. **Navigation:** 80 KB file; finding a method requires grep or IDE search
2. **Discoverability:** Users must read full API docs; can't discover submodules
3. **Testing:** Test fixtures are parametrized across ALL 48 methods (slow)
4. **Maintenance:** Any change in one method risks breaking unrelated functionality (high coupling)

---

## Recommendation

**v1.0:** Keep Option A (as-is). Rationale:
- Functionality is stable; focus is release, not refactoring
- Monolithic structure is working (fast, tested, minimal bugs)
- Constitution Principle I (library-first) is met — CHILmesh IS a standalone library
- Constitution says: "no YAGNI abstractions" — decomposition not yet needed

**Post-v1.0 (Phase 5):** Revisit Option B or C if:
- Users complain about discoverability
- Smoothing module grows (planned mutation API)
- I/O module needs to support more formats (GeoTIFF, Shapefile, etc.)

---

## MATLAB Port Fidelity

**Assessment:** Structure faithfully mirrors MATLAB QuADMesh+.  
**Pythonicity:** Medium. Class-based design is appropriate; method naming is OK.  
**Criticality:** Low. No bugs or correctness issues caused by monolithic structure.

---

## Related

- Issue #103: Code structure audit (this audit)
- Constitution Principle I: Library-first architecture (satisfied)
- Issue #101: Python 3.8+ compatibility (different concern)
- Issue #94-#92: Phase 5 design (future refactoring trigger point)
