# Downstream Project Integration Guide

**Version:** 2.0 (revised for CHILmesh v1.0.0)
**CHILmesh Version:** 1.0.0+
**Target Projects:** MADMESHR, ADMESH, ADMESH-Domains

---

## Overview

Guide for developers of downstream research projects integrating with CHILmesh v1.0.0+.

**TL;DR:** Existing code keeps working — `CHILmesh` is still importable. The new `Mesh` alias is the v1.0.0 preferred idiom; adopt it when convenient. Optional C++ backend gives 46× speedup with bit-identical output.

---

## What's New in v1.0.0?

### Preferred Idiom (recommended for new code)
```python
from chilmesh import Mesh

mesh = Mesh.read_from_fort14("ocean.14")
```

`Mesh` is an alias for `CHILmesh` — no behaviour change. Existing imports continue to work through the v1.x series. A `DeprecationWarning` lands in v1.1.0; removal no earlier than v2.0.0.

### Optional C++ Acceleration
Heavy initialization workloads can opt into the bundled C++ half-edge extension. Output is bit-identical to Python (verified by `tests/test_backend_equivalence.py`, 36 parametrized cases).

```python
import chilmesh
chilmesh.backend_info()
# {'available': ['cpp', 'rust', 'python'], 'selected': 'cpp', ...}
```

On WNAT_Hagen (52,774 verts · 98,365 elements):
- Python full init: 3.21 s
- C++ full init: **0.069 s** (46×)
- C++ skeletonization in isolation: **0.033 s vs Python 2.20 s (66×)**

Force a specific backend with the `CHILMESH_BACKEND` environment variable
(`python`, `cpp`, or `rust`).

### New Public Surface (cumulative since v0.4.1)
- `MeshAdapterForMADMESHR`, `MeshAdapterForADMESH`, `MeshAdapterForADMESHDomains` re-exported at package root.
- `CHILmesh.submesh(elem_ids)` — public sub-mesh factory (#138).
- `CHILmesh.rebuild_adjacencies()` / `invalidate_adjacencies()` — for mid-sweep mutation flows (#143).
- `CHILmesh.ccw_edges_around_vert(vert_id)` — counter-clockwise edge walk (#133).
- `compute_adjacencies` kwarg on the constructor and `read_from_fort14` / `read_from_2dm` / `from_admesh_domain` (#134).
- `chilmesh.backend_info()` — runtime backend introspection.

### Guarantees
- **Public API frozen under semver.** Breaking changes require a v2.x bump.
- `CHILmesh` import preserved through v1.x; DeprecationWarning in v1.1, earliest removal v2.0.
- C++ / Rust backends produce bit-identical output to Python on the official fixture set.
- Python `>=3.10` (3.8 / 3.9 dropped — EOL).

---

## Do I Need to Change My Code?

**No.** Existing code continues to work without modification.

**Should I update?** Yes, eventually. Benefits: better performance (1.5x improvement), clearer code (bridge adapters), better error messages, new quality metrics.

**Migration timeline:**
- **Now**: Update and test (should be drop-in replacement)
- **Later**: Optionally adopt bridge adapters for clearer code
- **No rush**: Take time to validate tests pass first

---

## Quick Start by Project

### MADMESHR Quick Start

#### Before (v0.1.1)
```python
from chilmesh import CHILmesh

mesh = CHILmesh.read_from_fort14("domain.fort.14")

# Finding neighbors manually (verbose)
def get_neighbors(mesh, elem_id):
    neighbors = set()
    edges = mesh.elem2edge(elem_id)
    for edge_id in edges:
        e1, e2 = mesh.edge2elem(edge_id)
        if e1 == elem_id and e2 >= 0:
            neighbors.add(e2)
        elif e2 == elem_id and e1 >= 0:
            neighbors.add(e1)
    return neighbors

neighbors = get_neighbors(mesh, 0)
```

#### After (v0.2.0+) - Still Works
```python
from chilmesh import CHILmesh

mesh = CHILmesh.read_from_fort14("domain.fort.14")

# Old code still works - no changes needed
def get_neighbors(mesh, elem_id):
    neighbors = set()
    edges = mesh.elem2edge(elem_id)
    for edge_id in edges:
        e1, e2 = mesh.edge2elem(edge_id)
        if e1 == elem_id and e2 >= 0:
            neighbors.add(e2)
        elif e2 == elem_id and e1 >= 0:
            neighbors.add(e1)
    return neighbors

neighbors = get_neighbors(mesh, 0)
```

#### Recommended Update (Optional)
```python
from chilmesh import CHILmesh
from chilmesh.bridge import MeshAdapterForMADMESHR

mesh = CHILmesh.read_from_fort14("domain.fort.14")
adapter = MeshAdapterForMADMESHR(mesh)

# Much clearer - adapter handles the logic
neighbors = adapter.get_element_neighbors(0)
```

### ADMESH Quick Start

#### Before (v0.1.1)
```python
from chilmesh import CHILmesh
import numpy as np

mesh = CHILmesh.read_from_fort14("domain.fort.14")

# Manual quality assessment
quality, angles, stats = mesh.elem_quality()
poor_count = np.sum(quality < 0.3)
print(f"Mean quality: {stats['mean']:.3f}")
print(f"Poor elements: {poor_count}")
```

#### After (v0.2.0+) - Still Works
```python
from chilmesh import CHILmesh
import numpy as np

mesh = CHILmesh.read_from_fort14("domain.fort.14")

# Same code works exactly as before
quality, angles, stats = mesh.elem_quality()
poor_count = np.sum(quality < 0.3)
print(f"Mean quality: {stats['mean']:.3f}")
print(f"Poor elements: {poor_count}")
```

#### Recommended Update (Optional)
```python
from chilmesh import CHILmesh
from chilmesh.bridge import MeshAdapterForADMESH

mesh = CHILmesh.read_from_fort14("domain.fort.14")
adapter = MeshAdapterForADMESH(mesh)

# Clearer intent - dedicated quality report method
report = adapter.get_mesh_quality_report()
print(f"Mean quality: {report['mean']:.3f}")
print(f"Poor elements: {report['poor_count']}")
```

### ADMESH-Domains Quick Start

#### Before (v0.1.1)
```python
from chilmesh import CHILmesh

mesh = CHILmesh.read_from_fort14("domain.fort.14")

# Manual boundary extraction
boundary_edges = mesh.boundary_edges()
boundary_verts = set()
for edge_id in boundary_edges:
    v1, v2 = mesh.edge2vert(edge_id)
    boundary_verts.add(v1)
    boundary_verts.add(v2)

print(f"Boundary vertices: {len(boundary_verts)}")
```

#### After (v0.2.0+) - Still Works
```python
from chilmesh import CHILmesh

mesh = CHILmesh.read_from_fort14("domain.fort.14")

# Same code works exactly as before
boundary_edges = mesh.boundary_edges()
boundary_verts = set()
for edge_id in boundary_edges:
    v1, v2 = mesh.edge2vert(edge_id)
    boundary_verts.add(v1)
    boundary_verts.add(v2)

print(f"Boundary vertices: {len(boundary_verts)}")
```

#### Recommended Update (Optional)
```python
from chilmesh import CHILmesh
from chilmesh.bridge import MeshAdapterForADMESHDomains

mesh = CHILmesh.read_from_fort14("domain.fort.14")
adapter = MeshAdapterForADMESHDomains(mesh)

# Clearer intent - dedicated domain boundary method
boundaries = adapter.get_domain_boundaries()
boundary_verts = boundaries[0]
print(f"Boundary vertices: {len(boundary_verts)}")
```

---

## Migration Strategy (Recommended)

### Phase 1: Update CHILmesh (1 hour)
```bash
pip install --upgrade 'chilmesh>=1.0.0'
```

### Phase 2: Test Compatibility (30 minutes)
Run existing test suite. Everything should pass without code changes.

**If tests fail:**
- Check error messages
- Report issue at: https://github.com/domattioli/CHILmesh/issues
- Include: CHILmesh version, code snippet, error traceback

### Phase 3: Adopt Bridge Adapters (Optional, 2-4 hours)
Gradually replace verbose patterns with adapter methods:
1. Find mesh interactions in your code
2. Check if adapter has convenience method for your use case
3. Refactor to use adapter (code becomes clearer)
4. Run tests to verify behavior unchanged

### Phase 4: Performance Profiling (Optional, 1-2 hours)
```bash
python -m cProfile -s cumtime your_script.py > profile.txt
```
Compare before/after upgrade. Bridge adapters often run faster than manual patterns.

---

## Common Patterns & Bridge Adapters

### Pattern 1: Finding Element Neighbors

**Old (v0.1.1):**
```python
def get_neighbors(mesh, elem_id):
    neighbors = set()
    edges = mesh.elem2edge(elem_id)
    for edge_id in edges:
        e1, e2 = mesh.edge2elem(edge_id)
        if e1 == elem_id and e2 >= 0:
            neighbors.add(e2)
        elif e2 == elem_id and e1 >= 0:
            neighbors.add(e1)
    return neighbors
```

**New (v0.2.0+ with adapter):**
```python
from chilmesh.bridge import MeshAdapterForMADMESHR
adapter = MeshAdapterForMADMESHR(mesh)
neighbors = adapter.get_element_neighbors(elem_id)
```

### Pattern 2: Quality Assessment

**Old (v0.1.1):**
```python
quality, angles, stats = mesh.elem_quality()
poor_count = np.sum(quality < 0.3)
mean_qual = stats['mean']
```

**New (v0.2.0+ with adapter):**
```python
from chilmesh.bridge import MeshAdapterForADMESH
adapter = MeshAdapterForADMESH(mesh)
report = adapter.get_mesh_quality_report()
poor_count = report['poor_count']
mean_qual = report['mean']
```

### Pattern 3: Boundary Extraction

**Old (v0.1.1):**
```python
boundary_edges = mesh.boundary_edges()
boundary_verts = set()
for edge_id in boundary_edges:
    v1, v2 = mesh.edge2vert(edge_id)
    boundary_verts.add(v1)
    boundary_verts.add(v2)
```

**New (v0.2.0+ with adapter):**
```python
from chilmesh.bridge import MeshAdapterForADMESHDomains
adapter = MeshAdapterForADMESHDomains(mesh)
boundaries = adapter.get_domain_boundaries()
boundary_verts = boundaries[0]
```

---

## Troubleshooting

### Error: "AttributeError: CHILmesh has no attribute 'get_vertex_edges'"

**Cause:** Using CHILmesh < 1.0.0

**Fix:**
```bash
pip install --upgrade 'chilmesh>=1.0.0'
```

### Error: "ModuleNotFoundError: No module named 'chilmesh.bridge'"

**Cause:** Using CHILmesh < 1.0.0, or bridge module not found

**Fix:**
```bash
pip install --upgrade 'chilmesh>=1.0.0'
python -c "from chilmesh.bridge import MeshAdapterForMADMESHR; print('OK')"
```

### Error: "KeyError: 'Vert2Edge'"

**Cause:** Accessing internal adjacency structures directly

**Old (broken in 0.2.0):**
```python
# DON'T do this - it breaks in 0.2.0
edges = mesh.adjacencies['Vert2Edge'][v]
```

**New (recommended):**
```python
# Use the public API instead
edges = mesh.get_vertex_edges(v)
```

### Performance Issues

**Problem:** Code slower after upgrading? (Rare — should be faster.)

1. **Verify 0.2.0+:**
   ```python
   import chilmesh
   print(chilmesh.__version__)  # Should be 0.2.0 or later
   ```

2. **Profile:**
   ```bash
   python -m cProfile -s cumtime your_script.py > profile.txt
   ```

3. **Report** with CHILmesh version, minimal code, timing before/after, mesh size.

---

## API Reference by Project

### MADMESHR: Mesh Adaptation Research

**Relevant CAI Methods:**
```python
# Direct methods
elem_neighbors = mesh.elem2edge(elem_id)
quality, angles, stats = mesh.elem_quality(elem_ids)
elements_at_edges = mesh.edge2elem(edge_id)

# New bridge methods (v0.2.0+)
neighbors = adapter.get_element_neighbors(elem_id)
quality_neighborhood = adapter.get_element_quality_neighborhood(elem_id)
refinement_region = adapter.get_refinement_region(seed_elems)
```

### ADMESH: Adaptive Mesh Refinement

**Relevant CAI Methods:**
```python
# Direct methods
quality, angles, stats = mesh.elem_quality()
interior_angles = mesh.interior_angles()

# New bridge methods (v0.2.0+)
report = adapter.get_mesh_quality_report()
angle_summary = adapter.get_element_angles_summary()
```

### ADMESH-Domains: Domain Decomposition

**Relevant CAI Methods:**
```python
# Direct methods
boundary_edges = mesh.boundary_edges()
vertices = mesh.points
connectivity = mesh.connectivity_list

# New bridge methods (v0.2.0+)
boundaries = adapter.get_domain_boundaries()
connectivity_info = adapter.get_mesh_connectivity_info()
```

---

## Getting Help

### Questions About Migration
1. Check this guide
2. Review `docs/CHILmesh_Access_Interface.md`
3. Look at `tests/test_bridge_adapters.py`

### Bug Reports
Include:
- CHILmesh version: `python -c "import chilmesh; print(chilmesh.__version__)"`
- Python version: `python --version`
- Minimal reproducible example
- Full error traceback
- Downstream project name (MADMESHR/ADMESH/ADMESH-Domains)

Report at: https://github.com/domattioli/CHILmesh/issues

---

## Version Compatibility

| CHILmesh | MADMESHR | ADMESH | ADMESH-Domains | Notes |
|----------|----------|--------|----------------|-------|
| 0.1.1    | ✅ Legacy | ✅ Legacy | ✅ Legacy | Old version, still works |
| 0.2.0    | ✅ Recommended | ✅ Recommended | ✅ Recommended | Current, use this |
| 0.2.x    | ✅ Recommended | ✅ Recommended | ✅ Recommended | Bug fixes, recommended |
| 1.0.0    | TBD | TBD | TBD | Future, full stability |

---

## FAQ

**Q: Will my code break when upgrading?**
No. CHILmesh 0.2.0 is backward compatible. All existing APIs work unchanged.

**Q: Should I use the bridge adapters?**
Optional but recommended. Code becomes clearer and often runs faster. Adopt at your own pace.

**Q: What if I need old behavior?**
Stick with v0.1.1. But upgrade recommended for performance benefits.

**Q: Can I mix old and new code?**
Yes. Both old patterns and bridge adapters work in same codebase during migration.

**Q: Is bridge adapter API stable?**
Yes. Adapter methods are part of CAI, guaranteed stable through v1.0.

**Q: How much code do I need to change?**
None required. Migration completely optional.

---

## Example Integration Projects

See `examples/` directory for complete working examples:
- `madmeshr_refinement.py` - MADMESHR mesh adaptation workflow
- `admesh_quality.py` - ADMESH quality assessment
- `admesh_domains_setup.py` - ADMESH-Domains domain initialization

---

## Support Timeline

**CHILmesh 0.2.0 Support:**
- Release date: April 2026
- Bug fix support: Until 0.3.0 release
- Security support: Until 1.0.0 release

**Migration Support:**
- Expect 1-4 hours updating your project
- Most time spent testing, not code changes

---

**Last Updated:** 2026-04-27
**Guide Version:** 1.0

For latest information, visit: https://github.com/domattioli/CHILmesh
