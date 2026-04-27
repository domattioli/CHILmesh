# Downstream Project Integration Guide

**Version:** 1.0  
**CHILmesh Version:** 0.2.0+  
**Target Projects:** MADMESHR, ADMESH, ADMESH-Domains  

---

## Overview

This guide is for developers of downstream research projects integrating with CHILmesh 0.2.0+. If you use CHILmesh in MADMESHR, ADMESH, or ADMESH-Domains, this guide explains what changed and how to adapt.

**TL;DR:** Your code probably still works. CHILmesh 0.2.0 is backward compatible. New bridge adapters make common operations clearer.

---

## What Changed in 0.2.0?

### Performance
- **1.5x+ faster** on large meshes (block_o: 30s → ~20s initialization)
- O(1) edge lookups instead of O(n) searches
- More efficient adjacency structures

### Architecture
- Adjacency structures modernized (internal implementation)
- **Public API unchanged** - your code still works
- New bridge adapters provide convenient shortcuts

### API Additions
- New public methods: `get_vertex_edges()`, `get_vertex_elements()`
- New bridge adapters for domain-specific operations
- Stable API contract via CAI (CHILmesh Access Interface)

### Guarantees
- Backward compatible (no breaking changes)
- Stable through v1.0 (method signatures guaranteed)
- Clear path for future evolution

---

## Do I Need to Change My Code?

### Short Answer
**No.** Your existing code will continue to work without modification.

### Should I Update?
**Yes, eventually.** Updating provides:
- Better performance (1.5x improvement on large meshes)
- Clearer code (bridge adapters eliminate verbose patterns)
- Better error messages and validation
- Access to new quality metrics and analysis methods

### Migration Timeline
- **Now**: Update and test (should be drop-in replacement)
- **Later**: Optionally adopt bridge adapters for clearer code
- **No rush**: Take time to validate your tests pass first

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
pip install --upgrade 'chilmesh>=0.2.0'
```

### Phase 2: Test Compatibility (30 minutes)
Run your existing test suite. Everything should pass without code changes.

**If tests fail:**
- Check error messages
- Report issue at: https://github.com/domattioli/CHILmesh/issues
- Include: CHILmesh version, your code snippet, error traceback

### Phase 3: Adopt Bridge Adapters (Optional, 2-4 hours)
Gradually replace verbose patterns with adapter methods:

1. **Find mesh interactions** in your code
2. **Check if adapter has convenience method** for your use case
3. **Refactor to use adapter** (your code becomes clearer)
4. **Run tests** to verify behavior unchanged

### Phase 4: Performance Profiling (Optional, 1-2 hours)
Profile your code to measure improvements:

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

**Cause:** You're using CHILmesh < 0.2.0

**Fix:** Upgrade to 0.2.0+
```bash
pip install --upgrade 'chilmesh>=0.2.0'
```

### Error: "ModuleNotFoundError: No module named 'chilmesh.bridge'"

**Cause:** You're using CHILmesh < 0.2.0, or bridge module not found

**Fix:** Upgrade and verify installation
```bash
pip install --upgrade 'chilmesh>=0.2.0'
python -c "from chilmesh.bridge import MeshAdapterForMADMESHR; print('OK')"
```

### Error: "KeyError: 'Vert2Edge'"

**Cause:** You're accessing internal adjacency structures directly

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

**Problem:** Your code is slower after upgrading?

This is rare (should be faster), but if it happens:

1. **Verify you're using 0.2.0+:**
   ```python
   import chilmesh
   print(chilmesh.__version__)  # Should be 0.2.0 or later
   ```

2. **Profile to find bottlenecks:**
   ```bash
   python -m cProfile -s cumtime your_script.py > profile.txt
   ```

3. **Report with reproducible example:**
   - CHILmesh version
   - Minimal code to reproduce
   - Timing before/after (with mesh size)

---

## API Reference by Project

### MADMESHR: Mesh Adaptation Research

**Key Operations:**
- Finding element neighbors for refinement decisions
- Assessing local mesh quality for adaptation criteria
- Building refinement regions from seed elements

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

**Key Operations:**
- Assessing overall mesh quality
- Identifying poor quality elements for coarsening/refinement
- Analyzing angle distributions for robustness

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

**Key Operations:**
- Extracting domain boundaries
- Understanding mesh connectivity for domain splitting
- Managing multi-domain scenarios

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
1. Check this guide (you might find the answer)
2. Review code examples in `docs/CHILmesh_Access_Interface.md`
3. Look at test examples in `tests/test_bridge_adapters.py`

### Bug Reports
Include:
- CHILmesh version: `python -c "import chilmesh; print(chilmesh.__version__)"`
- Your Python version: `python --version`
- Minimal reproducible example
- Full error traceback
- Your downstream project name (MADMESHR/ADMESH/ADMESH-Domains)

Report at: https://github.com/domattioli/CHILmesh/issues

### Feature Requests
Tell us what would make integration easier:
- What operations do you repeat most?
- Are there bridge adapter methods you'd like to see?
- Performance concerns specific to your use case?

### Performance Concerns
Include:
- CHILmesh version
- Mesh size (n_verts, n_elems)
- Operation that's slow
- Timing before/after upgrade
- Code snippet showing the slow operation

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

### Q: Will my code break when upgrading?
**A:** No. CHILmesh 0.2.0 is backward compatible. All existing APIs work unchanged.

### Q: Should I use the bridge adapters?
**A:** They're optional but recommended. They make code clearer and often run faster. Adopt them at your own pace.

### Q: What if I need the old behavior?
**A:** You can stick with v0.1.1. But we recommend upgrading for performance benefits.

### Q: Can I mix old and new code?
**A:** Yes. You can use both old patterns and bridge adapters in the same codebase during migration.

### Q: Is the bridge adapter API stable?
**A:** Yes. Adapter methods are part of the CAI and guaranteed stable through v1.0.

### Q: How much code do I need to change?
**A:** None required. Migration is completely optional. Adapters are conveniences, not requirements.

---

## Example Integration Projects

See `examples/` directory for complete working examples:

- `madmeshr_refinement.py` - MADMESHR mesh adaptation workflow
- `admesh_quality.py` - ADMESH quality assessment
- `admesh_domains_setup.py` - ADMESH-Domains domain initialization

Copy and adapt these examples for your own projects.

---

## Support Timeline

**CHILmesh 0.2.0 Support:**
- Release date: April 2026
- Bug fix support: Until 0.3.0 release
- Security support: Until 1.0.0 release

**Migration Support:**
- Expect to spend 1-4 hours updating your project
- Most time spent testing, not code changes
- Bridge adapters make future updates easier

---

**Last Updated:** 2026-04-27  
**Guide Version:** 1.0

For the latest information, visit: https://github.com/domattioli/CHILmesh
