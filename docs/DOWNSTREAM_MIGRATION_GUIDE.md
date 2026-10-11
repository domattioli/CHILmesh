# Downstream Project Integration Guide

**Version:** 2.1 (adds v1.4.0 #187 breaking-rename migration)
**CHILmesh Version:** 1.4.0+
**Target Projects:** an earlier internal prototype, ADMESH, Valence

---

## Overview

Most existing code keeps working with CHILmesh v1.0.0+; `CHILmesh` is still importable. The one exception is v1.4.0, which renamed several symbols without aliases. Read [v1.4.0 Breaking Renames](#-v140-breaking-renames-187) before you widen a `chilmesh` pin.

- The `Mesh` alias is the preferred idiom since v1.0.0. Adopt it when convenient.
- The optional C++ backend gives a 46× speedup with bit-identical output.
- This guide is for developers of downstream research projects.

> **⚠️ v1.4.0 breaks existing code.** The #187 lexicon
> ratification renamed several public and semi-public symbols **without compatibility
> aliases**. If your pin allows `chilmesh` 1.4.x and you still call an old name, the code breaks
> on upgrade (silently at install, `AttributeError` at runtime).

---

## ⚠️ v1.4.0 Breaking Renames (#187)

v1.4.0 (2026-07-11) ratified the CHILmesh skeleton/layer lexicon (#187). It renamed the symbols below with no compatibility aliases.

The release was a *minor* bump, so a `chilmesh>=1.x,<2` or caret pin does not exclude it. A consumer breaks on upgrade when both conditions hold:

1. The pin admits 1.4.x.
2. The code still calls a renamed symbol.

The `smooth_mesh` signature change (#251) shipped in the same release. It is a second break vector.

### Rename map

| Old (≤1.3.x) | New (≥1.4.0) | Surface |
|---|---|---|
| `mesh.skeletonize()` | `mesh.peel_layers()` | `CHILmesh` public |
| `mesh._skeletonize()` | `mesh._peel()` | `CHILmesh` private (**removed**; no alias) |
| `_layerize` | `_peel` | `CHILmesh` private |
| `MutableMesh.reskeletonize_local(...)` | `MutableMesh.repeel_local(...)` | `chilmesh.mutations` |
| `MutableMesh.skeletonize_diff(...)` | `MutableMesh.layers_diff(...)` | `chilmesh.mutations` |
| `smooth_mesh(...)` positional splat | `smooth_mesh(method, acknowledge_change=False, *, sdf=None, size_fn=None, **kwargs)` | `CHILmesh` public (#251) |

The names `skeletonize`, `_skeletonize`, `_layerize`, `reskeletonize_local`, and `skeletonize_diff` no longer exist. Grep your tree for those five groups.

### `smooth_mesh` signature (#251)

`method` and `acknowledge_change` stay positional-or-keyword. These calls are safe:

- `mesh.smooth_mesh("laplacian", True)`
- `mesh.smooth_mesh(method=..., acknowledge_change=True)`

`sdf` and `size_fn` are now **keyword-only**. Switch any call that passed them positionally to keyword form.

```python
# Before (≤1.3.x) — positional sdf breaks on 1.4.x
mesh.smooth_mesh("laplacian", True, my_sdf)

# After (≥1.4.0) — sdf/size_fn keyword-only
mesh.smooth_mesh("laplacian", True, sdf=my_sdf)
```

### Migrating a call site

```python
# Before (≤1.3.x)
layers = mesh.skeletonize()

# After (≥1.4.0)
layers = mesh.peel_layers()
```

To support both `chilmesh<1.4` and `>=1.4`, use this shim:

```python
peel = getattr(mesh, "peel_layers", None) or mesh.skeletonize
layers = peel()
```

### Pin guidance

- Pin `chilmesh<1.4`: safe for now. Plan the rename before you widen the pin.
- Admit 1.4.x and call any renamed symbol: broken on upgrade. Rename the call sites first, or add the `getattr` shim above.

---

## Do I Need to Change My Code?

**No**, apart from the v1.4.0 renames above. Existing code works without modification.

You should update eventually. The benefits are better performance (1.5x improvement), clearer code (bridge adapters), better error messages, and new quality metrics.

Migration timeline:

- **Now**: Update and test. The update should be a drop-in replacement.
- **Later**: Optionally adopt bridge adapters for clearer code.
- **No rush**: Confirm that your tests pass first.

---

## What's New in v1.0.0?

### Preferred Idiom (recommended for new code)
```python
from chilmesh import Mesh

mesh = Mesh.read_from_fort14("ocean.14")
```

`Mesh` is an alias for `CHILmesh` with no behaviour change. Existing imports work through the v1.x series. A `DeprecationWarning` lands in v1.1.0. Removal comes no earlier than v2.0.0.

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

Set the `CHILMESH_BACKEND` environment variable to force a backend (`python`, `cpp`, or `rust`). When it is unset, CHILmesh selects the fastest available backend (order: C++ → Rust → Python).

> **Downstream note:** the **Rust** backend is **frozen**. It stays output-equivalent but gets no further development. Build **C++** for acceleration.
> See [`RUST_EVALUATION.md`](RUST_EVALUATION.md).

### New Public Surface (cumulative since v0.4.1)
- the legacy research adapter, `MeshAdapterForADMESH`, `MeshAdapterForADMESHDomains` re-exported at package root.
- `CHILmesh.submesh(elem_ids)`: public sub-mesh factory (#138).
- `CHILmesh.rebuild_adjacencies()` / `invalidate_adjacencies()`: for mid-sweep mutation flows (#143).
- `CHILmesh.ccw_edges_around_vert(vert_id)`: counter-clockwise edge walk (#133).
- `compute_adjacencies` kwarg on the constructor and `read_from_fort14` / `read_from_2dm` / `from_admesh_domain` (#134).
- `chilmesh.backend_info()`: runtime backend introspection.

### Guarantees
- The public API is frozen under semver. Breaking changes require a v2.x bump.
- The `CHILmesh` import stays through v1.x. A DeprecationWarning arrives in v1.1. Earliest removal is v2.0.
- The C++ and Rust backends produce bit-identical output to Python on the official fixture set.
- Python `>=3.10` is required (3.8 and 3.9 are dropped; EOL).

---

## Quick Start by Project

Each project below follows the same three steps:

1. Run your old (v0.1.1) code. It still works.
2. Optionally switch to a bridge adapter.
3. Run your tests.

### Legacy research adapter Quick Start

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
from chilmesh.bridge import MeshAdapterForADMESH

mesh = CHILmesh.read_from_fort14("domain.fort.14")
adapter = MeshAdapterForADMESH(mesh)

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

### Valence Quick Start

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
Run your existing test suite. Everything should pass without code changes.

If tests fail:

1. Check the error messages.
2. Report the issue at: https://github.com/domattioli/CHILmesh/issues
3. Include the CHILmesh version, a code snippet, and the error traceback.

### Phase 3: Adopt Bridge Adapters (Optional, 2-4 hours)
Replace verbose patterns with adapter methods, one at a time:
1. Find mesh interactions in your code.
2. Check whether the adapter has a convenience method for your use case.
3. Refactor to use the adapter.
4. Run tests to verify that behavior is unchanged.

### Phase 4: Performance Profiling (Optional, 1-2 hours)
```bash
python -m cProfile -s cumtime your_script.py > profile.txt
```
Compare the profile before and after the upgrade. Bridge adapters often run faster than manual patterns.

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
from chilmesh.bridge import MeshAdapterForADMESH
adapter = MeshAdapterForADMESH(mesh)
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

## API Reference by Project

### Legacy research adapter: Mesh Adaptation Research

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

### Valence: Domain Decomposition

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

## Troubleshooting

### Error: "AttributeError: CHILmesh has no attribute 'get_vertex_edges'"

**Cause:** You use CHILmesh < 1.0.0.

**Fix:**
```bash
pip install --upgrade 'chilmesh>=1.0.0'
```

### Error: "ModuleNotFoundError: No module named 'chilmesh.bridge'"

**Cause:** You use CHILmesh < 1.0.0, or Python cannot find the bridge module.

**Fix:**
```bash
pip install --upgrade 'chilmesh>=1.0.0'
python -c "from chilmesh.bridge import MeshAdapterForADMESH; print('OK')"
```

### Error: "KeyError: 'Vert2Edge'"

**Cause:** Your code accesses internal adjacency structures directly.

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

**Problem:** Code runs slower after the upgrade. This is rare; it should run faster.

1. Verify that you run 0.2.0 or later:
   ```python
   import chilmesh
   print(chilmesh.__version__)  # Should be 0.2.0 or later
   ```

2. Profile:
   ```bash
   python -m cProfile -s cumtime your_script.py > profile.txt
   ```

3. Report the CHILmesh version, a minimal code sample, timing before and after, and the mesh size.

---

## Version Compatibility

| CHILmesh | an earlier internal prototype | ADMESH | Valence | Notes |
|----------|----------|--------|----------------|-------|
| 0.1.1    | ✅ Legacy | ✅ Legacy | ✅ Legacy | Old version, still works |
| 0.2.0    | ✅ Recommended | ✅ Recommended | ✅ Recommended | Current, use this |
| 0.2.x    | ✅ Recommended | ✅ Recommended | ✅ Recommended | Bug fixes, recommended |
| 1.0.0    | ✅ | ✅ | ✅ | `Mesh` alias, C++ backend |
| 1.4.0+   | ⚠️ | ⚠️ | ⚠️ | **#187 breaking renames; see [above](#-v140-breaking-renames-187)** |

---

## FAQ

**Q: Will my code break when upgrading?**
Not from 0.2.0 alone; CHILmesh 0.2.0 is backward compatible and all existing APIs work unchanged. Version 1.4.0 is the exception: see [v1.4.0 Breaking Renames](#-v140-breaking-renames-187).

**Q: Should I use the bridge adapters?**
Optional but recommended. Code becomes clearer and often runs faster. Adopt them at your own pace.

**Q: What if I need old behavior?**
Stay on v0.1.1. Upgrading is recommended for the performance benefits.

**Q: Can I mix old and new code?**
Yes. Old patterns and bridge adapters work in the same codebase during migration.

**Q: Is the bridge adapter API stable?**
Yes. Adapter methods are part of CAI and stay stable through v1.0.

**Q: How much code do I need to change?**
None is required. Migration is optional.

---

## Getting Help

### Questions About Migration
1. Check this guide.
2. Review `docs/CHILmesh_Access_Interface.md`.
3. Look at `tests/test_bridge_adapters.py`.

### Bug Reports
Include:
- CHILmesh version: `python -c "import chilmesh; print(chilmesh.__version__)"`
- Python version: `python --version`
- A minimal reproducible example
- The full error traceback
- The downstream project name (ADMESH/Valence)

Report at: https://github.com/domattioli/CHILmesh/issues

---

## Example Integration Projects

The `examples/` directory holds complete working examples:
- `refinement.py` - legacy research mesh adaptation workflow
- `admesh_quality.py` - ADMESH quality assessment
- `registry_setup.py` - Valence domain initialization

---

## Support Timeline

**CHILmesh 0.2.0 Support:**
- Release date: April 2026
- Bug fix support: Until 0.3.0 release
- Security support: Until 1.0.0 release

**Migration Support:**
- Expect 1-4 hours to update your project.
- Most of that time goes to testing, not code changes.

---

**Last Updated:** 2026-07-21
**Guide Version:** 2.1

For latest information, visit: https://github.com/domattioli/CHILmesh
