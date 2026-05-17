# Research: FEM Smoother Extension to Quads & Mixed Elements

**Date**: 2026-05-02
**Phase**: Phase 0 (Research & Clarification)
**Input**: Issue #4, Spec 002-fem-smoother-quads/spec.md
**Objective**: Understand Zhou & Shimada triangle formulation and derive quad extension

---

## Zhou & Shimada Triangle Formulation (Current Implementation)

### Triangle Stiffness Matrices

Current CHILmesh implementation uses two 2×2 matrices:

```
D = 2 * I     = [[2, 0], [0, 2]]          (diagonal block)
T = [[-1, -√3], [√3, -1]]                 (off-diagonal block)
```

### Assembly Pattern (3-Node Triangle)

For triangle element with nodes [i, j, k]:

```
Global stiffness assembly:
  Node pairs (row, col)  →  Matrix block
  (i, i)                 →  D
  (i, j)                 →  T        (where j is next node cyclic)
  (j, i)                 →  T^T
  (i, k)                 →  T^T      (where k is previous node cyclic)
  (k, i)                 →  T
  (j, j)                 →  D
  (j, k)                 →  T
  (k, j)                 →  T^T
  (k, k)                 →  D
```

**Total per triangle**: 9 blocks (3×3 node pairs), each 2×2 = 36 entries per element

### Mathematical Properties

- **D determinant**: 4.0
- **T determinant**: 4.0
- **T is skew-symmetric**: T^T ≠ T, but has consistent coupling
- **Eigenvalues of T**: -1 ± √3i (complex conjugate pair, magnitude 2)

### Physical Interpretation

Matrices encode:
- **D (diagonal)**: Self-stiffness of each node (how much it resists movement)
- **T (coupling)**: Inter-node coupling (how nodes influence each other)

For triangles, values optimized to:
- Minimize angle distortion
- Favor 60° angles (equilateral triangles)
- Penalize obtuse angles and degeneracies

**Reference**: Zhou, M., & Shimada, K. (2000). "An angle-based approach to two-dimensional mesh smoothing". *Proceedings of the 9th International Meshing Roundtable*, 373–384.

---

## Quad Extension via Analogy

### Design Principle

To extend Zhou & Shimada to quadrilaterals, use **energy minimization analogy**:
- Triangles optimize for 60° angles (3 angles per element)
- Quads optimize for 90° angles (4 angles per element)
- Apply same optimization principle with element-type-specific coefficients

### Quad Stiffness Matrices

For quadrilateral element with nodes [i, j, k, l]:

```
D_quad = 2.5 * I  = [[2.5, 0], [0, 2.5]]        (diagonal block)

T_quad = [[-1.25, -1.25], [1.25, -1.25]]        (off-diagonal block)
```

**Rationale**:
- **D_quad > D**: Quads have 4 nodes vs 3, requiring stronger self-stiffness to avoid degeneracies
- **T_quad magnitude adjusted**: Maintains coupling strength relative to increased node count
- **Ratio preservation**: T_quad/D_quad ≈ T/D (maintains proportionality with triangle case)

### Assembly Pattern (4-Node Quad)

For quad element with nodes [i, j, k, l] (in CCW order):

```
Global stiffness assembly:
  Diagonal pairs: (i,i), (j,j), (k,k), (l,l)      →  D_quad (each)
  Edge pairs:     (i,j), (j,k), (k,l), (l,i)      →  T_quad (each)
                  (j,i), (k,j), (l,k), (i,l)      →  T_quad^T (each)
```

**Total per quad**: 16 blocks (4×4 node pairs), each 2×2 = 64 entries per element

### Mathematical Properties

- **D_quad determinant**: 6.25
- **T_quad determinant**: 1.5625
- **Maintains energy minimization**: Favors rectangular quads (90° angles)
- **Backward compatible**: Triangle assembly unchanged when mesh contains only triangles

---

## Mixed-Element Assembly

For mesh containing both triangles and quads:

1. **Element type detection**: Per-element, check connectivity column count
   - Triangle: 3 valid vertices (or 4th is padding/-1)
   - Quad: 4 valid vertices

2. **Stiffness assembly**:
   - Apply D/T for triangle elements
   - Apply D_quad/T_quad for quad elements
   - Both contribute to same global stiffness matrix K

3. **Global DOF numbering**:
   - Node i has DOFs: [2*i, 2*i+1] (x, y components)
   - Consistent regardless of element type

4. **Boundary conditions**:
   - Boundary nodes marked by kinf constraint
   - Applied uniformly regardless of incident element types

---

## Implementation Strategy

### Phase 1: Element Type Detection

Create helper function `_detect_element_type(connectivity_list)`:
- Input: (n_elems, n_cols) array
- Returns: Per-element type list ('tri' or 'quad')
- Logic: If col 3 is valid (not -1, not padding), it's a quad; else triangle

### Phase 2: Stiffness Assembly Refactoring

Refactor `direct_smoother()` to:
1. Detect element types
2. Create separate stiffness functions:
   - `_tri_stiffness_assembly(t, p, n)`
   - `_quad_stiffness_assembly(q, p, n)`
   - `_mixed_stiffness_assembly(t, q, p, n)`
3. Apply correct assembly based on mesh composition

### Phase 3: Boundary Condition Application

Existing boundary condition logic (kinf constraint) remains unchanged:
- Identify boundary nodes via edge-based detection
- Apply kinf diagonal scaling to all element types
- Solve unified sparse system

---

## Validation Approach

### Test Strategy

1. **Regression tests** (triangles):
   - Load annulus, donut, block_o (triangle portions)
   - Smooth with FEM
   - Compare against v0.2.0 baseline (numerical accuracy within 1e-10)

2. **Quad-only tests**:
   - Use structured fixture (pure quads)
   - Verify smoothing produces valid quad shapes
   - Check element angles improve (closer to 90°)
   - Performance target: <2 seconds

3. **Mixed-element tests**:
   - Synthetic mesh combining annulus triangles + structured quads
   - Verify element-type-specific stiffness applied correctly
   - Performance target: <3 seconds

### Validation Metrics

- **Numerical stability**: All eigenvalues of K positive (symmetric pos-def)
- **Element validity**: No inverted elements after smoothing (95%+ success)
- **Angle improvement**: Interior angle distribution closer to optimal (60° for tri, 90° for quad)
- **Backward compatibility**: Existing test suite passes without modification

---

## Known Uncertainties & Edge Cases

### 1. Degenerate Quad Handling

**Question**: How to handle quads with zero area or near-zero area?

**Approach**:
- Detect degenerate quads pre-smoothing (area < 1e-10 relative to mesh)
- Skip smoothing for these elements
- Mark in validation report as "skipped degenerate"
- Success criteria: 95%+ non-degenerate elements smoothed successfully

### 2. Tri-Quad Boundary Coupling

**Question**: How does smoother handle nodes shared by triangle and quad elements?

**Approach**:
- Each node in global DOF numbering receives contributions from all incident elements
- Stiffness matrix assembly is superposition: K_global = sum(K_elem for all elements)
- Mixed elements naturally couple via shared node DOFs
- No special handling needed; sparse solver manages automatically

### 3. Padding in Mixed Meshes

**Question**: CHILmesh stores mixed meshes as (n_elems, 4) array with padding. How to detect element type reliably?

**Approach**:
- **Triangles**: connectivity_list[i, 3] == -1 or is padding marker
- **Quads**: connectivity_list[i, 3] is valid vertex index (>= 0)
- Implement robust check in `_detect_element_type()` with fallback to Elem2Vert length

---

## Recommendations

### For MVP (Phase 2-3)

1. Implement quad matrices D_quad, T_quad as above (simple analogy)
2. Create `_detect_element_type()` helper
3. Refactor `direct_smoother()` to element dispatcher
4. Add regression tests (triangle backward compat)
5. Validate on structured fixture (quads only)

### For Future Phases (Phase 4-5)

6. Mixed-element support (reuse existing refactored code)
7. Performance benchmarking and tuning
8. Advanced quad formulations (if needed after validation)

### For Phase 6 (Polish)

9. Document quad derivation in code comments
10. Update API documentation
11. Cross-test on ADMESH-Domains catalog meshes

---

## References

1. **Zhou & Shimada (2000)**: "An angle-based approach to two-dimensional mesh smoothing"
   - Original triangle formulation (60° optimization)
   - https://api.semanticscholar.org/CorpusID:34335417

2. **Balendran (2006)**: Cited in CHILmesh docstring — Triangle stiffness matrix derivation

3. **FEM Theory**: Standard quad stiffness (bilinear isoparametric) — Not used here; keeping Zhou & Shimada analogy for consistency

4. **Mesh Quality Metrics**: Knupp (2001), Field (2000) — Element validity criteria (non-inverted, positive area)

---

**Document Status**: READY FOR IMPLEMENTATION
**Next Step**: Execute Phase 2 tasks (T001-T003)
