# PHASE 1: Unified Feature Specification

**"Research & Design CHILmesh's Modern Graph Representation & Mesh Operations Framework"**

**Branch:** `planning-optimize_modernize`  
**Date:** 2026-04-26  
**Scope:** Consolidates Issues #35–38 into ONE coherent feature  
**Methodology:** Spec-Kit (Specify Phase)

---

## EXECUTIVE SUMMARY

Phase 1 is ONE unified research feature:

> **Design modern graph representation for CHILmesh enabling dynamic mesh operations, preserving skeletonization invariants, supporting downstream integration (MADMESHR advancing-front, ADMESH-Domains bulk I/O).**

Four research areas (#35–38) are interdependent threads:
1. **Graph Structure Choice** (#35) determines feasibility of **Dynamic Operations** (#37)
2. **Dynamic Operations** (#37) inform **Pinch-Point Detection** (#38) API design
3. **Skeletonization Analysis** (#36) validates chosen structure preserves **Skeletonization Invariants**
4. All three feed into **Decision Memo** (Phase 1 output)

---

## PART 1: CURRENT STATE

### 1.1 What Exists Today (CHILmesh 0.1.x)

```python
class CHILmesh:
  connectivity_list: np.ndarray (M×3 or M×4)    # Element → Vertex
  points: np.ndarray (N×3)                       # Vertex coordinates
  adjacencies: Dict[str, Any]                    # 6 structures:
    - Elem2Vert, Edge2Vert, Elem2Edge (List[List[int]])  # ← Bottleneck: O(n²)
    - Vert2Edge, Vert2Elem (List[List[int]])
    - Edge2Elem: np.ndarray (E×2)   # sentinel -1 for boundary
  layers: Dict                                   # OE, IE, OV, IV, bEdgeIDs
```

**Known Performance Hotspots:**
- `_identify_edges()`: O(n²) — enumerate all element vertex pairs
- `_build_elem2edge()`: O(n²) — match elements to edges via linear search
- Adjacency maintenance: dual-representation burden (Elem2Edge + Vert2Elem must sync)

**Correctness Issues (Audit B1–B11):**
- B3: Padded-triangle sentinel (vertex[3] == 0) silently demotes quads
- B4: Quad flip permutation breaks on padded triangles

### 1.2 Downstream Consumer Needs

**MADMESHR** (RL-based advancing-front mesh generator):
- Element insertion O(1) or O(log n) — 1000s episodes × 100+ placements each
- Boundary edge enumeration O(n), O(1) per-edge lookup
- Domain splitting (pinch-point detection); residual closure (≤4 vertex boundary auto-fill)

**ADMESH-Domains** (mesh registry/bulk I/O):
- Block_O (~15K nodes) init currently 30s — unacceptable for interactive bulk load
- Need: <1s initialization per mesh (50× speedup target)
- Metadata queries; batch export (already working)

**ADMESH** (GitHub 404, future): edge-swapping, node repositioning, refinement/coarsening

### 1.3 Hard Constraints

1. **Skeletonization Semantics Locked** (Decision A2) — Layer structure immutable; tests pin as regression suite
2. **API Backward Compatibility** — Public method signatures unchanged; internal adjacencies dict may refactor
3. **.fort.14 Roundtrip** (Audit B1) — Byte-identical I/O; graph structure must map cleanly
4. **Mixed-Element Support** (Audit B4) — vertex3 == vertex0 for triangles in 4-column connectivity
5. **No New External Dependencies** — scipy, numpy, matplotlib only

---

## PART 2: INTERCONNECTED RESEARCH THREADS

Cross-cutting concerns across all 4 issues:

| Concern | Affects | Implication |
|---------|---------|-------------|
| O(n²) bottleneck | #35, #36 | Graph structure MUST support O(n log n) or O(n) edge discovery |
| Skeletonization invariants | #35, #36, #37 | Any graph structure MUST preserve exact layer semantics |
| Mixed-element handling | #35, #36, #37 | Must clarify element-type handling (enum vs. padding rule) |
| Dynamic ops feasibility | #35, #37, #38 | Graph structure must support O(1) element addition (no full rebuild) |
| API stability | All | Phase 1 design frozen; Phase 2 implements exact signatures |

---

## PART 3: UNIFIED RESEARCH DELIVERABLES

```
research/
  ├── GRAPH_COMPARISON.md  — All 4 structures; benchmarks on 4 fixtures; recommendation: Compact Graph
  ├── SKELETONIZATION_ANALYSIS.md  — Algorithm pseudo-code; invariants; optimization opportunities
  ├── DYNAMIC_OPS_DESIGN.md  — 6 core operations; pre/post conditions; transactional model
  ├── PINCH_POINT_DETECTION.md  — 3+ definitions; algorithm survey; reuse skeletonization; API sketch
  └── DECISION_MEMO.md  — Compact Graph recommendation; validates constraints; top 3 risks; signed approval gate
```

### Success Criteria

#### Benchmark Framework (Issue #35):
- [ ] All 4 structures implemented + tested on 4 fixtures
- [ ] Block_O benchmark shows >10× improvement potential

#### Skeletonization Validation (Issue #36):
- [ ] Algorithm documented with pseudo-code
- [ ] Invariants verified on 4 test fixtures
- [ ] Analysis confirms chosen graph structure preserves output exactly

#### Dynamic Operations Design (Issue #37):
- [ ] 6 operations enumerated with complexity analysis and pre/post conditions
- [ ] Feasibility demonstrated with chosen graph structure

#### Pinch-Point Detection (Issue #38):
- [ ] 2+ definitions of "pinch point" with examples
- [ ] Proposes reusing skeletonization layers; API signature sketched

#### Decision Memo (Feature output):
- [ ] Recommends Compact Graph with benchmark evidence
- [ ] All constraints validated (skeletonization, .fort.14, mixed-element, API)
- [ ] Risk register (top 3 risks + mitigation)
- [ ] **Approval gate:** Signed off before Phase 2

#### Testing & Documentation:
- [ ] All existing tests pass (0 regressions)
- [ ] GitHub Issues #35–38 closed; Epic #39 remains open

---

## PART 4: TIMELINE

```
Week 1: Day 1–2: Graph benchmarking (#35)
        Day 2–3: Skeletonization analysis (#36)
        Day 3–4: Dynamic ops design (#37)
        Day 4–5: Pinch-point detection (#38)
Week 2: Day 6–7: Cross-linking & validation
        Day 8–9: Decision memo → sign-off gate
```

**Approval Gate:** Decision memo reviewed + approved before Phase 2 begins.

---

**Document Status:** UNIFIED FEATURE SPECIFICATION  
**Author:** Spec-Kit Specify Phase (Claude)  
**Date:** 2026-04-26  
**Branch:** `planning-optimize_modernize`  
**Related Issues:** #35, #36, #37, #38 (consolidated into ONE feature)  
**Parent Epic:** #39 (Data Structure Modernization)
