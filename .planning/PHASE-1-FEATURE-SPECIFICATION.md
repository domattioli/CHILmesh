# PHASE 1: Unified Feature Specification

**"Research & Design CHILmesh's Modern Graph Representation & Mesh Operations Framework"**

**Branch:** `planning-optimize_modernize`  
**Date:** 2026-04-26  
**Scope:** Consolidates Issues #35–38 into ONE coherent feature  
**Methodology:** Spec-Kit (Specify Phase) — Address Cross-Cutting Concerns

---

## EXECUTIVE SUMMARY

Phase 1 is **NOT four separate research tasks**. It's **ONE unified research feature**:

> **Design a modern graph representation for CHILmesh that enables dynamic mesh operations, preserves skeletonization invariants, and supports downstream integration (MADMESHR advancing-front, ADMESH-Domains bulk I/O).**

The four research areas (#35–38) are **interdependent threads** of this single feature:

1. **Graph Structure Choice** (#35) determines feasibility of **Dynamic Operations** (#37)
2. **Dynamic Operations** (#37) inform **Pinch-Point Detection** (#38) API design
3. **Skeletonization Analysis** (#36) validates that chosen structure preserves **Skeletonization Invariants**
4. **All three** feed into **Decision Memo** (output of Phase 1)

---

## PART 1: CURRENT STATE (Spec-Kit Specify Phase)

### 1.1 What Exists Today (CHILmesh 0.1.x)

**Core Graph Representation:**
```python
class CHILmesh:
  connectivity_list: np.ndarray (M×3 or M×4)    # Element → Vertex
  points: np.ndarray (N×3)                       # Vertex coordinates
  adjacencies: Dict[str, Any]                    # 6 structures:
    - Elem2Vert: connectivity reference
    - Edge2Vert: np.ndarray (E×2)
    - Elem2Edge: List[List[int]]    # ← Bottleneck: O(n²) discovery
    - Vert2Edge: List[List[int]]
    - Vert2Elem: List[List[int]]
    - Edge2Elem: np.ndarray (E×2)   # with sentinel -1 for boundary
  layers: Dict                                   # Skeletonization result
    - OE, IE, OV, IV, bEdgeIDs: List[np.ndarray]
```

**Known Performance Hotspots:**
- `_identify_edges()`: O(n²) — enumerate all element vertex pairs
- `_build_elem2edge()`: O(n²) — match elements to edges via linear search
- Skeletonization: Redundant set membership checks per iteration
- Adjacency maintenance: Dual-representation burden (Elem2Edge + Vert2Elem must sync)

**Known Correctness Issues (Audit B1–B11):**
- B3: Padded-triangle sentinel (vertex[3] == 0) silently demotes quads
- B4: Quad flip permutation breaks on padded triangles
- Recursion depth bug in smoothing (tests caught, code review missed)

---

### 1.2 Downstream Consumers & Their Needs

#### **MADMESHR (Advancing-Front RL Mesh Generator)**

**What MADMESHR does:**
- Uses RL (Soft Actor-Critic) to place quads/triangles on advancing-front boundary
- Iteratively removes boundary loops, inserts elements
- Detects pinch points (bottlenecks) and splits domain
- Performs residual closure (final <4 vertex region auto-fill)

**Critical Operations MADMESHR Will Perform on CHILmesh:**
1. **Element Insertion** (advancing-front placement)
   - Must happen in real-time during RL episode (1000s of episodes × 100+ placements each)
   - Cannot afford O(n²) rebuild per insertion
   - Need: O(1) or O(log n) element addition

2. **Boundary Edge Enumeration** (frontier detection)
   - Find all edges on current mesh boundary
   - Need: O(n) enumeration, O(1) per-edge lookup

3. **Domain Splitting** (pinch-point detection)
   - Identify narrow bottlenecks in domain
   - Split at narrowest edge
   - Need: Efficient bottleneck identification (can it reuse skeletonization?)

4. **Residual Closure** (final auto-fill)
   - When boundary shrinks to ≤4 vertices, auto-triangulate remainder
   - Need: Detect when boundary is small

**Implication for Phase 1 Research:**
- Graph structure must support O(1) or O(log n) element insertion
- Skeletonization analysis must determine if layer data can identify pinch points
- Dynamic ops design must specify element insertion pre/post conditions
- Decision memo must justify that chosen structure supports all 4 operations

#### **ADMESH-Domains (Mesh Registry & Bulk I/O)**

**What ADMESH-Domains does:**
- Hosts ~150+ coastal/riverine domain meshes (pypi + HuggingFace)
- Bulk loads meshes (up to 1M nodes) on demand
- Filters by metadata (size, quality, domain type)
- Exports subsets to `.fort.14`

**Critical Operations ADMESH-Domains Will Perform:**
1. **Fast Mesh Loading** (initialization bottleneck)
   - Currently: Block_O (~15K nodes) takes ~30s due to O(n²) edge discovery
   - Unacceptable for interactive use (bulk load 10 meshes = 5 minutes)
   - Need: <1s initialization per mesh (50× speedup target)

2. **Metadata Queries** (search/filter)
   - Size, quality, element type, boundary complexity
   - Need: Fast metadata extraction (skeletonization optional)

3. **Batch Export** (subset to `.fort.14`)
   - Export subset of mesh to ADCIRC format
   - Need: Efficient I/O (already working; no change)

**Implication for Phase 1 Research:**
- Graph structure O(n²) edge discovery is **urgent blocker** (motivates Phase 1 timeline)
- Phase 3 optimization (O(n log n)) directly solves this pain point
- Decision memo must show >10× improvement potential to justify effort
- Skeletonization optional for ADMESH-Domains (can defer layer computation)

#### **ADMESH (Hypothetical Mesh Adaptation, GitHub 404)**

**Status:** Not yet public. Assume future availability.

**Likely Operations (from thesis context):**
- Edge-swapping (improve mesh quality, reduce high-valence vertices)
- Node repositioning (smoothing + quality improvement)
- Refinement/coarsening (local mesh modification)

**Implication for Phase 1 Research:**
- Dynamic ops API (#37) must support `swap_edge()`, `remove_vertex()`, `add_element()`
- These operations are Phase 2 deliverables, but design must anticipate ADMESH needs
- No direct ADMESH integration in Phase 4 (GitHub 404), but APIs future-proof

---

### 1.3 Hard Constraints (MUST Preserve)

1. **Skeletonization Semantics Locked** (Audit Q3 Decision A2)
   - Layer structure (OE, IE, OV, IV, bEdgeIDs) immutable
   - Disjoint cover + monotone-shrinking invariants must hold exactly
   - Tests pin this as regression suite
   - *Implication for Phase 1:* Graph structure choice cannot change skeletonization output

2. **API Backward Compatibility (0.1.x → 0.2.0)**
   - Public methods (`signed_area()`, `elem_quality()`, `plot()`) keep signatures
   - Internal: adjacencies dict may refactor (0.2.0 breaking change acceptable)
   - *Implication for Phase 1:* Design must distinguish public vs. internal APIs

3. **.fort.14 Roundtrip** (Audit B1)
   - Byte-identical I/O roundtrip (mesh → .fort.14 → mesh)
   - Cannot change connectivity storage format
   - *Implication for Phase 1:* Graph structure must map cleanly to/from `.fort.14`

4. **Mixed-Element Support** (Audit B4)
   - Triangles + quads in single mesh
   - Padded-triangle semantic (vertex3 == vertex0 for triangles)
   - *Implication for Phase 1:* Graph structure must handle 3-vertex and 4-vertex elements explicitly

5. **No New External Dependencies**
   - scipy, numpy, matplotlib only (no NetworkX as runtime dependency)
   - *Implication for Phase 1:* Rule out options requiring heavy libraries

---

### 1.4 Soft Constraints (Nice to Have)

- Backwards compatibility with `_mesh_layers()` (deprecated wrapper)
- Plotting support (matplotlib, non-critical path)
- Type hints throughout (Python 3.10+)
- Single-threaded execution (no parallelization)

---

## PART 2: INTERCONNECTED RESEARCH THREADS

### 2.1 How Issues #35–38 Are Interdependent

```
┌─────────────────────────────────────────────────────────┐
│ PHASE 1 FEATURE: Modern Graph Representation            │
│                                                         │
│  #35: Graph Structure      ─┐                          │
│        (NetworkX vs        │                           │
│         Compact vs    ◄────┴─── Informs #37 API       │
│         CSR vs        ├─── Validates #36 Invariants   │
│         Half-Edge)    └─── Enables #38 Algorithm      │
│                                                         │
│  #36: Skeletonization     ─┐                          │
│        Analysis     ◄──────┴─── Validates #35 Choice  │
│                     └─────────── Informs #38 Reuse    │
│                                                         │
│  #37: Dynamic Ops  ◄──┬────────── Requires #35 Choice │
│        Design      │   └────────── Uses #36 Invariant │
│                    └─────────────┐                    │
│                                  ▼                    │
│  #38: Pinch-Point ◄──────────── Reuses #36 Layers    │
│        Detection     ◄────────── Requires #37 API     │
│                                                         │
│  DECISION MEMO ◄──────────────────────────────────────│
│        (Output of Phase 1)                             │
│  • Recommend graph structure (#35 choice)              │
│  • Validate skeletonization preserved (#36 proof)      │
│  • Confirm dynamic ops feasible (#37 design)           │
│  • Sketch pinch-point API (#38 algorithm)              │
└─────────────────────────────────────────────────────────┘
```

### 2.2 Cross-Cutting Concerns

**Concern 1: Performance Bottleneck (O(n²) Edge Discovery)**
- **Affects:** #35 (graph structure choice), #36 (skeletonization analysis)
- **Root Cause:** Current `_identify_edges()` enumerates all element pairs
- **Implication:** Graph structure MUST support O(n log n) or O(n) edge discovery
- **Metric:** Block_O init <1s (vs ~30s) validates improvement

**Concern 2: Skeletonization Invariants**
- **Affects:** #35 (graph structure), #36 (algorithm analysis), #37 (dynamic ops)
- **Root Cause:** Audit Q3 confirmed disjoint cover + monotone-shrinking as valid invariants
- **Implication:** Any graph structure MUST preserve exact layer semantics
- **Test:** Skeletonization output byte-identical after dynamic ops

**Concern 3: Mixed-Element Handling**
- **Affects:** #35 (graph structure), #36 (skeletonization), #37 (dynamic ops)
- **Root Cause:** Padded-triangle semantic (vertex3 == vertex0) is fragile in Python
- **Implication:** Graph structure must clarify element-type handling (enum vs. padding rule)
- **Test:** Synthetic test case (triangle in 4-column with vertex 0 in slot 3)

**Concern 4: Dynamic Operations Feasibility**
- **Affects:** #35 (graph structure choice), #37 (API design), #38 (pinch-point API)
- **Root Cause:** MADMESHR needs real-time element insertion; static-only API insufficient
- **Implication:** Graph structure must support O(1) element addition (no full rebuild)
- **Test:** Benchmarks show element insertion feasible in Phase 3

**Concern 5: API Stability**
- **Affects:** All issues (#35–38), MADMESHR integration (Phase 4)
- **Root Cause:** Once we commit to dynamic ops API, MADMESHR builds on it
- **Implication:** Phase 1 design must be frozen; Phase 2 implements exact signatures
- **Test:** MADMESHR team reviews Phase 1 decision memo

---

## PART 3: UNIFIED RESEARCH DELIVERABLES

Instead of 4 separate research artifacts (Issues #35–38), deliver **ONE coherent feature specification**:

### 3.1 Unified Deliverable Structure

```
research/
  ├── GRAPH_COMPARISON.md
  │   ├── All 4 structures (NetworkX, Compact, CSR, Half-Edge)
  │   ├── Benchmarks on 4 fixtures (annulus, block_o, structured, synthetic 100K)
  │   ├── Comparison table (construction time, lookup, traversal, memory)
  │   └── Recommendation: Compact Graph with rationale
  │
  ├── SKELETONIZATION_ANALYSIS.md
  │   ├── Current algorithm (pseudo-code)
  │   ├── Complexity analysis
  │   ├── Invariant verification (disjoint cover, monotone-shrinking)
  │   ├── Identifies optimization opportunities
  │   └── Validates that chosen graph structure preserves output exactly
  │
  ├── DYNAMIC_OPS_DESIGN.md
  │   ├── 6 core operations (add/remove element, vertex, split/swap edge)
  │   ├── For each: pre/post conditions, consistency rules, complexity
  │   ├── Transactional model (batch operations)
  │   ├── Demonstrates feasibility of chosen graph structure
  │   └── Sketches integration with skeletonization
  │
  ├── PINCH_POINT_DETECTION.md
  │   ├── Defines "pinch point" (3+ interpretations with examples)
  │   ├── Surveys detection algorithms (medial axis, bottleneck sampling, distance field)
  │   ├── Proposes reusing skeletonization (can layer data identify pinches?)
  │   ├── Sketches API signature for Phase 2/3
  │   └── Justifies feasibility with chosen graph + dynamic ops APIs
  │
  └── DECISION_MEMO.md
      ├── Recommends Compact Graph (Option B) based on benchmarks
      ├── Validates skeletonization invariants preserved
      ├── Confirms dynamic ops feasible with chosen structure
      ├── Sketches pinch-point detection as reuse of skeletonization
      ├── Top 3 risks + mitigation for Phase 2
      └── Signed approval gate (required before Phase 2)
```

### 3.2 Unified Success Criteria

**Feature "Research & Design Modern Graph Representation" is COMPLETE when:**

#### **Benchmark Framework (Issue #35 core):**
- [ ] `research/GRAPH_COMPARISON.md` written
- [ ] All 4 structures implemented + tested on 4 fixtures
- [ ] Comparison table shows Compact Graph beats others on ≥2 metrics
- [ ] Block_O benchmark shows >10× improvement potential (validated by mock Compact Graph)

#### **Skeletonization Validation (Issue #36 core):**
- [ ] `research/SKELETONIZATION_ANALYSIS.md` written
- [ ] Current algorithm documented with pseudo-code
- [ ] Invariants verified on 4 test fixtures
- [ ] Optimization opportunities identified (3+ with tradeoffs)
- [ ] **Cross-link:** Analysis confirms chosen graph structure preserves output exactly

#### **Dynamic Operations Design (Issue #37 core):**
- [ ] `research/DYNAMIC_OPS_DESIGN.md` written
- [ ] 6 operations enumerated with complexity analysis
- [ ] Pre/post conditions clarified for each
- [ ] Transactional model sketched
- [ ] **Cross-link:** Design demonstrates feasibility with chosen graph structure

#### **Pinch-Point Detection (Issue #38 core):**
- [ ] `research/PINCH_POINT_DETECTION.md` written
- [ ] 2+ definitions of "pinch point" with examples
- [ ] 3+ algorithms compared (complexity, accuracy, reusability)
- [ ] **Cross-link:** Proposes reusing skeletonization layers (from analysis in #36)
- [ ] API signature sketched for Phase 2/3 implementation

#### **Decision Memo (Feature output):**
- [ ] `research/DECISION_MEMO.md` written
- [ ] Recommends Compact Graph with benchmark evidence
- [ ] Validates all constraints preserved (skeletonization, .fort.14, mixed-element, API)
- [ ] Confirms Phase 2 dynamic ops feasible
- [ ] Confirms Phase 3 pinch-point detection feasible
- [ ] Risk register (top 3 risks + mitigation)
- [ ] **Approval gate:** Decision memo reviewed + signed off by team/lead before Phase 2

#### **Testing:**
- [ ] All existing tests pass (0 regressions)
- [ ] Benchmark fixtures load without error
- [ ] No new bugs introduced

#### **Documentation:**
- [ ] All 5 research docs linked in CLAUDE.md
- [ ] GitHub Issues #35–38 closed (all addressed)
- [ ] Epic #39 remains open (tracks Phase 1 + 2 + 3 + 4)

---

## PART 4: INTEGRATED CONSTRAINTS & TRADEOFFS

### 4.1 How Each Constraint Affects All Four Research Areas

| Constraint | Affects #35 (Graph) | Affects #36 (Skeletonization) | Affects #37 (Dynamic Ops) | Affects #38 (Pinch-Points) |
|-----------|-------------------|-------------------------------|-------------------------|--------------------------|
| **Skeletonization immutable** | Must choose structure that preserves exact output | Algorithm analysis validates structure choice | Dynamic ops cannot change layer semantics | Pinch-point detection must reuse layers exactly |
| **Mixed-element support** | Structure must handle 3-vert + 4-vert elements | Algorithm works with mixed types; verify | Element insertion must respect padding rule | Pinch-point detection across mixed types |
| **O(n²) bottleneck** | Primary selection criterion for structure | Analysis reveals why O(n²) emerges; optimization path | Dynamic ops cannot reintroduce O(n²) per insertion | Pinch-point detection must be O(n log n) or O(n) |
| **.fort.14 roundtrip** | Structure must map cleanly to/from format | Skeletonization doesn't affect I/O; verify independence | Dynamic ops don't change connectivity storage | Pinch-point detection doesn't affect I/O |
| **No new dependencies** | Rules out NetworkX (has runtime overhead) | Scipy/numpy acceptable for analysis | Transactional ops use only numpy | Pinch-point detection uses only numpy |

### 4.2 Risk Mitigation via Unified Specification

**Risk: Graph structure choice is wrong** → Mitigation: Benchmark all 4 structures (#35), validate choice against all constraints (skeletonization, dynamic ops, pinch-points)

**Risk: Skeletonization invariants broken** → Mitigation: Analyze algorithm, verify on test fixtures (#36), test output after dynamic ops (#37)

**Risk: Dynamic ops API doesn't work** → Mitigation: Design before implementation (#37), validate feasibility with chosen structure (#35), sketch pinch-point usage (#38)

**Risk: Pinch-point detection too slow** → Mitigation: Propose reusing skeletonization (#38), validate with O(n log n) graph structure (#35)

---

## PART 5: CONSOLIDATED TIMELINE

**ALL research work (Issues #35–38) in ONE 2–3 week feature cycle:**

```
Week 1 (2026-04-26 → 2026-05-03):
  Day 1–2: Graph benchmarking (#35 core)
    - Implement 4 graph structures (minimal, proof-of-concept)
    - Benchmark on 4 fixtures
    - Generate comparison table
    
  Day 2–3: Skeletonization analysis (#36 core)
    - Document current algorithm
    - Verify invariants on fixtures
    - Identify optimization opportunities
    - Cross-validate with #35 choice
    
  Day 3–4: Dynamic ops design (#37 core)
    - Enumerate 6 operations
    - Specify pre/post/consistency rules
    - Sketch transactional model
    - Validate feasibility with #35 structure
    
  Day 4–5: Pinch-point detection (#38 core)
    - Define "pinch point" precisely
    - Survey algorithms
    - Propose reusing #36 skeletonization
    - Sketch API for Phase 2/3

Week 2 (2026-05-03 → 2026-05-10):
  Day 6–7: Cross-linking & validation
    - Verify all 4 research areas consistent
    - Validate constraints satisfied
    - Identify edge cases
    - Risk assessment
    
  Day 8–9: Decision memo
    - Write DECISION_MEMO.md
    - Summarize findings from all 4 areas
    - Recommend Compact Graph
    - Sign-off gate

Timeline: 2 weeks → Phase 1 research complete
```

---

## PART 6: GOVERNANCE & ACCOUNTABILITY

### 6.1 Unified Success Criteria (Feature-Level, Not Issue-Level)

**Feature:** "Research & Design Modern Graph Representation & Mesh Operations API"

**Status:** ACCEPTED when:
- [ ] All 5 research documents written + linked
- [ ] Benchmarks validate >10× improvement potential
- [ ] Skeletonization invariants verified preserved
- [ ] Dynamic ops API sketched + feasibility confirmed
- [ ] Pinch-point detection algorithm proposed
- [ ] Decision memo recommends Compact Graph with full rationale
- [ ] Risk register identifies + mitigates top 3 risks
- [ ] All existing tests pass (0 regressions)
- [ ] Team/lead signs off decision memo

**Approval Gate:** Decision memo must be reviewed + approved before Phase 2 begins.

### 6.2 Dependency on Downstream Validation

**Phase 1 → Phase 4 Feedback Loop:**
- Phase 1 decision memo proposes APIs for dynamic ops, pinch-point detection
- Phase 4 includes "Coordinate with MADMESHR team" to validate APIs match needs
- If MADMESHR feedback contradicts Phase 1 design, escalate to architecture review

---

## CONCLUSION

**Phase 1 is NOT four research tasks. It's ONE unified feature:** designing a modern graph representation for CHILmesh.

The four research areas (#35–38) are **interdependent threads**:
- Graph structure choice (#35) is validated by skeletonization analysis (#36)
- Dynamic operations design (#37) is feasible only with right graph structure (#35)
- Pinch-point detection (#38) reuses skeletonization data (#36) and dynamic ops API (#37)
- All feed into decision memo (Phase 1 output)

**Benefits of unified specification:**
- ✅ Cross-cutting concerns identified early
- ✅ Constraints validated across all areas
- ✅ Risks mitigated holistically (not in silos)
- ✅ Timeline reduced (parallelizable within feature)
- ✅ Decision memo integrates all findings

**Next Step:** Execute Phase 1 as unified feature (2-week timeline). Deliver decision memo + 5 research documents. Gate: all documents must be internally consistent + approved before Phase 2.

---

**Document Status:** UNIFIED FEATURE SPECIFICATION  
**Author:** Spec-Kit Specify Phase (Claude)  
**Date:** 2026-04-26  
**Branch:** `planning-optimize_modernize`  
**Related Issues:** #35, #36, #37, #38 (consolidated into ONE feature)  
**Parent Epic:** #39 (Data Structure Modernization)
