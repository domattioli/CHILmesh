# CHILmesh Speckit Constitution

**Version:** 1.0  
**Effective Date:** 2026-05-14  
**Scope:** Specification structure and governance for CHILmesh feature development

---

## Purpose

This document defines the standard structure, governance, and lifecycle for all feature specifications in the CHILmesh project. It ensures consistency across specs, maintains quality standards, and facilitates clear communication between development team and stakeholders.

---

## 1. Spec Directory Structure

Every feature specification MUST have the following canonical structure:

```
specs/[NNN]-[feature-slug]/
├── spec.md              (MANDATORY: User scenarios, requirements, acceptance criteria)
├── plan.md              (MANDATORY: Implementation strategy, milestones, dependencies)
├── tasks.md             (MANDATORY: Actionable task list with effort estimates)
├── research.md          (OPTIONAL: Background, design decisions, alternatives considered)
├── data-model.md        (CONDITIONAL: Required if spec involves data structure changes)
├── checklists/
│   └── progress.md      (OPTIONAL: Checkpoint progress tracking)
├── contracts/           (OPTIONAL: Integration contracts, API signatures)
└── [other-docs]/        (OPTIONAL: Domain-specific docs, quickstart.md, ANALYSIS.md, etc.)
```

### Mandatory Files

#### 1.1 spec.md
- **Purpose:** Define user scenarios, requirements, success criteria
- **Template:** Use `.specify/templates/spec-template.md`
- **Key sections:**
  - User Scenarios & Testing (prioritized, independently testable stories)
  - Requirements (functional, edge cases, key entities)
  - Success Criteria (measurable outcomes)
  - Assumptions
- **Lifecycle:** Created at spec inception; refined during design phases; finalized before implementation

#### 1.2 plan.md
- **Purpose:** Document implementation strategy, phases, milestones, blockers
- **Template:** Use `.specify/templates/plan-template.md`
- **Key sections:**
  - High-level approach
  - Phase breakdown (if multi-phase)
  - Milestones with success metrics
  - Known blockers and dependencies
  - Risk assessment
- **Lifecycle:** Drafted after spec approval; updated as phases complete

#### 1.3 tasks.md
- **Purpose:** Actionable task list for implementation team
- **Template:** Use `.specify/templates/tasks-template.md`
- **Key sections:**
  - Task breakdown (one logical unit per task)
  - Effort estimates (t-shirt sizes or hours)
  - Dependencies (task X blocks task Y)
  - Acceptance criteria (how to verify task completion)
  - Test scenarios (what to test; where to test)
- **Lifecycle:** Created from plan.md; updated as work progresses; closed when all tasks complete

### Conditional Files

#### 2.1 data-model.md
- **Trigger:** Required if spec changes data structures, adjacency formats, I/O format, or introduces new key entities
- **Content:**
  - Entity definitions (attributes, types, invariants)
  - Relationships (how entities connect)
  - Access patterns (how data is queried/updated)
  - Backward compatibility (what breaks, what's preserved)
- **Example:** specs/001-admesh-domains-io-parity/data-model.md

#### 2.2 research.md
- **Trigger:** Recommended if spec involves algorithmic changes, novel techniques, or significant design decisions
- **Content:**
  - Related work (published papers, prior art)
  - Design alternatives considered (pros/cons)
  - Decision rationale
  - References and citations
- **Example:** specs/002-fem-smoother-quads/research.md

### Optional Files

#### 3.1 contracts/
- **Purpose:** Explicit API signatures, integration contracts, communication protocols
- **Content:** JSON/YAML schemas, type signatures, message formats
- **Example:** specs/001-admesh-domains-io-parity/contracts/

#### 3.2 checklists/
- **Purpose:** Progress tracking across milestones
- **Files:**
  - `progress.md` — Status dashboard (Phase 1: design ✅, implementation 🔄, testing ⏳)
  - Phase-specific checklists (e.g., `review-checklist.md`, `release-checklist.md`)
- **Example:** specs/001-admesh-domains-io-parity/checklists/progress.md

#### 3.3 Domain-Specific Files
- **quickstart.md** — Quick-start guide for users/developers
- **ANALYSIS.md** — Detailed technical analysis or failure analysis
- **[component]-design.md** — Component-level design details
- **[component]-api.md** — Detailed API reference for a component

---

## 2. File Content Standards

### 2.1 Markdown Style
- **Format:** GitHub-flavored Markdown (GFM)
- **Headings:** H1 for document title, H2 for top-level sections, H3+ for subsections
- **Links:** Relative paths for internal references; absolute URLs for external
- **Code blocks:** Language-tagged fences (```python, ```json, etc.)
- **Tables:** Use GFM table syntax for structured data

### 2.2 Metadata Headers
Every spec MUST include metadata header in spec.md:

```markdown
# Feature Specification: [FEATURE NAME]

**Feature Branch**: `[###-feature-name]`  
**Created**: [DATE]  
**Status**: [Draft | In Review | Approved | In Progress | Complete]  
**Issue**: [#NNN]  
**Input**: [User description or requirements source]
```

### 2.3 Cross-References
- Reference related specs: `See specs/###-related-feature/`
- Reference GitHub issues: `#NNN` auto-links on GitHub
- Reference sections: `See § 3.1 [Section Name]`
- Reference requirements: `FR-001`, `SC-003` (from spec.md)

### 2.4 Decision Documentation
Major decisions MUST be documented with rationale:

```markdown
## Design Decision: [DECISION TITLE]

**Decision:** [What was decided]  
**Rationale:** [Why this choice over alternatives]  
**Alternatives Considered:**
- Option A: [pros/cons]
- Option B: [pros/cons]  

**Owner:** [Who made the decision]  
**Date:** [YYYY-MM-DD]  
**Reversible:** [Yes/No] — [Explanation if irreversible]
```

---

## 3. Spec Lifecycle

### Phase 1: Inception (Days 1-3)
- **Input:** GitHub issue with user request/problem statement
- **Deliverables:**
  - spec.md (draft) — user scenarios + requirements
  - plan.md (high-level) — approach sketch
- **Gate:** Issue author + maintainer approve direction
- **Status:** `Draft` in spec.md metadata

### Phase 2: Design Review (Days 4-7)
- **Activities:**
  - Flesh out plan.md with phases, dependencies, risks
  - Create data-model.md if data structures change
  - Add research.md if algorithmic/novel
  - Solicit feedback on design
- **Gate:** 3+ reviewers approve; concerns addressed or documented
- **Status:** `In Review` → `Approved`

### Phase 3: Implementation Planning (Days 8-10)
- **Activities:**
  - Create tasks.md from plan.md
  - Estimate effort per task
  - Set up checklists/ for progress tracking
  - Create contracts/ for integrations
- **Gate:** Team agrees on task breakdown; effort estimates reasonable
- **Status:** `In Progress`

### Phase 4: Implementation (Days 11+)
- **Activities:**
  - Execute tasks.md in order
  - Update checklists/progress.md regularly
  - Document blockers, decisions, changes in spec
  - Add ANALYSIS.md if design changes discovered mid-flight
- **Gate:** All tasks complete; all acceptance criteria met; all tests pass
- **Status:** Update to `Complete` when shipped

### Phase 5: Post-Mortem (Day after ship)
- **Activity:** Update research.md with lessons learned, update constitution if needed
- **Deliverable:** Lessons learned section added to appropriate file

---

## 4. Quality Standards for Specs

### 4.1 Completeness Checklist

Every spec at `Approved` status MUST have:

- [ ] spec.md with ≥2 user stories, each independently testable
- [ ] All requirements labeled (FR-###, SC-###, etc.)
- [ ] Success criteria are measurable (not "good", "fast", "correct")
- [ ] Assumptions documented
- [ ] Edge cases identified
- [ ] plan.md with phase breakdown + milestones
- [ ] tasks.md with effort estimates and acceptance criteria
- [ ] data-model.md (if data changes) OR explicit note "No data model changes"
- [ ] research.md (if algorithmic) OR explicit note "No design decision needed"
- [ ] Cross-references to related specs/issues/docs

### 4.2 Consistency Requirements

- No conflicting requirements (FR-## vs FR-## must align)
- Tasks in tasks.md must cover all FR/SC items
- Success criteria in spec.md must align with tasks.md acceptance criteria
- Assumptions in spec.md must be tested or documented in plan.md
- Milestones in plan.md must map to task groupings in tasks.md

### 4.3 Style & Language

- **Tone:** Professional, precise, assume technical audience
- **Tense:** Present tense for requirements (e.g., "System MUST") and conditionals (e.g., "If X, then Y")
- **Clarity:** Every sentence should be unambiguous; if unsure, add example
- **Completeness:** Avoid "TBD" or "TK" unless explicitly tracked; prefer `[NEEDS CLARIFICATION: ...]`

---

## 5. Governance

### 5.1 Authority & Approval

| Decision | Authority | Process |
|----------|-----------|---------|
| Create new spec | Any contributor | Open GitHub issue; link to specs/ dir |
| Approve spec.md | Issue author + Maintainer | Review comment with `:+1:` or `Approved` label |
| Approve plan.md | Maintainer | Review for feasibility, dependencies, risks |
| Approve tasks.md | Implementation lead | Estimate effort; verify all FR/SC covered |
| Change scope mid-flight | Maintainer + downstream stakeholders | Update spec.md with new requirement; replan if needed |

### 5.2 Breaking Changes in Mid-Flight

If implementation discovers a breaking change needed:

1. Document the change in spec.md (add `## Mid-Flight Scope Changes` section)
2. Flag for maintainer review (GitHub comment with `@maintainer`)
3. Get approval before implementing breaking change
4. Update plan.md, tasks.md, and affected downstream specs

### 5.3 Deprecation of Specs

If a spec is abandoned or superseded:

1. Mark as `Deprecated` in spec.md metadata
2. Link to replacement spec (if applicable)
3. Leave in version control for historical reference
4. Do NOT delete

---

## 6. Integration with CHILmesh Governance

### 6.1 Alignment with constitution.md

CHILmesh's global constitution (`.planning/constitution.md`) takes precedence over speckit-constitution. If conflict:

- **Global governance principles** (testing, backward compatibility, peer review) → apply as written
- **Spec structure/format** (file names, metadata) → follow speckit-constitution

### 6.2 Alignment with CLAUDE.md

CHILmesh's development guide (`.claude/CLAUDE.md`) dictates:

- **Branch policy:** All work on `daily-issue-fixing` branch (not `spec-branch` variable in spec.md)
- **Testing baseline:** All tests must pass before committing spec changes
- **Code review:** Peer review required per constitution.md § 3.3

### 6.3 Alignment with project_plan.md

The overall project roadmap (`.planning/project_plan.md`) defines phase milestones. Specs must:

- Reference which project phase they belong to
- Document dependencies on other phases
- Confirm no backward compatibility breaking before Phase N ships

---

## 7. Templates & Tools

### 7.1 Creating a New Spec

```bash
# From repo root:
bash .specify/scripts/bash/create-new-feature.sh "001-my-feature-name"
# OR manually:
mkdir -p specs/001-my-feature-name/checklists
cp .specify/templates/spec-template.md specs/001-my-feature-name/spec.md
cp .specify/templates/plan-template.md specs/001-my-feature-name/plan.md
cp .specify/templates/tasks-template.md specs/001-my-feature-name/tasks.md
```

### 7.2 Checking Spec Completeness

Maintainers SHOULD run before merging:

```bash
# Check that all specs in /specs have mandatory files:
for spec in specs/*/; do
  test -f "$spec/spec.md" || echo "MISSING: $spec/spec.md"
  test -f "$spec/plan.md" || echo "MISSING: $spec/plan.md"
  test -f "$spec/tasks.md" || echo "MISSING: $spec/tasks.md"
done
```

### 7.3 Lint Spec Markdown

- [ ] TODO: Add spec linting script to `.specify/scripts/` (optional enhancement)

---

## 8. Common Spec Patterns

### Pattern 1: Bug Fix Spec

Use when spec is fixing a bug or regression.

**Key differences from feature spec:**
- User Story 1: "X used to work; now it breaks"  
- Requirements: "System MUST restore behavior to [version]"  
- Success Criteria: "Regression test passes; baseline behavior restored"  
- data-model.md: Usually N/A unless bug exposed data structure issue  
- research.md: Usually N/A unless root cause analysis needed

### Pattern 2: Refactoring Spec

Use when spec is internal cleanup without user-visible changes.

**Key differences:**
- User Story 1: "[Maintainer] wants to refactor [component] for clarity"  
- Requirements: "System MUST preserve all public APIs; internal structure may change"  
- Success Criteria: "All tests pass; code complexity reduced (per metric); performance ≥ baseline"  
- data-model.md: REQUIRED (internal structure change)  
- research.md: Should document why refactoring is needed

### Pattern 3: Performance Optimization Spec

Use when spec improves speed/memory without changing functionality.

**Key differences:**
- User Story 1: "[User] has [problem] that requires faster [operation]"  
- Requirements: Quantified (e.g., "System MUST reduce latency from 10s to <1s on Block_O mesh")  
- Success Criteria: Benchmark before/after on all fixtures  
- research.md: REQUIRED (document algorithm change, tradeoffs)  
- data-model.md: Maybe (if internal caching structure added)

### Pattern 4: Integration Spec

Use when spec bridges CHILmesh with downstream projects (MADMESHR, ADMESH, etc.).

**Key differences:**
- User Story 1: "[Downstream project] needs [capability] to [goal]"  
- Requirements: Include downstream stakeholder feedback  
- Success Criteria: Integration tests pass; downstream PR merges  
- contracts/: REQUIRED (explicit API, message formats, version pins)  
- research.md: Document integration architecture

---

## 9. Document Control

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2026-05-14 | Claude Code Session | Initial constitution: mandatory files (spec/plan/tasks), conditional (data-model/research), optional (contracts/checklists), lifecycle, approval gates, templates, quality standards |

---

## Appendix A: Example Compliant Spec Layout

```
specs/003-example-feature/
├── spec.md
│   └── User Stories (3 P1 stories, each independently testable)
│   └── Requirements (FR-001 through FR-008)
│   └── Success Criteria (SC-001 through SC-004)
│   └── Assumptions
│   └── Metadata: Status=Approved
│
├── plan.md
│   └── Phase 1: Design & Validation (Days 1-5)
│   └── Phase 2: Implementation (Days 6-15)
│   └── Phase 3: Testing & Review (Days 16-20)
│   └── Blockers: Depends on PR #99 landing first
│
├── tasks.md
│   ├── [Phase 1] Design API contract (2h) → acceptance: contracts/ reviewed
│   ├── [Phase 2] Implement core logic (8h) → acceptance: tests pass
│   ├── [Phase 2] Add integration tests (4h) → acceptance: ≥85% coverage
│   └── [Phase 3] Documentation (3h) → acceptance: docstrings reviewed
│
├── data-model.md (if data structures change)
│   └── New Entity: MeshBoundary (attributes, invariants)
│   └── Updated Entity: Mesh (added boundary field)
│   └── Access patterns: How to query/update
│   └── Backward compat: Old code using mesh.adjacencies still works
│
├── research.md (if algorithmic changes)
│   └── Related work: Cite papers on boundary handling
│   └── Algorithm choice: Why option B over option A
│   └── Design alternatives: What was considered and rejected
│
├── contracts/
│   ├── api-schema.json (explicit types/signatures)
│   └── example-usage.py (integration example)
│
└── checklists/
    └── progress.md
        - Phase 1: ✅ Design complete (2026-05-15)
        - Phase 2: 🔄 Implementation in progress (50%)
        - Phase 3: ⏳ Testing not started
```

---

**This constitution is a living document. Amendments follow the process in § 5 of `.planning/constitution.md`.**
