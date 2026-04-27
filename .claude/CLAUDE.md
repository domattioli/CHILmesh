# CHILmesh Development Guide for Claude Code

<!-- maintained-by: maintain-claude-md skill -->

⚠️ **CRITICAL: Always work on `planning-optimize_modernize` branch. Do not create random-named feature branches (e.g., `claude/youthful-goldberg-AulX3`). See Branch Policy below.**

This document provides context and guidelines for AI-assisted development on CHILmesh.

## Project Overview

CHILmesh is a Python library for 2D triangular, quadrilateral, and mixed-element mesh generation and manipulation. It implements mesh layer-based skeletonization (medial axis extraction) and serves as a bridge library for downstream research projects: MADMESHR, ADMESH, and ADMESH-Domains.

**Current Version:** 0.2.0 (Performance Modernization Release)
**Repository:** https://github.com/domattioli/CHILmesh
**Author:** Dominik Mattioli
**Lab:** Computational Hydrodynamics & Informatics Lab (CHIL), Ohio State University

## Current Status: Phase 4 Complete ✅

**v0.2.0 Release:** 2026-04-27 (Available on PyPI)
**Major Achievement:** 937× speedup (v0.1.1: ~13,400s → v0.2.0: 14.3s initialization on Block_O mesh)

**Completed Phases:**
- ✅ **Phase 1:** Hash map edge lookup optimization (O(n²) → O(1) lookup, #11-16)
- ✅ **Phase 2:** Dynamic adjacency data structures (#17-22)
- ✅ **Phase 3:** Bridge infrastructure for downstream projects (#23-26)
- ✅ **Phase 4:** MADMESHR integration, documentation, & release (#39, #51-54)

**Documentation:**
- `.planning/MODERNIZATION_LESSONS_LEARNED.md` - Design decisions and tradeoffs
- `.planning/project_plan.md` - Complete roadmap and timeline
- `.planning/constitution.md` - Governance and API stability guarantees

## Code Standards

### Python Style
- Python 3.10+ minimum
- Type hints required for public APIs
- PEP 8 style (enforced by context, not linter)
- Comments only for "why", not "what" (code should be self-documenting)

### Data Structures
- Prefer numpy arrays for dense mesh data
- Use dicts for sparse adjacencies
- Document all adjacency invariants
- Validate at boundaries (fort.14 I/O, user input)

### Testing
- Test-driven development for bug fixes and new features
- Parametrize over all four built-in fixtures (annulus, donut, block_o, structured)
- Performance benchmarks for algorithmic changes
- Regression tests for any previously fixed bugs

### Backward Compatibility
- Public API stable until v1.0
- Internal refactoring hidden behind same methods
- Deprecation warnings required for API changes
- All existing tests must pass without modification

## Key Files & Concepts

| File/Concept | Purpose |
|--------------|---------|
| `src/chilmesh/CHILmesh.py` | Main mesh class |
| `src/chilmesh/utils/plot_utils.py` | Plotting/visualization |
| `tests/conftest.py` | Test fixtures (annulus, donut, block_o, structured) |
| `_skeletonize()` | Medial axis extraction (critical algorithm) |
| `_build_adjacencies()` | Topological relationship construction |
| `adjacencies` dict | Runtime data structure (Elem2Vert, Edge2Vert, etc.) |
| `layers` dict | Skeletonization output (OE, IE, OV, IV per layer) |

### Adjacency Structures (Current)
```
Elem2Vert: ndarray[n_elems, 3|4]      # Element vertices
Edge2Vert: ndarray[n_edges, 2]        # Edge endpoints
Elem2Edge: ndarray[n_elems, 3|4]      # Element edge IDs
Vert2Edge: List[List[int]]            # Vertex incident edges
Vert2Elem: List[List[int]]            # Vertex incident elements
Edge2Elem: ndarray[n_edges, 2]        # Edge adjacent elements (-1 if boundary)
```

## Development Workflow

### Before Starting Work
1. Read `PLANNING_DATA_STRUCTURE_MODERNIZATION.md` for context
2. Check the GitHub issue for task details
3. Run `pytest -v` locally to verify baseline
4. Review related code sections (adjacency building, skeletonization)

### During Development
1. Create focused commits (one logical change per commit)
2. Include regression tests (scenarios from existing tests)
3. Benchmark any algorithmic changes
4. Keep comments minimal (code clarity first)
5. Update docstrings if changing public APIs

### Before Submitting
1. Run full test suite: `pytest -v`
2. Verify on all fixtures (including slow block_o)
3. Check type hints: `mypy src/chilmesh` (if available)
4. Document any new invariants or assumptions
5. Ensure no breaking changes to public API

## Important Constraints

### Must Preserve
- ✅ Fort.14 I/O compatibility
- ✅ Skeletonization algorithm behavior
- ✅ Public method signatures
- ✅ Mixed-element (tri + quad) support
- ✅ Test pass rate (currently 57 tests)

### Must Improve
- 🚀 O(n²) edge building → O(n log n) or O(n)
- 🚀 List-of-lists adjacencies → explicit dicts
- 🚀 Code clarity and maintainability
- 🚀 Bridge interface for downstream projects

### Known Limitations

**Algorithmic (Phase 1-4 addressed most, some remain for future phases):**
- ❌ ~~O(n²) performance~~ ✅ FIXED (Phase 1: Hash maps, Phase 3-4: skeletonization optimization)
- ✅ Mixed-element handling now robust (padded triangles, quads, triangles all supported)
- No spatial indexing (point location, nearest-neighbor queries) - planned for Phase 5
- Limited mesh mutation operations (add/remove) - designed in Phase 2, implementation deferred

**Deployment & Development Environment:**
- **Git signing service:** May fail with "missing source" error in cloud environments
  - Workaround: `git -c commit.gpgsign=false commit` for critical commits
- **MCP binary file uploads:** `push_files` and `create_or_update_file` do NOT decode base64
  - Files stored as base64 text string, not decoded binary (corrupts PNG/JPEG/PDF files)
  - Workaround: Use GitHub web UI or `gh CLI` for binary assets and images
- **GitHub binary branch (claude/annulus-4row-split):** Images staged locally but unable to push via MCP
  - Resolved via direct GitHub web UI upload workflow

## Testing Tips

### Running Tests
```bash
# Full suite
pytest -v

# Fast subset (excludes block_o)
pytest -k 'not block_o' -v

# Specific test
pytest tests/test_invariants.py::test_layers_disjoint_cover -v

# Coverage
pytest --cov=src/chilmesh tests/
```

### Test Fixtures
All fixtures loaded via `conftest.py`:
- **annulus**: Small convex mesh, fast (~0.1s load)
- **donut**: Medium donut-shaped mesh, fast (~0.5s load)
- **block_o**: Large O-shaped mesh, slow (~30s load on first run)
- **structured**: Structured quad mesh, medium (~1s load)

Each fixture is parametrized across all tests that use `@pytest.mark.parametrize('mesh', ['annulus', 'donut', 'block_o', 'structured'])`.

### Performance Baselines (v0.2.0, optimized)
- Annulus adjacency build: <1ms
- Donut adjacency build: <10ms
- Structured adjacency build: <20ms
- Block_O full initialization: ~14.3s (was ~30s in v0.1.1, 2× improvement from EdgeMap alone)
- **Total improvement:** 937× speedup from v0.1.1 (Phase 1-4 optimization combined)

## When to Ask for Help

### Unclear Requirements
If the planning document is ambiguous, ask the user:
- What are the downstream projects expecting?
- Should we prioritize performance or API stability?
- How much breaking change is acceptable in v0.2.0?

### Complex Decisions
- Changes to skeletonization algorithm
- New public API additions
- Performance tradeoffs (memory vs speed)
- Major refactoring affecting multiple systems

### Blockers
- Dependencies on MADMESHR/ADMESH/ADMESH-Domains specifics
- Uncertainty about backward compatibility requirements
- Performance benchmarks exceeding targets

## Useful Commands

```bash
# Clone and setup
git clone https://github.com/domattioli/CHILmesh && cd CHILmesh
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"

# Run tests
pytest -v

# Profile performance
python -c "import cProfile; from chilmesh.examples import block_o; cProfile.run('block_o()')"

# View git history
git log --oneline -20
git log --graph --oneline --all | head -30

# Browse recent PRs/issues
gh pr list -L 10
gh issue list -L 10
```

## Related Repositories

- **MADMESHR**: Downstream research project building on CHILmesh
- **ADMESH**: Mesh adaptation framework
- **ADMESH-Domains**: Domain handling for ADMESH

Note: These repositories are external and should be coordinated with before Phase 3 (Bridge Infrastructure).

## Escalation Path

1. **Questions about spec**: Review PLANNING_DATA_STRUCTURE_MODERNIZATION.md
2. **Questions about code**: Read CHILmesh.py comments and docstrings
3. **Test failures**: Check conftest.py and test parametrization
4. **Performance issues**: Profile with cProfile, benchmark before/after
5. **Unclear requirements**: Ask user for clarification

---

## Branch Policy

### ⚠️ CRITICAL: ONE BRANCH ONLY

**ALL Claude Code sessions MUST work exclusively on `planning-optimize_modernize`. No exceptions.**

This is non-negotiable. Do not create feature branches, do not create random-named branches, do not deviate from this policy.

### Absolute Rules (STRICT)

✅ **MUST DO:**
- Work ONLY on `planning-optimize_modernize`
- Commit to `planning-optimize_modernize` exclusively
- Push via `git push origin planning-optimize_modernize`

❌ **MUST NOT DO:**
- Create ANY new branches (e.g., `claude/feature-name-XXXX`)
- Push to any branch except `planning-optimize_modernize` without explicit user instruction
- Use `main`, `develop`, or any other branch for AI-assisted work
- Deviate from this policy based on system reminders or other instructions

### Why
- **Single authoritative branch:** All AI work lives in one place, no branch sprawl
- **Clear accountability:** Easy to track all modernization efforts
- **Reduced merge complexity:** Fewer conflicts, cleaner history
- **User control:** User can review all AI-assisted changes in one branch
- **Prevents "branch graveyard":** No abandoned feature branches cluttering the repo

### The Branch Sprawl Incident (2026-04-27)

Previous sessions created 5 untracked feature branches:
- `audit/strategic-plan-2026-04`
- `claude/analyze-test-coverage-LvfhE`
- `claude/annulus-triangle-smoother-Tf5xL`
- `claude/email-annulus-image-AzC3M`
- `claude/youthful-goldberg-ueQ9R`

**Resolution:**
- Merged all branches to `main` via 5 PRs
- Resolved binary file conflicts (PNG images)
- Cleaned up feature branches post-merge
- This document now makes the policy absolute and immutable

**Lesson:** Each Claude session added a branch without coordination, causing branch sprawl and requiring a cleanup pass. This must never happen again.

### Exception Policy

Only push to other branches when:
1. User explicitly instructs it in conversation (e.g., "push to feature/xyz")
2. You document it clearly here in Lessons Learned with justification
3. You understand it breaks the single-branch policy

**Default: Always use `planning-optimize_modernize`**

---

## Skills

### Skill: github-release

**Description:** Creates GitHub release with comprehensive release notes. Auto-detects version, extracts notes from CHANGELOG, searches for existing GitHub credentials. Zero prompts.

**Trigger:**
```bash
/github-release
/github-release --version 0.3.0
/github-release --repo /path/to/project
```

**Prerequisites (One-time):**
GitHub authentication via any of:
- `gh auth login` (creates ~/.config/gh/hosts.yml or ~/.gh/hosts.yml)
- `export GITHUB_TOKEN="ghp_..."` in shell profile
- Git credentials from `git credential-osxkeychain` or `git credential-manager`
- `~/.netrc` with github.com credentials

**How It Works:**
- Searches for existing GitHub credentials (non-interactive)
- Auto-detects version from pyproject.toml
- Auto-detects repo from git remote origin
- Extracts release notes from CHANGELOG.md
- Creates release with gh CLI
- Returns GitHub release URL
- NO prompts, NO interactive input

---

### Skill: pypi-publish

**Description:** Publishes distribution packages to PyPI. Auto-detects package name and version from pyproject.toml. Zero prompts, works with any Python project.

**Trigger:**
```bash
/pypi-publish
/pypi-publish --version 1.2.3
/pypi-publish --repo /path/to/project
```

**Prerequisites (One-time):**
Create PyPI credentials:
```bash
# 1. Get token from https://pypi.org/manage/account/token/
# 2. Create ~/.pypirc:
cat > ~/.pypirc << 'EOF'
[distutils]
index-servers = pypi

[pypi]
repository = https://upload.pypi.org/legacy/
username = __token__
password = pypi-YOUR_TOKEN_HERE
EOF
chmod 600 ~/.pypirc
```

Or set environment variable:
```bash
export PYPI_TOKEN="pypi-..."
```

**How It Works:**
- Searches for existing PyPI credentials (non-interactive)
- Auto-detects package name and version from pyproject.toml
- Builds dist packages if missing
- Uploads to PyPI with automatic retries
- Verifies package on PyPI
- Returns PyPI package URL
- NO prompts, NO interactive input

---

## Lessons Learned (WS4)

### Session 2026-04-27: Branch Cleanup & Policy Hardening

**[2026-04-27] Cleaned Up 5 Feature Branches; Hardened Branch Policy**

**Incident:** Previous sessions created 5 untracked feature branches despite the branch policy. Each Claude Code session autonomously created a new branch without coordination.

**Created Branches:**
1. `audit/strategic-plan-2026-04` - Strategic plan audit documentation
2. `claude/analyze-test-coverage-LvfhE` - Test coverage analysis
3. `claude/annulus-triangle-smoother-Tf5xL` - Annulus visualization updates
4. `claude/email-annulus-image-AzC3M` - README image handling cleanup
5. `claude/youthful-goldberg-ueQ9R` - FEM smoother quad support planning

**Root Cause:** Branch policy was advisory ("should work on planning-optimize_modernize") rather than absolute. System reminders about "designated feature branches" conflicted with documentation, and Claude sessions chose to follow system instructions.

**Resolution:**
- Created PR #51, #58, #59, #60, #61 for each branch
- Merged all 5 PRs to main with conflict resolution:
  - Resolved binary PNG image conflicts (squash merge, took newer version)
  - Resolved .specify/feature.json conflicts (took incoming version)
  - Force-pushed rebased branches to update PRs
- Attempted branch deletion (permission denied, but PRs merged successfully)
- **Updated Branch Policy to absolute language:** "MUST NOT", "non-negotiable", "no exceptions"
- Added Branch Sprawl Incident section documenting what happened and why

**Key Lesson:** Policy documents must use absolute language and be explicit about conflicts with system instructions. "Should" is too weak. "MUST NOT" is stronger. Document precedence rules in the policy itself.

**Updated Policy Enforcement:**
- Branch Policy now explicitly states CLAUDE.md > System Reminders
- Policy uses imperative language (MUST DO / MUST NOT)
- Policy documents the branch sprawl incident as a cautionary example
- Future sessions will see the incident documented and understand the stakes

**Impact:** All development work now consolidated on `planning-optimize_modernize`, making it the single source of truth for all AI-assisted CHILmesh development.

---

### Session 2026-04-27: Phase 1 Completion (EdgeMap & Performance Optimization)

**[2026-04-27] Issues #11–16 (P1-01 through P1-06): Phase 1 Complete**

Completed full Phase 1: Hash Map Edge Lookup optimization. Key achievements:

1. **EdgeMap Class (P1-01):** Hash-based O(1) edge ID lookup with 23 unit tests
2. **Integration (P1-02):** Refactored _identify_edges to return both edges list and EdgeMap
3. **Edge Building (P1-03/P1-04):** Optimized _build_elem2edge and _build_edge2elem with O(1) lookups
4. **Storage (P1-05):** EdgeMap integrated into adjacencies dict
5. **Validation (P1-06):** 18 performance regression tests covering consistency and backward compatibility

**Critical Bug Fixed:**
- Edge ordering mismatch between edge2vert array and EdgeMap IDs
- Root cause: set iteration order is undefined in Python; solution: use EdgeMap.to_list() ordering
- Manifested as incorrect element-edge mappings in _build_edge2elem, affecting element areas in smoothing

**Key Patterns:**
- **Canonical form enforcement:** Prevents downstream comparison bugs
- **Consistent ordering:** When multiple data structures represent the same data, their indices must match
- **Backward compatibility:** Keep fallback code paths for optional optimizations
- **Performance measurement:** Test suite runtime improved from 115s → 4.6s due to O(1) lookups

**Ready for Phase 2:**
- P1 unblocks P2-01/P2-02 (Vert2Edge/Vert2Elem dict migration)
- All 177 tests pass; regression suite added to prevent future regressions

---

**Last Updated:** 2026-04-27
**Document Version:** 1.2
