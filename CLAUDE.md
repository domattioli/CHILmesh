# CHILmesh Development Guide for Claude Code

⚠️ **CRITICAL: Always work on `planning-optimize_modernize` branch. Do not create random-named feature branches (e.g., `claude/youthful-goldberg-AulX3`). See Branch Policy below.**

This document provides context and guidelines for AI-assisted development on CHILmesh.

## Project Overview

CHILmesh is a Python library for 2D triangular, quadrilateral, and mixed-element mesh generation and manipulation. It implements mesh layer-based skeletonization (medial axis extraction) and serves as a bridge library for downstream research projects: MADMESHR, ADMESH, and ADMESH-Domains.

**Current Version:** 0.1.1 (Alpha)
**Repository:** https://github.com/domattioli/CHILmesh
**Author:** Dominik Mattioli
**Lab:** Computational Hydrodynamics & Informatics Lab (CHIL), Ohio State University

## Current Initiative: Data Structure Modernization

**Status:** Planning Phase (No implementation yet)
**Branch:** `claude/zen-fermi-NGYbR`
**Documentation:** `PLANNING_DATA_STRUCTURE_MODERNIZATION.md`

The project is modernizing its internal data structures to support:
1. Better performance (eliminate O(n²) edge building)
2. Efficient mesh alteration (add/remove nodes/edges)
3. Clean bridge interface for downstream projects
4. Improved maintainability and testability

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
- O(n²) performance on large meshes (block_o: ~30s build time)
- Mixed-element handling adds code complexity (padded triangles as 4-col arrays)
- No spatial indexing (point location, nearest-neighbor queries)
- Limited mesh mutation operations (no add/remove yet)

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

### Performance Baselines (for reference)
- Annulus adjacency build: <1ms
- Donut adjacency build: <10ms
- Structured adjacency build: <20ms
- Block_O adjacency build: ~30s (O(n²) bottleneck)

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

**ALL Claude Code sessions must work ONLY on the `planning-optimize_modernize` branch.**

### Rules
- ✅ DO: Commit all work to `planning-optimize_modernize`
- ✅ DO: Push to `git push origin planning-optimize_modernize`
- ❌ DON'T: Create random-named branches (e.g., `claude/youthful-goldberg-AulX3`)
- ❌ DON'T: Push to feature branches like `main` or `develop` without explicit permission

### Why
- Single authoritative branch for all AI-assisted work
- Prevents branch proliferation and confusion
- Enables clear history of modernization efforts
- Allows easy rebasing and integration with main

### Exception
- Only push to `main` when explicitly directed by user
- Document any deviations in Lessons Learned

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
