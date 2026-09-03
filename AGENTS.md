# CHILmesh

CHILmesh is a Python library for 2D triangular, quadrilateral, and mixed-element mesh generation and manipulation. It implements mesh layer-based skeletonization (medial axis extraction) and serves as a bridge library for downstream research projects: MADMESHR, ADMESH, and ADMESH-Domains. **Current Version:** 1.2.2 (Production/Stable — Hybrid Python/C++ Release).

## Hard rules

**API stability & public method signatures:**
- Public method signatures are stable; no breaking changes without deprecation warnings
- All existing tests must pass without modification
- Fort.14 I/O compatibility is non-negotiable
- Skeletonization algorithm behavior must be preserved
- Mixed-element (tri + quad) support is required

**Code quality:**
- Type hints required for public APIs
- Comments document the "why", not the "what" (code should be self-documenting)
- Document all adjacency invariants
- Validate inputs at boundaries (fort.14 I/O, user input)

**No secrets in commits:**
- Never commit `.env`, `*token*`, `*secret*`, `*.pem`, `*credentials*`
- No force-push to `main` or shared branches

## Repository layout

```
CHILmesh/
├── src/chilmesh/               # Main package
│   ├── CHILmesh.py            # Main mesh class
│   ├── utils/plot_utils.py    # Plotting/visualization
│   └── ...                     # Supporting modules
├── tests/
│   ├── conftest.py            # Test fixtures (annulus, donut, block_o, structured)
│   ├── test_invariants.py     # Topology & skeletonization tests
│   └── ...                     # Per-module test files
├── .specify/memory/
│   └── constitution.md        # Project governance rules
├── .planning/
│   ├── project_plan.md        # Roadmap and milestones
│   └── MODERNIZATION_LESSONS_LEARNED.md  # Design decisions
├── pyproject.toml             # Package metadata (Python 3.10+, dependencies)
└── README.md                  # Project summary
```

**Key files:**
- `src/chilmesh/CHILmesh.py` — Main mesh class
- `src/chilmesh/utils/plot_utils.py` — Visualization utilities
- `_skeletonize()` — Medial axis extraction (critical algorithm)
- `_build_adjacencies()` — Topological relationship construction
- `adjacencies` dict — Runtime data structure (Elem2Vert, Edge2Vert, Elem2Edge, Vert2Edge, Vert2Elem, Edge2Elem)
- `layers` dict — Skeletonization output (OE, IE, OV, IV per layer)

**Adjacency data structures:**
```
Elem2Vert: ndarray[n_elems, 3|4]      # Element vertices
Edge2Vert: ndarray[n_edges, 2]        # Edge endpoints
Elem2Edge: ndarray[n_elems, 3|4]      # Element edge IDs
Vert2Edge: List[List[int]]            # Vertex incident edges
Vert2Elem: List[List[int]]            # Vertex incident elements
Edge2Elem: ndarray[n_edges, 2]        # Edge adjacent elements (-1 if boundary)
```

## Commands

```bash
# Clone and setup
git clone https://github.com/domattioli/CHILmesh && cd CHILmesh
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"

# Run tests
pytest -v

# Fast subset (excludes slow block_o fixture)
pytest -k 'not block_o' -v

# Specific test
pytest tests/test_invariants.py::test_layers_disjoint_cover -v

# Coverage report
pytest --cov=src/chilmesh tests/

# Profile performance
python -c "import cProfile; from chilmesh.examples import block_o; cProfile.run('block_o()')"

# View git history
git log --oneline -20
git log --graph --oneline --all | head -30

# Browse recent PRs/issues
gh pr list -L 10
gh issue list -L 10
```

## Conventions

### Python style
- Python 3.10+ minimum (required by chilmesh dependency)
- Type hints required for public APIs
- PEP 8 style (enforced by context, not linter)
- Comments only for "why", not "what" (code should be self-documenting)

### Data structures
- Prefer numpy arrays for dense mesh data
- Use dicts for sparse adjacencies
- Document all adjacency invariants
- Validate at boundaries (fort.14 I/O, user input)

### Backward compatibility
- Public API stable until v1.0 (now stable as of v1.0.0)
- Internal refactoring hidden behind same methods
- Deprecation warnings required for API changes
- All existing tests must pass without modification

### Docstring template
```python
def function_name(arg1, arg2):
    """Short summary.

    Longer description of behavior, assumptions, and invariants.

    Parameters
    ----------
    arg1 : type
        Description.
    arg2 : type
        Description.

    Returns
    -------
    result : type
        Description.
    """
```

## Testing

### Test organization
- One test module per main module: `tests/test_<module>.py`
- Parametrize across all four built-in fixtures (annulus, donut, block_o, structured)
- Performance benchmarks for algorithmic changes
- Regression tests for any previously fixed bugs

### Test fixtures
Loaded via `conftest.py`:
- **annulus**: Small convex mesh, fast (~0.1s load)
- **donut**: Medium donut-shaped mesh, fast (~0.5s load)
- **block_o**: Large O-shaped mesh, slow (~30s load on first run)
- **structured**: Structured quad mesh, medium (~1s load)

Each fixture is parametrized across tests that use `@pytest.mark.parametrize('mesh', ['annulus', 'donut', 'block_o', 'structured'])`.

### Performance baselines (v0.2.0 optimized, v1.0.0 with C++ backend)
- Annulus adjacency build: <1ms
- Donut adjacency build: <10ms
- Structured adjacency build: <20ms
- Block_O full initialization: ~14.3s (was ~30s in v0.1.1 before Phase 1 optimization; v1.0.0 C++ backend adds ~15× speedup on top)
- **Historical total improvement:** 937× speedup from v0.1.1 to v0.2.0 (pure-Python Phase 1-4 optimization combined)

### What "tested" means for a PR
- New feature: regression tests + feature-specific tests
- Bug fix: reproduction test (scenario that previously failed) + fix verification
- Refactor: all existing tests pass without modification
- Optimization: performance benchmarks showing improvement + all tests pass

## Branch & commit policy

**Default branch:** `main` (production-ready)
**Development branch:** `development` (AI-session staging branch)

### Branching rules
- Work on `development`, then create PR `development → main`
- No direct pushes to `main`; all changes flow through PR review + CI
- If creating feature branches: use `<type>/<slug>` naming where `<type>` ∈ {fix, feat, docs, chore, refactor, test}
- Delete branches after merge to main
- No long-lived feature branches; PR and merge promptly

### Commit format
```
<type>: <imperative summary>

Optional longer explanation of why this change is needed.

<type> ∈ {fix, feat, docs, chore, refactor, test}
```

### PR rules
- Single-purpose PRs only
- Fill out PR template
- All CI checks must pass (no admin merge or force-push to bypass)
- Squash or rebase merge (not plain merge) to keep main history clean

### Hard stops (refuse outright)
- Force-push to `main` or any shared branch
- Committing `*.env`, `*token*`, `*secret*`, `*.pem`, `*credentials*`
- Merging without CI green

## Edit hygiene

The number of tokens used to edit files is best minimized, all else being equal. Therefore, when it will not affect the end result, opt first for surgical edits rather than rewriting entire existing files.

**Stream timeout prevention:**
- One task per turn; confirm before next
- ≤150 lines per file write; split if longer
- Grep short; use `-l`, `--include` flags to limit output
- If timeout: retry same step, shorter form
- On >20 tool calls: start a fresh session

## Related repositories

| Repo | Purpose |
|------|---------|
| **MADMESHR** | Downstream research project building on CHILmesh |
| **ADMESH** | Mesh adaptation framework using CHILmesh |
| **ADMESH-Domains** | Domain handling for ADMESH |

## Reference docs

**Project governance & roadmap:**
- `.specify/memory/constitution.md` — Canonical project governance rules, principles, API stability contract
- `.planning/project_plan.md` — Roadmap, milestones, "where we are today" status
- `.planning/MODERNIZATION_LESSONS_LEARNED.md` — Design decisions and optimization tradeoffs

**Session handoff & state:**
- `.planning/.continue-here.md` — Most recent session notes and in-flight context

**Test & CI audit:**
- `.planning/TEST-AUDIT.md` — What is/isn't tested, test coverage surface
- `.planning/HOOKS-AUDIT.md` — Pre-commit hook validation

**Repo-local labels:**

These labels have no DomI canonical equivalent and are kept repo-local.

| Label | Meaning |
|---|---|
| `admesh` | Cross-repo dependency: ADMESH |
| `admesh-domains` | Cross-repo dependency: ADMESH-Domains |
| `API` | Public API surface changes |
| `benchmark` | Benchmarking / performance measurement work |
| `bridge` | Phase 3 bridge infrastructure for downstream projects |
| `code-quality` | Code quality and maintainability work |
| `coordination` | Cross-repo coordination with sibling repos |
| `data-structures` | Adjacency / topology data structure design |
| `design` | Algorithmic or API design decisions |
| `domi-sync` | DomI sync-contract issues (opened by `notify-downstream.yml`) |
| `FEM-smoother` | FEM smoother specific work |
| `integration` | Cross-repo integration work |
| `io` | fort.14 / .2dm file I/O |
| `mixed-element` | Tri/quad mixed-element support |
| `optimization` | Runtime optimization work |
| `performance` | Performance measurement and improvement |
| `phase-1` through `phase-5` | Historical milestone tracking (Phases 1–5) |
| `portability` | Cross-platform portability concerns |
| `quality` | Mesh quality metrics |
| `report` | Analysis / benchmark reports |
| `rust-backend` | Rust backend investigations |
| `transfer-candidate` | Candidate for transfer to a sibling repo |
| `types` | Type annotation work |
| `validation` | Validation logic |
