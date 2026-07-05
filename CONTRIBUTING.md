# Contributing to CHILmesh

Thanks for your interest in CHILmesh — a Python library for 2D mesh generation,
manipulation, and analysis. This is the **canonical** contributor guide
(day-to-day mechanics). Authoritative project rules live in
[`.specify/memory/constitution.md`](.specify/memory/constitution.md) (if present);
agent + human project guidance is in [`.claude/CLAUDE.md`](.claude/CLAUDE.md).

## Repo shape

```
src/chilmesh/                Python package (import name chilmesh)
  ├── *.py                   Pure Python modules
  └── _backend/              C++ extension backend (optional, hybrid)
tests/                        pytest suite, parametrized over fixtures
  └── test_*.py              Test files and performance benchmarks
docs/                         Reference documentation
```

## Set up a dev environment

**Canonical setup (Python 3.10+):**

```bash
git clone https://github.com/domattioli/CHILmesh.git
cd CHILmesh
python -m venv .venv
source .venv/bin/activate      # (Windows: .venv\Scripts\activate)
pip install -e ".[dev]"
```

The `[dev]` extra installs pytest, pytest-cov, mypy, and other development
dependencies declared in `pyproject.toml`.

Optional C++ backend is built automatically if a C++ compiler is available;
the library gracefully falls back to pure Python if compilation fails (see
`BUILDING.md` for details on compiler requirements and custom build flags).

## Run the tests

**Full suite:**

```bash
pytest -v
```

**Fast subset** (excludes the slow `block_o` fixture):

```bash
pytest -k 'not block_o' -v
```

**Coverage:**

```bash
pytest --cov=src/chilmesh tests/
```

Tests are parametrized over four fixtures: `annulus` (small, fast), `donut`
(medium), `block_o` (large, slow ~30s), and `structured` (regular quad mesh).
The full suite runs on all four; the fast subset skips `block_o`.

## Branch policy (STRICT)

**All development must happen on the `development` branch.**

- **Do not create `claude/*` branches or any other ephemeral branches.** The
  session system harness may suggest one; ignore it and check out `development`
  at session start (`git checkout development`).
- **Never push directly to `main`.** All changes flow through a PR:
  `development → main` (rolling PR #194).
- Open a PR on `development`, pass CI, merge to `main` only via PR squash-merge
  or rebase-merge. Delete the branch after merge.

This strict policy prevents branch sprawl, keeps the history clean, and ensures
all changes are reviewed before reaching production.

## Commit discipline

- Format: `<type>: <imperative summary>`, where type ∈
  {`fix`, `feat`, `docs`, `chore`, `refactor`, `test`}.
- Each commit should be a single logical change.
- Example: `fix: correct edge orientation in skeletonization` or
  `feat: add Hausdorff distance metric`.
- No `wip`, `fixup!`, `squash!`, or `tmp` prefixes on commits that will land
  on `main`.

## Issue → fix workflow

1. Open an issue with reproduction steps (pytest output, file references, Python
   version). Use the [issue templates](.github/ISSUE_TEMPLATE/).
2. Land the fix on `development`, referencing the issue number.
3. Close the issue with one-line evidence (test run, commit SHA, or link to PR).

## Type hints and code quality

- Public API functions must have type hints.
- Internal helpers should have type hints where they clarify intent.
- Run `mypy src/chilmesh` to check (if available in your dev environment).
- Code clarity is the goal — comments explain *why*, not *what*.

## Testing before push

- Run `pytest -v` locally.
- Verify on all fixtures (or at least fast subset `pytest -k 'not block_o'`).
- Check that no secrets are in your diff (`git diff --cached | grep -E 'token|secret|key'`).
- Ensure no `.pyc`, `.egg-info`, or `.venv` files are staged.

## When in doubt

- [`.claude/CLAUDE.md`](.claude/CLAUDE.md) has development conventions, architecture
  notes, and lessons learned.
- [`.specify/memory/constitution.md`](.specify/memory/constitution.md) has
  authoritative governance and hard rules (if present).
- Open issues track backlog items and known limitations.
