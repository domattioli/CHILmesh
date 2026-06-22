# CHILmesh semver integrity — gate + `write_fort14` break

**Date:** 2026-06-20
**Status:** Design approved; implementation local-only (no push/tag/PyPI until operator says so)
**Repo:** CHILmesh (primary). valence consumer update = tracked follow-on (separate repo/PR).

## Background

`write_fort14`'s public signature changed from a free-function form
`write_fort14(path, points, elements, name)` (raw arrays) to an object form
`write_fort14(mesh, filename)` in PR #184 (commit `6080077`, 2026-06-08, "feat: Gmsh .msh I/O").
It is exported in `__all__` (`src/chilmesh/CHILmesh.py:17`, `src/chilmesh/__init__.py`).

This breaking change shipped as **1.2.2** — a patch bump — with **no CHANGELOG entry**
for the break (the `[1.2.2]` section records only registry renames). The project claims
SemVer adherence (CHANGELOG header). So: a silent, undocumented breaking change in a patch release.

Downstream impact: `valence-domains`' generator calls the **old** 4-arg form
(`src/valence_domains/generator.py:329`), pinned `chilmesh>=1.2.1`. In any environment with an
editable/newer CHILmesh (the `[seg]` mesh-refinement venv), the generator + CLI tests fail with
`TypeError: write_fort14() takes 2 positional arguments but 4 were given` (14 failures). Published-wheel
CI for valence stays green only because it resolves an older CHILmesh.

## Goals

1. **Prevent recurrence** — CI gate that fails when the public Python API breaks without a major bump.
2. **Own the existing break** — make the `write_fort14` change semver-honest (keep the new API; correct version + changelog).
3. **(Follow-on) Restore the consumer** — update valence's generator to the new API + re-pin.

Non-goal: covering the Rust-core (`chilmesh_core`, pyo3) ABI. The gate inspects the Python surface
(`src/chilmesh`) only — sufficient for the class of break that occurred here.

## Component 1 — Semver CI gate

New job `api-semver-gate` in `.github/workflows/python-package.yml`, runs on `pull_request`.

**Tool:** [`griffe`](https://mkdocstrings.github.io/griffe/) `check` — purpose-built Python API-diff;
static (no install/Rust build), can load an old git ref directly.

**Steps:**
1. `actions/checkout@v4` with `fetch-depth: 0` (needs tags + history).
2. Set up Python 3.12; `pip install griffe`.
3. `baseline=$(git describe --tags --abbrev=0)`.
4. Run `griffe check chilmesh -s src -a "$baseline"`; capture exit code (nonzero ⇒ breaking public-API change vs baseline). Don't let a nonzero abort the step — capture and branch in step 6.
5. Extract versions: `OLD=$(git show "$baseline":pyproject.toml | grep -m1 '^version')`, `NEW=$(grep -m1 '^version' pyproject.toml)`; parse the major component of each.
6. **Gate decision:**
   - breaking **and** `major(NEW) == major(OLD)` → **FAIL**: print the griffe diff + "Public API broke without a major bump — bump major (`pyproject.toml`) or revert."
   - breaking **and** major bumped → pass.
   - not breaking → pass.

YAGNI: not enforcing minor-bump-on-addition; the failure mode that bit us is breaking-without-major.

**Placement:** sibling job to `test`; independent (no `needs`). PR-only (the break must be caught before merge).

## Component 2 — Own the `write_fort14` break (CHILmesh)

The new `write_fort14(mesh, filename)` is the intended design (object form, emits NOPE/NBOU boundary
sections) — **keep it**. Make semver honest:

- `pyproject.toml`: `version` `1.2.2` → **`2.0.0`**.
- `CHANGELOG.md`: new `[2.0.0]` section with a **BREAKING** note documenting
  `write_fort14(path, points, elements, name)` → `write_fort14(mesh, filename)`, and stating the
  change actually landed in 1.2.2 (#184) without a major bump; 2.0.0 corrects the semver record.
- **Deferred / operator-gated:** tag `v2.0.0` + GitHub release → `publish-pypi.yml`. **Not executed**
  while changes are local-only. Until 2.0.0 is on PyPI, Component 3's re-pin can't go green — expected.

## Component 3 — valence consumer update (follow-on, separate repo/PR)

Tracked, not part of this spec's implementation. Ordering: **CHILmesh 2.0.0 on PyPI first**, then valence.

- `_write_fort14_via_chilmesh` (`src/valence_domains/generator.py`): build
  `CHILmesh(connectivity=elements, points=points, grid_name=mesh_name)` then `write_fort14(mesh, tmp_path)`.
  (Constructor: `CHILmesh.py:243` — `connectivity` first, `points` second.)
- Re-pin `chilmesh>=2.0.0` (3 sites in `pyproject.toml`).
- Add a contract test: import + call `chilmesh.write_fort14` exactly as the generator does, so future
  producer drift fails valence CI immediately regardless of the producer's version hygiene.
- valence has its own DomI-drift hard-stop governance — its PR follows valence's branch discipline.

## Testing

- **Component 1:** verify the gate locally before committing — run the griffe command against the
  current tree vs the latest tag and confirm exit semantics; sanity-check the major-parse + branch logic.
  (A red/green fixture pair is optional; primary signal is the real `write_fort14` break showing as breaking.)
- **Component 2:** no code/test change — version + changelog only. Existing suite must stay green.
- **Component 3 (follow-on):** the new contract test + existing generator/CLI tests pass against
  CHILmesh 2.0.0.

## Risks / notes

- griffe sees only `src/chilmesh` Python surface — Rust ABI breaks invisible (accepted, see Non-goal).
- `git describe --tags` needs tags fetched (`fetch-depth: 0`); if no tag exists yet the gate must
  no-op gracefully (skip when `git describe` fails).
- Local-only constraint means the PyPI publish (and therefore valence's green re-pin) is deferred;
  the spec is complete but Components 2-publish and 3 land later.
