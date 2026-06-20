# CHILmesh Semver Integrity Gate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a CI gate that fails when CHILmesh's public Python API breaks without a major version bump, and own the already-shipped `write_fort14` break by correcting the version + changelog.

**Architecture:** The gate logic lives in a small, unit-testable Python helper (`scripts/api_semver_gate.py`): pure functions (`parse_major`, `decide`) carry the policy and are TDD'd; a thin `main()` shells out to `griffe check` and `git` and is exercised by a real run. A new PR-only workflow job calls the helper. Component 2 is a version + CHANGELOG correction (no code change).

**Tech Stack:** Python 3.12, pytest, [griffe](https://mkdocstrings.github.io/griffe/) (API differ), GitHub Actions.

## Global Constraints

- **Local-only.** No `git push`, no tag, no GitHub release, no PyPI publish. All commits stay on local branch `feature/semver-integrity-gate`.
- Branch discipline: no direct commits to `main`; work on `feature/semver-integrity-gate` (already checked out).
- Gate inspects the **Python** surface (`src/chilmesh`) only — Rust-core (`chilmesh_core`) ABI is out of scope.
- Gate is **pull_request-only** (must catch breaks before merge).
- `griffe` is installed ad hoc in the CI job (`pip install griffe`); it is NOT added to package deps. Unit tests must not import griffe.
- Component 3 (valence consumer update) is OUT OF SCOPE — tracked follow-on in the valence repo, blocked on CHILmesh 2.0.0 reaching PyPI (deferred).

## File Structure

- `scripts/api_semver_gate.py` (create) — gate helper: `parse_major`, `decide`, `main`.
- `tests/test_api_semver_gate.py` (create) — unit tests for the pure logic.
- `.github/workflows/python-package.yml` (modify) — add `api-semver-gate` job.
- `pyproject.toml` (modify) — `version` 1.2.2 → 2.0.0.
- `CHANGELOG.md` (modify) — new `[2.0.0]` BREAKING section.

---

### Task 1: Semver gate helper + CI job

**Files:**
- Create: `scripts/api_semver_gate.py`
- Test: `tests/test_api_semver_gate.py`
- Modify: `.github/workflows/python-package.yml` (add `api-semver-gate` job after the `test` job)

**Interfaces:**
- Produces: `parse_major(version: str) -> int`; `decide(breaking: bool, old_major: int, new_major: int) -> tuple[bool, str]`; `main() -> int` (exit code).

- [ ] **Step 1: Write the failing tests**

Create `tests/test_api_semver_gate.py`:

```python
from __future__ import annotations

import pathlib
import sys

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
import api_semver_gate as g  # noqa: E402


def test_parse_major():
    assert g.parse_major("1.2.2") == 1
    assert g.parse_major("2.0.0") == 2


def test_parse_major_rejects_garbage():
    with pytest.raises(ValueError):
        g.parse_major("not-a-version")


def test_decide_breaking_no_major_bump_fails():
    ok, msg = g.decide(True, 1, 1)
    assert ok is False
    assert "without a major bump" in msg


def test_decide_breaking_with_major_bump_passes():
    ok, _ = g.decide(True, 1, 2)
    assert ok is True


def test_decide_non_breaking_passes():
    ok, _ = g.decide(False, 1, 1)
    assert ok is True
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd ~/Developer/CHILmesh && python -m pytest tests/test_api_semver_gate.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'api_semver_gate'`.

- [ ] **Step 3: Write the helper**

Create `scripts/api_semver_gate.py`:

```python
"""Fail when CHILmesh's public Python API breaks without a major version bump.

Run in CI on pull_request. Compares the public API of ``src/chilmesh`` at the
latest git tag against the working tree via ``griffe check``, then gates the
result against the major component of ``pyproject.toml``'s version.
"""
from __future__ import annotations

import re
import subprocess
import sys
from typing import Optional, Tuple


def parse_major(version: str) -> int:
    """Return the integer major component of a version like ``1.2.2``."""
    m = re.search(r"(\d+)\.(\d+)\.(\d+)", version)
    if not m:
        raise ValueError(f"unparseable version: {version!r}")
    return int(m.group(1))


def decide(breaking: bool, old_major: int, new_major: int) -> Tuple[bool, str]:
    """Gate decision. Returns ``(ok, message)``."""
    if not breaking:
        return True, "No breaking public-API change vs baseline."
    if new_major > old_major:
        return True, f"Breaking change present; major bumped {old_major}->{new_major}. OK."
    return False, (
        f"Public API broke without a major bump (major still {old_major}). "
        "Bump the major in pyproject.toml or revert the breaking change."
    )


def _run(cmd: list) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True)


def _latest_tag() -> Optional[str]:
    p = _run(["git", "describe", "--tags", "--abbrev=0"])
    return p.stdout.strip() if p.returncode == 0 else None


def _version_at(ref: Optional[str]) -> str:
    if ref is None:
        with open("pyproject.toml") as f:
            text = f.read()
    else:
        text = _run(["git", "show", f"{ref}:pyproject.toml"]).stdout
    m = re.search(r'^version\s*=\s*"([^"]+)"', text, re.M)
    if not m:
        raise ValueError("version not found in pyproject.toml")
    return m.group(1)


def main() -> int:
    tag = _latest_tag()
    if tag is None:
        print("No baseline tag found; skipping semver gate.")
        return 0
    griffe = _run(["griffe", "check", "chilmesh", "-s", "src", "-a", tag])
    breaking = griffe.returncode != 0
    if griffe.stdout:
        print(griffe.stdout)
    if griffe.stderr:
        print(griffe.stderr, file=sys.stderr)
    old_major = parse_major(_version_at(tag))
    new_major = parse_major(_version_at(None))
    ok, msg = decide(breaking, old_major, new_major)
    print(msg)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd ~/Developer/CHILmesh && python -m pytest tests/test_api_semver_gate.py -v`
Expected: PASS (5 passed).

- [ ] **Step 5: Add the CI job**

In `.github/workflows/python-package.yml`, add this job after the `test` job (sibling, no `needs`), keeping existing indentation:

```yaml
  api-semver-gate:
    name: API semver gate
    if: github.event_name == 'pull_request'
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
        with:
          fetch-depth: 0

      - name: Set up Python
        uses: actions/setup-python@v5
        with:
          python-version: "3.12"

      - name: Install griffe
        run: |
          python -m pip install --upgrade pip
          python -m pip install griffe

      - name: Check public API vs semver
        run: python scripts/api_semver_gate.py
```

- [ ] **Step 6: Local integration smoke (real griffe run)**

Run:
```bash
cd ~/Developer/CHILmesh
python -m pip install griffe >/dev/null 2>&1
python scripts/api_semver_gate.py; echo "exit=$?"
```
Expected: prints either "No baseline tag found; skipping semver gate." (exit=0) or a griffe diff + a decision line. Either way the script runs to completion without a traceback. (Forward-looking gate — a clean exit here is correct; it catches *future* breaks.)

- [ ] **Step 7: Commit**

```bash
cd ~/Developer/CHILmesh
git add scripts/api_semver_gate.py tests/test_api_semver_gate.py .github/workflows/python-package.yml
git commit -m "ci: add api-semver-gate — fail breaking API change without major bump

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01W1SjM79TLgQEeG9D4NvaD2"
```

---

### Task 2: Own the `write_fort14` break — version + CHANGELOG

**Files:**
- Modify: `pyproject.toml` (`version`)
- Modify: `CHANGELOG.md` (new `[2.0.0]` section at top of the entries)

**Interfaces:**
- Produces: `pyproject.toml` version `2.0.0`; a `[2.0.0]` CHANGELOG section documenting the `write_fort14` break. No code/API change (the new `write_fort14(mesh, filename)` is already the intended API).

- [ ] **Step 1: Bump the version**

In `pyproject.toml`, change the project version line from:
```toml
version = "1.2.2"
```
to:
```toml
version = "2.0.0"
```

- [ ] **Step 2: Add the CHANGELOG section**

In `CHANGELOG.md`, insert immediately above the existing `## [1.2.2] — 2026-06-15` line:

```markdown
## [2.0.0] — 2026-06-20

### Changed — BREAKING

- **`write_fort14` signature** — `write_fort14(path, points, elements, name)`
  (free function over raw arrays) → `write_fort14(mesh, filename)` (CHILmesh
  object form; also writes ADCIRC NOPE/NBOU boundary sections). This change
  actually landed in 1.2.2 ([#184](https://github.com/domattioli/CHILmesh/issues/184))
  without a major bump; 2.0.0 corrects the SemVer record. Callers passing raw
  arrays must construct a `CHILmesh(connectivity=..., points=..., grid_name=...)`
  first, then call `write_fort14(mesh, filename)`.

```

- [ ] **Step 3: Verify the edits**

Run:
```bash
cd ~/Developer/CHILmesh
python -c "import tomllib; assert tomllib.load(open('pyproject.toml','rb'))['project']['version']=='2.0.0'; print('version OK')"
grep -n "## \[2.0.0\]" CHANGELOG.md
grep -n "write_fort14 signature" CHANGELOG.md
```
Expected: prints `version OK`, and two grep hits (the heading + the bullet).

- [ ] **Step 4: Verify the suite still passes**

Run: `cd ~/Developer/CHILmesh && python -m pytest -n auto -m "not slow"`
Expected: PASS (no test depends on the version string; this confirms nothing regressed).

- [ ] **Step 5: Commit**

```bash
cd ~/Developer/CHILmesh
git add pyproject.toml CHANGELOG.md
git commit -m "release: 2.0.0 — own write_fort14 breaking change (semver correction)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01W1SjM79TLgQEeG9D4NvaD2"
```

> **Deferred (operator-gated, NOT in this plan):** tag `v2.0.0`, GitHub release → `publish-pypi.yml`, and the valence consumer update (Component 3). All blocked by the local-only constraint.

---

## Notes / Risks

- `griffe check` exits nonzero both on detected breakages **and** on tool/load errors (e.g. package not found). First cut treats any nonzero as "breaking"; if false positives appear, refine `main()` to distinguish (inspect stderr / use `--format`). Out of scope for this plan.
- The gate is forward-looking: if the latest tag already contains the `write_fort14` break, griffe vs HEAD shows no break — correct. The gate prevents the *next* silent break; the *existing* one is owned by Task 2.
