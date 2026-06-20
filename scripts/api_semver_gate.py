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
    griffe = _run([sys.executable, "-m", "griffe", "check", "chilmesh", "-s", "src", "-a", tag])
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
