"""SemVer decisions for the release-integrity gate (0.x-aware)."""
from __future__ import annotations

import re
from typing import Tuple


def parse_version(v: str) -> Tuple[int, int, int]:
    """Return (major, minor, patch) from a version like '1.2.2' or 'v1.2.2'."""
    m = re.search(r"(\d+)\.(\d+)\.(\d+)", v)
    if not m:
        raise ValueError(f"unparseable version: {v!r}")
    return int(m.group(1)), int(m.group(2)), int(m.group(3))


def extract_version(pyproject_text: str) -> str:
    """Pull the project version string out of pyproject.toml text."""
    m = re.search(r'^version\s*=\s*"([^"]+)"', pyproject_text, re.M)
    if not m:
        raise ValueError("version not found in pyproject.toml")
    return m.group(1)


def decide(breaking: bool, old_version: str, new_version: str) -> Tuple[bool, str]:
    """Gate decision. Returns (ok, message).

    Post-1.0: a breaking change requires a MAJOR bump.
    Pre-1.0 (0.x): a breaking change requires at least a MINOR bump.
    """
    if not breaking:
        return True, "No breaking public-API change vs baseline."
    old = parse_version(old_version)
    new = parse_version(new_version)
    if old[0] >= 1:
        if new[0] > old[0]:
            return True, f"Breaking change; major bumped {old[0]}->{new[0]}. OK."
        return False, (
            f"Public API broke without a major bump (major still {old[0]}). "
            "Bump the major in pyproject.toml or revert the breaking change."
        )
    # 0.x semantics
    if new[0] > old[0] or new[1] > old[1]:
        return True, f"Breaking change on 0.x; minor bumped {old[1]}->{new[1]}. OK."
    return False, (
        f"Public API broke on 0.x without a minor bump (still 0.{old[1]}.x). "
        "Bump the minor in pyproject.toml or revert the breaking change."
    )
