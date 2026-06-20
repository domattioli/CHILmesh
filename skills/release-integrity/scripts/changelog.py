"""CHANGELOG.md checks for the release-integrity gate (Keep a Changelog)."""
from __future__ import annotations

import re


def has_dated_section(text: str, version: str) -> bool:
    """True if a `## [version] — YYYY-MM-DD` heading exists (em dash or hyphen)."""
    pattern = rf'^##\s*\[{re.escape(version)}\]\s*[—-]\s*\d{{4}}-\d{{2}}-\d{{2}}'
    return re.search(pattern, text, re.M) is not None


def has_unreleased_bullet(text: str) -> bool:
    """True if the `## [Unreleased]` section contains at least one `- ` bullet."""
    m = re.search(r'^##\s*\[Unreleased\]\s*$', text, re.M)
    if not m:
        return False
    rest = text[m.end():]
    nxt = re.search(r'^##\s', rest, re.M)
    section = rest[: nxt.start()] if nxt else rest
    return re.search(r'^\s*-\s+\S', section, re.M) is not None
