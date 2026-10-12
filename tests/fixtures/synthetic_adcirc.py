"""Synthetic ADCIRC-format writers for tests (skeleton; extended in Stage B2/B3)."""
from __future__ import annotations

from pathlib import Path


def write_text(path: Path, lines: list[str]) -> Path:
    """Write ``lines`` to ``path`` with LF endings and return the path."""
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
