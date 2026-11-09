"""SC-005: no GPL competitor imports outside adapters/, no GPL text in the repo."""
from __future__ import annotations

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
ADAPTERS = REPO / "benchmarks" / "competitors" / "adapters"
SCAN_DIRS = ["src", "scripts", "tests", "benchmarks", "docs"]
SKIP_PARTS = {"__pycache__", ".git", "build", "target", ".venv", "venv", "node_modules"}
IMPORT_RE = re.compile(r"^\s*(?:import|from)\s+(oceanmesh|adcircpy)\b", re.M)
# Built by concatenation so this file does not match itself.
GPL_MARKERS = ["GNU GENERAL " + "PUBLIC LICENSE", "GNU AFFERO " + "GENERAL PUBLIC", "GNU LESSER " + "GENERAL PUBLIC"]
TEXT_SUFFIXES = {".py", ".md", ".txt", ".toml", ".cfg", ".m", ".rst", ".yml", ".yaml", ".in"}


def _files(suffixes):
    for d in SCAN_DIRS:
        root = REPO / d
        if not root.exists():
            continue
        for p in root.rglob("*"):
            if p.is_file() and p.suffix in suffixes and not (set(p.parts) & SKIP_PARTS):
                yield p


def test_no_gpl_imports_outside_adapters():
    bad = []
    for p in _files({".py"}):
        if ADAPTERS in p.parents or p == Path(__file__).resolve():
            continue
        if IMPORT_RE.search(p.read_text(encoding="utf-8", errors="ignore")):
            bad.append(str(p.relative_to(REPO)))
    assert not bad, f"GPL competitor imported outside adapters/: {bad}"


def test_no_gpl_license_text_in_repo():
    bad = []
    for p in _files(TEXT_SUFFIXES):
        text = p.read_text(encoding="utf-8", errors="ignore").upper()
        if any(m in text for m in GPL_MARKERS):
            bad.append(str(p.relative_to(REPO)))
    assert not bad, f"GPL text found: {bad}"


def test_adapters_are_standalone_scripts_not_package_imports():
    """Harness code must never import adapters; they run as subprocesses only."""
    bad = []
    for p in _files({".py"}):
        if ADAPTERS in p.parents or p == Path(__file__).resolve():
            continue
        if re.search(r"^\s*(?:import|from)\s+[\w.]*adapters\b", p.read_text(encoding="utf-8", errors="ignore"), re.M):
            bad.append(str(p.relative_to(REPO)))
    assert not bad, f"adapters imported in-process: {bad}"


def test_benchmarks_pruned_from_sdist():
    assert "prune benchmarks" in (REPO / "MANIFEST.in").read_text().splitlines()
