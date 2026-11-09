"""SC-006: no 'Python 3.10' left in README; other docs listed as known gaps."""
from __future__ import annotations

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
PY310 = re.compile(r"python\s*3\.10|py3\.10|cp310|>=\s*3\.10", re.I)


def test_readme_has_no_python_310_mention():
    text = (REPO / "README.md").read_text(encoding="utf-8")
    hits = [m.group(0) for m in PY310.finditer(text)]
    assert not hits, f"README still mentions Python 3.10: {hits}"
    assert "Python 3.11" in text


def test_readme_requires_311_consistently():
    text = (REPO / "README.md").read_text(encoding="utf-8")
    assert "python-3.11%2B" in text  # badge
    assert re.search(r"Python 3\.11 or newer", text)


def test_pyproject_requires_311():
    assert 'requires-python = ">=3.11"' in (REPO / "pyproject.toml").read_text()


def test_contributing_has_no_python_310_mention():
    text = (REPO / "CONTRIBUTING.md").read_text(encoding="utf-8")
    assert not PY310.search(text) and "3.10+" not in text


def test_workflows_do_not_target_310():
    bad = []
    for wf in (REPO / ".github" / "workflows").glob("*.yml"):
        if re.search(r"""python-version:.*["']3\.10["']""", wf.read_text(encoding="utf-8")):
            bad.append(wf.name)
    assert not bad
