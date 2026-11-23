"""Run the README `.chil` example in a temporary directory."""
from __future__ import annotations

import re
from pathlib import Path

import chilmesh
from chilmesh import CHILmesh

README = Path(__file__).resolve().parent.parent / "README.md"


def _example_code() -> str:
    text = README.read_text(encoding="utf-8")
    m = re.search(r"<!-- chil-example -->\s*```python\n(.*?)```", text, re.S)
    assert m, "README .chil example marker not found"
    return m.group(1)


def test_readme_chil_example_runs(tmp_path, monkeypatch, capsys):
    chilmesh.examples.annulus().write_to_fort14(str(tmp_path / "ocean.14"))
    monkeypatch.chdir(tmp_path)
    exec(compile(_example_code(), "README.md", "exec"), {})
    out = capsys.readouterr().out
    assert "EPSG:4326" in out and "fort.14" in out
    assert out.strip().splitlines()[-1] == "True"
    assert (tmp_path / "ocean.chil").is_file()
    assert CHILmesh.load(tmp_path / "ocean.chil").n_elems > 0
