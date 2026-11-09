"""CLI: real runs are opt-in behind CHILMESH_RUN_BENCH=1; README runbook covers the essentials."""
from __future__ import annotations

from pathlib import Path

import pytest

from benchmarks.competitors import cli

REPO = Path(__file__).resolve().parents[2]


def test_real_run_refused_without_opt_in(monkeypatch, capsys):
    monkeypatch.delenv("CHILMESH_RUN_BENCH", raising=False)
    rc = cli.main(["w2", "--mesh", "STOFS", "--out", "/tmp/never-written"])
    assert rc == 2
    assert "CHILMESH_RUN_BENCH=1" in capsys.readouterr().err


def test_w1_also_gated(monkeypatch, capsys):
    monkeypatch.delenv("CHILMESH_RUN_BENCH", raising=False)
    assert cli.main(["w1", "--mesh", "Shinnecock", "--out", "/tmp/never-written"]) == 2


def test_unknown_command_exits_nonzero():
    with pytest.raises(SystemExit):
        cli.main(["w3"])


def test_report_subcommand_works_without_opt_in(tmp_path, monkeypatch, capsys):
    monkeypatch.delenv("CHILMESH_RUN_BENCH", raising=False)
    (tmp_path / "records").mkdir()
    rc = cli.main(["report", "--records", str(tmp_path / "records"), "--workflow", "w2"])
    assert rc == 0
    assert "no records" in capsys.readouterr().out.lower()


def test_no_w3_command():
    assert "w3" not in cli.COMMANDS


def test_runbook_mentions_required_steps():
    text = (REPO / "benchmarks" / "competitors" / "README.md").read_text(encoding="utf-8")
    for needle in ("VALENCE_TOKEN", "CHILMESH_RUN_BENCH=1", "cold", "drop_caches", "purge",
                   "requirements-", "python -m benchmarks.competitors w2", "python -m benchmarks.competitors w1",
                   "rights.json", "second person"):
        assert needle in text, needle
    # README for the harness must be python 3.11+, not 3.10
    assert "3.10" not in text


def test_env_requirements_exist_and_are_not_in_pyproject():
    envs = list((REPO / "benchmarks" / "competitors" / "envs").glob("requirements-*.txt"))
    names = {p.stem for p in envs}
    assert {"requirements-adcircpy", "requirements-thalassa", "requirements-oceanmesh",
            "requirements-xarray"} <= names
    pyproject = (REPO / "pyproject.toml").read_text().lower()
    for dep in ("oceanmesh", "adcircpy", "thalassa"):
        assert dep not in pyproject
