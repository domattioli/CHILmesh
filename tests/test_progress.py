"""Tests for the private progress line (spec 008 FR-8..FR-14)."""
import sys
import time

import pytest

from chilmesh import _progress
from chilmesh import CHILmesh
from chilmesh.CHILmesh import write_fort14


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    monkeypatch.delenv("CHILMESH_PROGRESS", raising=False)
    # Default: pretend tqdm is absent so the stdlib renderer is tested.
    monkeypatch.setattr(_progress, "_load_tqdm", lambda: None)


def _run(n=3, sleep=0.0, **kw):
    with _progress.track("Peeling layers", "layer", past="Peeled", **kw) as t:
        for _ in range(n):
            if sleep:
                time.sleep(sleep)
            t.tick()


def test_silent_when_not_tty_in_auto(capsys, monkeypatch):
    monkeypatch.setattr(_progress, "_is_tty", lambda: False)
    monkeypatch.setattr(_progress, "_DELAY", 0.0)
    _run(sleep=0.05)
    out = capsys.readouterr()
    assert out.err == "" and out.out == ""


def test_silent_under_one_second_on_tty(capsys, monkeypatch):
    monkeypatch.setattr(_progress, "_is_tty", lambda: True)
    _run()
    out = capsys.readouterr()
    assert out.err == "" and out.out == ""


def test_auto_tty_draws_after_delay(capsys, monkeypatch):
    monkeypatch.setattr(_progress, "_is_tty", lambda: True)
    monkeypatch.setattr(_progress, "_DELAY", 0.05)
    monkeypatch.setattr(_progress, "_FRAME_INTERVAL", 0.02)
    _run(n=5, sleep=0.1)
    out = capsys.readouterr()
    assert out.out == ""
    assert "Peeling layers…" in out.err
    assert "Chillin'" in out.err


def test_env_one_forces_output(capsys, monkeypatch):
    monkeypatch.setenv("CHILMESH_PROGRESS", "1")
    monkeypatch.setattr(_progress, "_is_tty", lambda: False)
    _run(n=4)
    out = capsys.readouterr()
    assert out.out == ""
    assert "✓ Chillin'. Peeled 4 layers in " in out.err


def test_env_zero_suppresses(capsys, monkeypatch):
    monkeypatch.setenv("CHILMESH_PROGRESS", "0")
    monkeypatch.setattr(_progress, "_is_tty", lambda: True)
    monkeypatch.setattr(_progress, "_DELAY", 0.0)
    _run(sleep=0.05)
    out = capsys.readouterr()
    assert out.err == "" and out.out == ""


def test_private_kwarg_overrides_env(capsys, monkeypatch):
    monkeypatch.setenv("CHILMESH_PROGRESS", "0")
    _run(n=2, _progress=True)
    assert "Chillin'" in capsys.readouterr().err
    monkeypatch.setenv("CHILMESH_PROGRESS", "1")
    _run(n=2, _progress=False)
    assert capsys.readouterr().err == ""


def test_percent_only_when_total_known(capsys, monkeypatch):
    monkeypatch.setenv("CHILMESH_PROGRESS", "1")
    with _progress.track("Smoothing", "iteration", 10, past="Smoothed") as t:
        t.tick(5)
        time.sleep(0.3)
    err = capsys.readouterr().err
    assert "50%" in err
    _run(n=2)
    assert "%" not in capsys.readouterr().err


def test_no_done_line_on_exception(capsys, monkeypatch):
    monkeypatch.setenv("CHILMESH_PROGRESS", "1")
    with pytest.raises(RuntimeError):
        with _progress.track("Peeling layers", "layer"):
            raise RuntimeError("boom")
    assert "Chillin'" not in capsys.readouterr().err


def test_tqdm_path(capsys, monkeypatch):
    calls = {}

    class FakeBar:
        def __init__(self, **kw):
            calls["init"] = kw
            calls["updates"] = 0

        def update(self, n):
            calls["updates"] += n

        def close(self):
            calls["closed"] = True

    monkeypatch.setattr(_progress, "_load_tqdm", lambda: FakeBar)
    monkeypatch.setenv("CHILMESH_PROGRESS", "1")
    _run(n=3)
    assert calls["init"]["desc"] == "Peeling layers"
    assert calls["updates"] == 3
    assert calls["closed"] is True
    out = capsys.readouterr()
    assert out.out == ""
    assert "Chillin'" in out.err


def test_no_tqdm_path_when_import_fails(capsys, monkeypatch):
    monkeypatch.undo()  # drop the autouse stub so the real import runs
    monkeypatch.setitem(sys.modules, "tqdm", None)  # import tqdm -> ImportError
    assert _progress._load_tqdm() is None
    monkeypatch.setenv("CHILMESH_PROGRESS", "1")
    _run(n=2)
    err = capsys.readouterr().err
    assert "Peeled 2 layers" in err


def test_wired_ops_silent_stdout_and_forced_stderr(capsys, monkeypatch, tmp_path):
    import numpy as np

    pts = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0.5, 0.5, 0]], float)
    conn = np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]])
    monkeypatch.setenv("CHILMESH_PROGRESS", "1")
    mesh = CHILmesh(conn, pts)
    path = str(tmp_path / "m.14")
    mesh.write_to_fort14(path)
    CHILmesh.read_from_fort14(path)
    mesh.elem_quality()
    mesh.angle_based_smoother(n_iter=2)
    out = capsys.readouterr()
    assert out.out == ""
    for verb in ("Peeling layers", "Writing fort.14", "Reading fort.14",
                 "Analyzing quality", "Smoothing"):
        assert verb in out.err or verb.split()[0] in out.err
    assert out.err.count("Chillin'") >= 5
    # explicit off wins over env
    mesh.write_to_fort14(path, _progress=False)
    assert capsys.readouterr().err == ""
    write_fort14(mesh, path, _progress=False)
    assert capsys.readouterr().err == ""
