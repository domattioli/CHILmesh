"""End-to-end tests for the ``chilmesh`` CLI.

Spawns a subprocess so the installed entry point — or ``python -m chilmesh``
when the script isn't on PATH — is exercised exactly the way a shell user
would invoke it. Tests run against the bundled ``annulus`` fixture to keep
runtime modest.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from chilmesh import examples


def _run(args, **kw):
    """Invoke the CLI via ``python -m chilmesh`` for path-independence."""
    cmd = [sys.executable, "-m", "chilmesh", *args]
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


@pytest.fixture
def annulus_path():
    return examples.fixture_path("annulus_200pts.fort.14")


# ----------------------------------------------------------------------
# Top-level help / discovery
# ----------------------------------------------------------------------


def test_help_lists_subcommands():
    result = _run(["--help"])
    assert result.returncode == 0
    for sub in ("info", "convert", "smooth", "plot"):
        assert sub in result.stdout, f"--help missing subcommand '{sub}'"


def test_version_prints():
    result = _run(["--version"])
    assert result.returncode == 0
    assert "chilmesh" in result.stdout.lower()


@pytest.mark.parametrize("sub", ["info", "convert", "smooth", "plot"])
def test_each_subcommand_has_help(sub):
    result = _run([sub, "--help"])
    assert result.returncode == 0, result.stderr
    # argparse renders an Example: block in each subcommand description.
    assert "Example" in result.stdout


def test_missing_command_errors():
    result = _run([])
    assert result.returncode != 0


# ----------------------------------------------------------------------
# info
# ----------------------------------------------------------------------


def test_info_happy_path(annulus_path):
    result = _run(["info", str(annulus_path)])
    assert result.returncode == 0, result.stderr
    out = result.stdout
    assert "Vertices:" in out
    assert "Elements:" in out
    assert "Edges:" in out
    assert "Layers:" in out
    assert "Quality (skew):" in out


def test_info_no_layers_skips_layer_count(annulus_path):
    result = _run(["info", str(annulus_path), "--no-layers", "--no-quality"])
    assert result.returncode == 0, result.stderr
    assert "(skipped, --no-layers)" in result.stdout
    assert "Quality (skew)" not in result.stdout


def test_info_missing_file_errors(tmp_path):
    bogus = tmp_path / "does_not_exist.fort.14"
    result = _run(["info", str(bogus)])
    assert result.returncode != 0
    # Either FileNotFoundError or open-failure; stderr should mention it.
    assert result.stderr or result.stdout


def test_info_unknown_format_errors(tmp_path):
    bogus = tmp_path / "mesh.xyz"
    bogus.write_text("not a real mesh")
    result = _run(["info", str(bogus)])
    assert result.returncode != 0
    assert "format" in (result.stderr + result.stdout).lower()


# ----------------------------------------------------------------------
# convert
# ----------------------------------------------------------------------


def test_convert_fort14_roundtrip(annulus_path, tmp_path):
    out_path = tmp_path / "out.fort.14"
    result = _run(["convert", str(annulus_path), str(out_path)])
    assert result.returncode == 0, result.stderr
    assert out_path.exists()
    # File should be non-empty and reload cleanly via the same CLI.
    assert out_path.stat().st_size > 0
    info = _run(["info", str(out_path), "--no-quality", "--no-layers"])
    assert info.returncode == 0, info.stderr


def test_convert_unknown_output_format_errors(annulus_path, tmp_path):
    out_path = tmp_path / "out.bogus"
    result = _run(["convert", str(annulus_path), str(out_path)])
    assert result.returncode != 0
    assert not out_path.exists()


# ----------------------------------------------------------------------
# smooth
# ----------------------------------------------------------------------


def test_smooth_angle_based_writes_output(annulus_path, tmp_path):
    out_path = tmp_path / "smoothed.fort.14"
    result = _run([
        "smooth", str(annulus_path),
        "-o", str(out_path),
        "--method", "angle-based",
        "--iter", "1",
    ])
    assert result.returncode == 0, result.stderr
    assert out_path.exists()
    assert "median quality" in result.stdout


def test_smooth_rejects_unknown_method(annulus_path, tmp_path):
    out_path = tmp_path / "smoothed.fort.14"
    result = _run([
        "smooth", str(annulus_path),
        "-o", str(out_path),
        "--method", "not-a-method",
    ])
    assert result.returncode != 0
    assert "not-a-method" in (result.stderr + result.stdout)


# ----------------------------------------------------------------------
# plot
# ----------------------------------------------------------------------


def test_plot_writes_png(annulus_path, tmp_path):
    out_path = tmp_path / "mesh.png"
    result = _run(["plot", str(annulus_path), "-o", str(out_path)])
    assert result.returncode == 0, result.stderr
    assert out_path.exists()
    assert out_path.stat().st_size > 0


def test_plot_quality_overlay(annulus_path, tmp_path):
    out_path = tmp_path / "quality.png"
    result = _run([
        "plot", str(annulus_path), "-o", str(out_path), "--quality",
    ])
    assert result.returncode == 0, result.stderr
    assert out_path.exists()


def test_plot_layers_overlay(annulus_path, tmp_path):
    out_path = tmp_path / "layers.png"
    result = _run([
        "plot", str(annulus_path), "-o", str(out_path), "--layers",
    ])
    assert result.returncode == 0, result.stderr
    assert out_path.exists()


def test_plot_requires_output(annulus_path):
    result = _run(["plot", str(annulus_path)])
    assert result.returncode != 0
