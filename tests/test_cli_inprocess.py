"""In-process tests for the ``chilmesh`` CLI.

These tests call ``main(argv)`` and the subcommand functions directly
(no subprocess), providing honest coverage of error-handling branches and
CLI logic that subprocess-based tests cannot reach. Complements test_cli.py
by exercising the in-process code paths and capturing stdout/stderr via
pytest's capsys fixture.

Key difference from test_cli.py:
- test_cli.py spawns subprocesses (child process invisible to coverage)
- test_cli_inprocess.py calls main() and cmd_* functions directly,
  giving full coverage of cli.py error paths and exit codes.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from chilmesh import examples
from chilmesh.cli import (
    _detect_format,
    cmd_convert,
    cmd_info,
    cmd_plot,
    cmd_smooth,
    cmd_summary,
    main,
)


# ========================================================================
# Fixtures
# ========================================================================


@pytest.fixture
def quad_2x2_path():
    """Path to the tiny 2×2 quad fixture (9 nodes, 4 quads, fast)."""
    return examples.fixture_path("quad_2x2.fort.14")


@pytest.fixture
def annulus_path():
    """Path to the small annulus fixture (~580 elements, small)."""
    return examples.fixture_path("annulus_200pts.fort.14")


# ========================================================================
# Top-level: --version and empty argv
# ========================================================================


def test_main_version_raises_systemexit(capsys):
    """--version prints version string and raises SystemExit(0)."""
    with pytest.raises(SystemExit) as exc_info:
        main(["--version"])
    assert exc_info.value.code == 0
    captured = capsys.readouterr()
    assert "chilmesh" in captured.out.lower()


def test_main_empty_argv_raises_systemexit(capsys):
    """Empty argv triggers argparse error, raises SystemExit(2)."""
    with pytest.raises(SystemExit) as exc_info:
        main([])
    # argparse raises SystemExit(2) for missing required subcommand
    assert exc_info.value.code == 2


# ========================================================================
# _detect_format private function
# ========================================================================


def test_detect_format_fort14_variants():
    """_detect_format recognizes .14, .fort.14, .grd as fort14."""
    assert _detect_format(Path("mesh.14")) == "fort14"
    assert _detect_format(Path("mesh.fort.14")) == "fort14"
    assert _detect_format(Path("mesh.grd")) == "fort14"
    # Case-insensitive
    assert _detect_format(Path("MESH.14")) == "fort14"
    assert _detect_format(Path("Mesh.FORT.14")) == "fort14"


def test_detect_format_2dm():
    """_detect_format recognizes .2dm."""
    assert _detect_format(Path("mesh.2dm")) == "2dm"
    assert _detect_format(Path("MESH.2DM")) == "2dm"


def test_detect_format_unknown_suffix_raises_systemexit():
    """_detect_format raises SystemExit for unknown suffixes."""
    with pytest.raises(SystemExit) as exc_info:
        _detect_format(Path("mesh.xyz"))
    assert "format" in str(exc_info.value).lower()
    assert "xyz" in str(exc_info.value)


# ========================================================================
# cmd_info
# ========================================================================


def test_cmd_info_happy_path(quad_2x2_path, capsys):
    """cmd_info on quad_2x2: returns 0, prints mesh stats."""
    import argparse
    args = argparse.Namespace(
        input=str(quad_2x2_path),
        no_layers=False,
        no_quality=False,
    )
    result = cmd_info(args)
    assert result == 0
    captured = capsys.readouterr()
    out = captured.out
    assert "Vertices:" in out
    assert "Elements:" in out
    assert "Edges:" in out
    assert "Layers:" in out
    assert "Quality (skew):" in out
    # Spot-check some values
    assert "min:" in out
    assert "median:" in out


def test_cmd_info_no_layers_no_quality(quad_2x2_path, capsys):
    """cmd_info with --no-layers --no-quality skips those sections."""
    import argparse
    args = argparse.Namespace(
        input=str(quad_2x2_path),
        no_layers=True,
        no_quality=True,
    )
    result = cmd_info(args)
    assert result == 0
    captured = capsys.readouterr()
    out = captured.out
    assert "(skipped, --no-layers)" in out
    assert "Quality (skew):" not in out


def test_main_info_missing_file_returns_2(tmp_path, capsys):
    """main(['info', nonexistent_path]) catches FileNotFoundError, returns 2."""
    bogus = tmp_path / "does_not_exist.fort.14"
    result = main(["info", str(bogus)])
    assert result == 2
    captured = capsys.readouterr()
    # Error message printed to stderr
    error_out = captured.err + captured.out
    assert "error:" in error_out


# ========================================================================
# cmd_summary
# ========================================================================


def test_cmd_summary_happy_path(quad_2x2_path, capsys):
    """cmd_summary on quad_2x2: returns 0, prints lazy metadata."""
    import argparse
    args = argparse.Namespace(
        input=str(quad_2x2_path),
        deep=False,
    )
    result = cmd_summary(args)
    assert result == 0
    captured = capsys.readouterr()
    out = captured.out
    assert "Summary:" in out
    assert "Nodes:" in out
    assert "Elements:" in out


def test_cmd_summary_deep_flag(quad_2x2_path, capsys):
    """cmd_summary with --deep includes Type and Bbox."""
    import argparse
    args = argparse.Namespace(
        input=str(quad_2x2_path),
        deep=True,
    )
    result = cmd_summary(args)
    assert result == 0
    captured = capsys.readouterr()
    out = captured.out
    assert "Type:" in out
    assert "Bbox:" in out


def test_main_summary_missing_file_returns_1(tmp_path, capsys):
    """main(['summary', nonexistent_path]) raises generic Exception, returns 1.

    Note: summary() calls stat() which raises a different exception than
    FileNotFoundError, so main's generic Exception handler (not the
    FileNotFoundError handler) catches it and returns 1.
    """
    bogus = tmp_path / "does_not_exist.fort.14"
    result = main(["summary", str(bogus)])
    assert result == 1
    captured = capsys.readouterr()
    error_out = captured.err + captured.out
    assert "error:" in error_out


# ========================================================================
# cmd_convert
# ========================================================================


def test_cmd_convert_happy_path(quad_2x2_path, tmp_path, capsys):
    """cmd_convert quad_2x2 → tmp out.fort.14: returns 0, writes file."""
    import argparse
    out_path = tmp_path / "out.fort.14"
    args = argparse.Namespace(
        input=str(quad_2x2_path),
        output=str(out_path),
    )
    result = cmd_convert(args)
    assert result == 0
    assert out_path.exists()
    assert out_path.stat().st_size > 0
    captured = capsys.readouterr()
    assert "Wrote" in captured.out


def test_main_convert_unsupported_output_format_raises_systemexit(
    quad_2x2_path, tmp_path
):
    """main(['convert', input, output.2dm]) raises SystemExit (2dm unsupported)."""
    out_path = tmp_path / "out.2dm"
    with pytest.raises(SystemExit) as exc_info:
        main(["convert", str(quad_2x2_path), str(out_path)])
    # _write_mesh raises SystemExit; main re-raises it
    assert "2dm" in str(exc_info.value).lower() or str(exc_info.value).startswith(
        "error:"
    )


# ========================================================================
# cmd_smooth
# ========================================================================


def test_cmd_smooth_angle_based_writes_output(quad_2x2_path, tmp_path, capsys):
    """cmd_smooth with angle-based: returns 0, writes smoothed mesh."""
    import argparse
    out_path = tmp_path / "smoothed.fort.14"
    args = argparse.Namespace(
        input=str(quad_2x2_path),
        output=str(out_path),
        method="angle-based",
        iter=1,
    )
    result = cmd_smooth(args)
    assert result == 0
    assert out_path.exists()
    assert out_path.stat().st_size > 0
    captured = capsys.readouterr()
    assert "Smoothed" in captured.out
    assert "median quality" in captured.out


def test_cmd_smooth_laplacian_maps_to_fem(quad_2x2_path, tmp_path, capsys):
    """cmd_smooth with laplacian: returns 0 (maps internally to fem)."""
    import argparse
    out_path = tmp_path / "smoothed.fort.14"
    args = argparse.Namespace(
        input=str(quad_2x2_path),
        output=str(out_path),
        method="laplacian",
        iter=1,
    )
    result = cmd_smooth(args)
    assert result == 0
    assert out_path.exists()
    captured = capsys.readouterr()
    assert "Smoothed" in captured.out


def test_cmd_smooth_fem_explicit(quad_2x2_path, tmp_path, capsys):
    """cmd_smooth with fem: returns 0."""
    import argparse
    out_path = tmp_path / "smoothed.fort.14"
    args = argparse.Namespace(
        input=str(quad_2x2_path),
        output=str(out_path),
        method="fem",
        iter=1,
    )
    result = cmd_smooth(args)
    assert result == 0
    assert out_path.exists()


# ========================================================================
# cmd_plot
# ========================================================================


def test_cmd_plot_default_writes_png(quad_2x2_path, tmp_path, capsys):
    """cmd_plot default (no flags) writes PNG, returns 0."""
    import argparse
    out_path = tmp_path / "mesh.png"
    args = argparse.Namespace(
        input=str(quad_2x2_path),
        output=str(out_path),
        layers=False,
        quality=False,
        dpi=150,
    )
    result = cmd_plot(args)
    assert result == 0
    assert out_path.exists()
    assert out_path.stat().st_size > 0
    captured = capsys.readouterr()
    assert "Wrote" in captured.out


def test_cmd_plot_quality_overlay(quad_2x2_path, tmp_path, capsys):
    """cmd_plot with --quality: returns 0, writes quality-colored PNG."""
    import argparse
    out_path = tmp_path / "quality.png"
    args = argparse.Namespace(
        input=str(quad_2x2_path),
        output=str(out_path),
        layers=False,
        quality=True,
        dpi=150,
    )
    result = cmd_plot(args)
    assert result == 0
    assert out_path.exists()
    assert out_path.stat().st_size > 0


def test_cmd_plot_layers_overlay(quad_2x2_path, tmp_path, capsys):
    """cmd_plot with --layers: returns 0, writes layer-colored PNG."""
    import argparse
    out_path = tmp_path / "layers.png"
    args = argparse.Namespace(
        input=str(quad_2x2_path),
        output=str(out_path),
        layers=True,
        quality=False,
        dpi=150,
    )
    result = cmd_plot(args)
    assert result == 0
    assert out_path.exists()
    assert out_path.stat().st_size > 0


def test_cmd_plot_custom_dpi(quad_2x2_path, tmp_path, capsys):
    """cmd_plot with custom --dpi: returns 0."""
    import argparse
    out_path = tmp_path / "high_dpi.png"
    args = argparse.Namespace(
        input=str(quad_2x2_path),
        output=str(out_path),
        layers=False,
        quality=False,
        dpi=300,
    )
    result = cmd_plot(args)
    assert result == 0
    assert out_path.exists()


# ========================================================================
# Full main() integration tests
# ========================================================================


def test_main_info_flow(quad_2x2_path, capsys):
    """main(['info', path]) dispatches to cmd_info, returns 0."""
    result = main(["info", str(quad_2x2_path)])
    assert result == 0
    captured = capsys.readouterr()
    assert "Vertices:" in captured.out


def test_main_summary_flow(quad_2x2_path, capsys):
    """main(['summary', path]) dispatches to cmd_summary, returns 0."""
    result = main(["summary", str(quad_2x2_path)])
    assert result == 0
    captured = capsys.readouterr()
    assert "Summary:" in captured.out


def test_main_convert_flow(quad_2x2_path, tmp_path, capsys):
    """main(['convert', in, out]) dispatches, returns 0."""
    out_path = tmp_path / "converted.fort.14"
    result = main(["convert", str(quad_2x2_path), str(out_path)])
    assert result == 0
    assert out_path.exists()


def test_main_smooth_flow(quad_2x2_path, tmp_path, capsys):
    """main(['smooth', in, '-o', out, '--method', 'angle-based']) returns 0."""
    out_path = tmp_path / "smoothed.fort.14"
    result = main(
        [
            "smooth",
            str(quad_2x2_path),
            "-o",
            str(out_path),
            "--method",
            "angle-based",
            "--iter",
            "1",
        ]
    )
    assert result == 0
    assert out_path.exists()


def test_main_plot_flow(quad_2x2_path, tmp_path, capsys):
    """main(['plot', in, '-o', out]) dispatches to cmd_plot, returns 0."""
    out_path = tmp_path / "plot.png"
    result = main(["plot", str(quad_2x2_path), "-o", str(out_path)])
    assert result == 0
    assert out_path.exists()


def test_main_plot_with_quality_flag(quad_2x2_path, tmp_path):
    """main(['plot', in, '-o', out, '--quality']) returns 0."""
    out_path = tmp_path / "quality_plot.png"
    result = main(["plot", str(quad_2x2_path), "-o", str(out_path), "--quality"])
    assert result == 0
    assert out_path.exists()


def test_main_plot_with_layers_flag(quad_2x2_path, tmp_path):
    """main(['plot', in, '-o', out, '--layers']) returns 0."""
    out_path = tmp_path / "layers_plot.png"
    result = main(["plot", str(quad_2x2_path), "-o", str(out_path), "--layers"])
    assert result == 0
    assert out_path.exists()


# ========================================================================
# Annulus tests (slightly larger fixture for quality/layer tests)
# ========================================================================


def test_main_info_annulus_with_quality(annulus_path, capsys):
    """main(['info', annulus]) shows detailed quality stats."""
    result = main(["info", str(annulus_path)])
    assert result == 0
    captured = capsys.readouterr()
    out = captured.out
    # Annulus has more structure, should show all quality fields
    assert "min:" in out
    assert "max:" in out
    assert "std:" in out


def test_main_plot_annulus_layers(annulus_path, tmp_path):
    """main(['plot', annulus, '-o', out, '--layers']) works on larger mesh."""
    out_path = tmp_path / "annulus_layers.png"
    result = main(["plot", str(annulus_path), "-o", str(out_path), "--layers"])
    assert result == 0
    assert out_path.exists()
    assert out_path.stat().st_size > 0


# ========================================================================
# Error handling edge cases
# ========================================================================


def test_main_unknown_format_input_returns_systemexit(tmp_path):
    """main(['info', file.xyz]) raises SystemExit (unknown format)."""
    bogus = tmp_path / "mesh.xyz"
    bogus.write_text("not a real mesh")
    with pytest.raises(SystemExit):
        main(["info", str(bogus)])


def test_main_convert_unsupported_input_format_raises_systemexit(tmp_path):
    """main(['convert', file.xyz, out.fort.14]) raises SystemExit."""
    bogus = tmp_path / "mesh.xyz"
    bogus.write_text("not a real mesh")
    out_path = tmp_path / "out.fort.14"
    with pytest.raises(SystemExit):
        main(["convert", str(bogus), str(out_path)])


def test_cmd_info_nonexistent_path_raises_fileerror(tmp_path):
    """Calling cmd_info directly with nonexistent path raises FileNotFoundError."""
    import argparse
    bogus = tmp_path / "does_not_exist.fort.14"
    args = argparse.Namespace(
        input=str(bogus),
        no_layers=False,
        no_quality=False,
    )
    # cmd_info does not catch FileNotFoundError; it propagates up to main
    with pytest.raises(FileNotFoundError):
        cmd_info(args)


# ========================================================================
# Argument validation (via main, testing argparse integration)
# ========================================================================


def test_main_smooth_missing_output_flag_raises_systemexit(quad_2x2_path):
    """smooth without -o/--output raises SystemExit (required arg missing)."""
    with pytest.raises(SystemExit):
        main(["smooth", str(quad_2x2_path)])


def test_main_plot_missing_output_flag_raises_systemexit(quad_2x2_path):
    """plot without -o/--output raises SystemExit."""
    with pytest.raises(SystemExit):
        main(["plot", str(quad_2x2_path)])


def test_main_info_no_input_raises_systemexit():
    """info with no input path raises SystemExit."""
    with pytest.raises(SystemExit):
        main(["info"])
