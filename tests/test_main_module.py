"""Tests for chilmesh CLI entry point (python -m chilmesh).

Verifies that the module can be invoked as a script and that argument
parsing works correctly.
"""
from __future__ import annotations

import subprocess
import sys

import pytest


class TestMainModuleInvocation:
    """Test chilmesh as a runnable module."""

    def test_help_succeeds(self):
        """python -m chilmesh --help returns exit code 0 with usage text."""
        result = subprocess.run(
            [sys.executable, "-m", "chilmesh", "--help"],
            capture_output=True,
            text=True,
            env={"MPLBACKEND": "Agg"},
        )

        assert result.returncode == 0
        assert "usage" in result.stdout.lower()
        assert "chilmesh" in result.stdout.lower()

    def test_no_args_fails(self):
        """python -m chilmesh with no arguments returns non-zero exit code."""
        result = subprocess.run(
            [sys.executable, "-m", "chilmesh"],
            capture_output=True,
            text=True,
            env={"MPLBACKEND": "Agg"},
        )

        # argparse requires a subcommand; should exit with code 2
        assert result.returncode == 2
        # stderr should contain the error message
        assert "required: COMMAND" in result.stderr or "required" in result.stderr

    def test_version_flag_succeeds(self):
        """python -m chilmesh --version returns exit code 0."""
        result = subprocess.run(
            [sys.executable, "-m", "chilmesh", "--version"],
            capture_output=True,
            text=True,
            env={"MPLBACKEND": "Agg"},
        )

        assert result.returncode == 0
        # Should print version (e.g., "chilmesh 1.0.0" or similar)
        assert "chilmesh" in result.stdout.lower()

    def test_unknown_command_fails(self):
        """python -m chilmesh nonexistent returns non-zero exit code."""
        result = subprocess.run(
            [sys.executable, "-m", "chilmesh", "nonexistent"],
            capture_output=True,
            text=True,
            env={"MPLBACKEND": "Agg"},
        )

        # Unknown subcommand should fail
        assert result.returncode != 0

    def test_info_without_input_fails(self):
        """chilmesh info (missing required input arg) fails."""
        result = subprocess.run(
            [sys.executable, "-m", "chilmesh", "info"],
            capture_output=True,
            text=True,
            env={"MPLBACKEND": "Agg"},
        )

        # Missing required positional argument
        assert result.returncode != 0
        assert "required" in result.stderr.lower() or "required" in result.stdout.lower()

    def test_convert_without_args_fails(self):
        """chilmesh convert (missing required args) fails."""
        result = subprocess.run(
            [sys.executable, "-m", "chilmesh", "convert"],
            capture_output=True,
            text=True,
            env={"MPLBACKEND": "Agg"},
        )

        # Missing required arguments
        assert result.returncode != 0

    def test_smooth_without_output_fails(self):
        """chilmesh smooth (missing required -o) fails."""
        result = subprocess.run(
            [sys.executable, "-m", "chilmesh", "smooth", "dummy.14"],
            capture_output=True,
            text=True,
            env={"MPLBACKEND": "Agg"},
        )

        # Missing required -o / --output argument
        assert result.returncode != 0

    def test_plot_without_output_fails(self):
        """chilmesh plot (missing required -o) fails."""
        result = subprocess.run(
            [sys.executable, "-m", "chilmesh", "plot", "dummy.14"],
            capture_output=True,
            text=True,
            env={"MPLBACKEND": "Agg"},
        )

        # Missing required -o / --output argument
        assert result.returncode != 0

    def test_summary_help_succeeds(self):
        """chilmesh summary --help returns exit code 0."""
        result = subprocess.run(
            [sys.executable, "-m", "chilmesh", "summary", "--help"],
            capture_output=True,
            text=True,
            env={"MPLBACKEND": "Agg"},
        )

        assert result.returncode == 0
        assert "summary" in result.stdout.lower()

    def test_info_help_succeeds(self):
        """chilmesh info --help returns exit code 0."""
        result = subprocess.run(
            [sys.executable, "-m", "chilmesh", "info", "--help"],
            capture_output=True,
            text=True,
            env={"MPLBACKEND": "Agg"},
        )

        assert result.returncode == 0
        assert "info" in result.stdout.lower()

    def test_convert_help_succeeds(self):
        """chilmesh convert --help returns exit code 0."""
        result = subprocess.run(
            [sys.executable, "-m", "chilmesh", "convert", "--help"],
            capture_output=True,
            text=True,
            env={"MPLBACKEND": "Agg"},
        )

        assert result.returncode == 0
        assert "convert" in result.stdout.lower()

    def test_smooth_help_succeeds(self):
        """chilmesh smooth --help returns exit code 0."""
        result = subprocess.run(
            [sys.executable, "-m", "chilmesh", "smooth", "--help"],
            capture_output=True,
            text=True,
            env={"MPLBACKEND": "Agg"},
        )

        assert result.returncode == 0
        assert "smooth" in result.stdout.lower()

    def test_plot_help_succeeds(self):
        """chilmesh plot --help returns exit code 0."""
        result = subprocess.run(
            [sys.executable, "-m", "chilmesh", "plot", "--help"],
            capture_output=True,
            text=True,
            env={"MPLBACKEND": "Agg"},
        )

        assert result.returncode == 0
        assert "plot" in result.stdout.lower()
