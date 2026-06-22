from __future__ import annotations

import pathlib
import sys

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
import semver_core as s  # noqa: E402


def test_parse_version():
    assert s.parse_version("1.2.2") == (1, 2, 2)
    assert s.parse_version("v2.0.0") == (2, 0, 0)


def test_parse_version_rejects_garbage():
    with pytest.raises(ValueError):
        s.parse_version("not-a-version")


def test_extract_version():
    assert s.extract_version('name = "x"\nversion = "1.2.2"\n') == "1.2.2"


def test_non_breaking_always_ok():
    ok, _ = s.decide(False, "1.0.0", "1.0.0")
    assert ok is True


def test_post1_breaking_without_major_fails():
    ok, msg = s.decide(True, "1.2.0", "1.3.0")
    assert ok is False
    assert "major" in msg.lower()


def test_post1_breaking_with_major_passes():
    ok, _ = s.decide(True, "1.2.0", "2.0.0")
    assert ok is True


def test_pre1_breaking_with_minor_passes():
    ok, _ = s.decide(True, "0.4.0", "0.5.0")
    assert ok is True


def test_pre1_breaking_patch_only_fails():
    ok, msg = s.decide(True, "0.4.0", "0.4.1")
    assert ok is False
    assert "minor" in msg.lower()
