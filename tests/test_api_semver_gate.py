from __future__ import annotations

import pathlib
import sys

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
import api_semver_gate as g  # noqa: E402


def test_parse_major():
    assert g.parse_major("1.2.2") == 1
    assert g.parse_major("2.0.0") == 2


def test_parse_major_rejects_garbage():
    with pytest.raises(ValueError):
        g.parse_major("not-a-version")


def test_decide_breaking_no_major_bump_fails():
    ok, msg = g.decide(True, 1, 1)
    assert ok is False
    assert "without a major bump" in msg


def test_decide_breaking_with_major_bump_passes():
    ok, _ = g.decide(True, 1, 2)
    assert ok is True


def test_decide_non_breaking_passes():
    ok, _ = g.decide(False, 1, 1)
    assert ok is True
