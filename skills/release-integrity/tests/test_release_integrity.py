from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
import release_integrity as ri  # noqa: E402


def test_api_gate_warn_mode_never_fails():
    code, msg = ri.run_api_gate(True, "1.2.0", "1.3.0", mode="warn")
    assert code == 0
    assert "WARN" in msg


def test_api_gate_enforce_mode_fails_on_violation():
    code, msg = ri.run_api_gate(True, "1.2.0", "1.3.0", mode="enforce")
    assert code == 1


def test_api_gate_enforce_passes_when_major_bumped():
    code, _ = ri.run_api_gate(True, "1.2.0", "2.0.0", mode="enforce")
    assert code == 0


def test_release_check_passes_when_consistent():
    code, _ = ri.run_release_check(
        tag="v1.2.0", version="1.2.0",
        changelog_text="## [1.2.0] — 2026-05-24\n\n- x\n", twine_ok=True)
    assert code == 0


def test_release_check_fails_on_tag_mismatch():
    code, msg = ri.run_release_check(
        tag="v1.2.1", version="1.2.0",
        changelog_text="## [1.2.0] — 2026-05-24\n\n- x\n", twine_ok=True)
    assert code == 1
    assert "tag" in msg.lower()


def test_release_check_fails_on_missing_changelog():
    code, msg = ri.run_release_check(
        tag="v1.2.0", version="1.2.0", changelog_text="# Changelog\n", twine_ok=True)
    assert code == 1
    assert "changelog" in msg.lower()


def test_release_check_fails_on_twine():
    code, msg = ri.run_release_check(
        tag="v1.2.0", version="1.2.0",
        changelog_text="## [1.2.0] — 2026-05-24\n\n- x\n", twine_ok=False)
    assert code == 1
    assert "twine" in msg.lower()
