from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
import changelog as c  # noqa: E402

DATED = "# Changelog\n\n## [1.2.0] — 2026-05-24\n\n### Added\n- a thing\n"
UNRELEASED_FILLED = "## [Unreleased]\n\n### Added\n- new bullet\n\n## [1.2.0] — 2026-05-24\n"
UNRELEASED_EMPTY = "## [Unreleased]\n\n## [1.2.0] — 2026-05-24\n"


def test_has_dated_section_true():
    assert c.has_dated_section(DATED, "1.2.0") is True


def test_has_dated_section_missing_version():
    assert c.has_dated_section(DATED, "9.9.9") is False


def test_has_dated_section_requires_date():
    assert c.has_dated_section("## [1.2.0]\n\n- x\n", "1.2.0") is False


def test_has_unreleased_bullet_true():
    assert c.has_unreleased_bullet(UNRELEASED_FILLED) is True


def test_has_unreleased_bullet_empty():
    assert c.has_unreleased_bullet(UNRELEASED_EMPTY) is False


def test_has_unreleased_bullet_absent():
    assert c.has_unreleased_bullet(DATED) is False
