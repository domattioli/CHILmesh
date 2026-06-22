"""Labeled release scenarios for the release-integrity eval. Expected verdicts
are the gate's correct classification; the runner scores against these."""

CASES = [
    {"name": "additive-minor", "kind": "api", "breaking": False,
     "old": "1.2.0", "new": "1.3.0", "mode": "enforce", "expect_code": 0},
    {"name": "breaking-major", "kind": "api", "breaking": True,
     "old": "1.2.0", "new": "2.0.0", "mode": "enforce", "expect_code": 0},
    {"name": "breaking-no-major-enforce", "kind": "api", "breaking": True,
     "old": "1.2.0", "new": "1.3.0", "mode": "enforce", "expect_code": 1},
    {"name": "breaking-no-major-warn", "kind": "api", "breaking": True,
     "old": "1.2.0", "new": "1.3.0", "mode": "warn", "expect_code": 0},
    {"name": "0x-breaking-minor", "kind": "api", "breaking": True,
     "old": "0.4.0", "new": "0.5.0", "mode": "enforce", "expect_code": 0},
    {"name": "0x-breaking-patch", "kind": "api", "breaking": True,
     "old": "0.4.0", "new": "0.4.1", "mode": "enforce", "expect_code": 1},
    {"name": "release-consistent", "kind": "release", "tag": "v1.2.0", "version": "1.2.0",
     "changelog": "## [1.2.0] — 2026-05-24\n\n- x\n", "twine_ok": True, "expect_code": 0},
    {"name": "release-tag-mismatch", "kind": "release", "tag": "v1.2.1", "version": "1.2.0",
     "changelog": "## [1.2.0] — 2026-05-24\n\n- x\n", "twine_ok": True, "expect_code": 1},
    {"name": "release-missing-changelog", "kind": "release", "tag": "v1.2.0", "version": "1.2.0",
     "changelog": "# Changelog\n", "twine_ok": True, "expect_code": 1},
]
