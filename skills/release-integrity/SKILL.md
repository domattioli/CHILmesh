---
name: release-integrity
version: 1.0
description: Enforce package-release integrity — fail breaking public-API changes without a (major, or 0.x minor) bump and malformed releases (tag≠version, missing changelog, twine). Use in consumer CI via the spec-010 release-integrity workflow.
benchmark: malformed_releases_per_10_publishes
---

# release-integrity

Mechanical enforcement of DomI's Package Release & SemVer Governance article.

## Components
- `scripts/semver_core.py` — 0.x-aware SemVer `decide`.
- `scripts/changelog.py` — Keep-a-Changelog checks.
- `scripts/release_integrity.py` — CLI: `api-gate` (warn-first), `release-check` (enforce), `pr-changelog` (warn).
- `tests/eval/` — labeled scenario corpus scored for classification accuracy.

## Use
Consumers receive `templates/workflows/release-integrity.yml` via sync. Set repo
variables `RELEASE_INTEGRITY_PACKAGE` (import name) and `RELEASE_INTEGRITY_SRC`
(default `src`). The API gate is warn-first; a repo graduates to `--mode enforce`
when its `false_positive_gate_failures_per_release_cycle` is acceptably low.

## Version History
- **1.0 (2026-06-20):** initial — gate + release-consistency + eval; promoted from the CHILmesh `api_semver_gate` spike.
