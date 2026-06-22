# release-integrity benchmark

**Primary:** `malformed_releases_per_10_publishes` (lower-better) — releases failing any
release-consistency gate (tag≠version, missing/stale changelog, twine), per 10 publishes.
**Secondary (qualitative):** `false_positive_gate_failures_per_release_cycle` — spurious gate
failures; the signal for warn→enforce graduation. C6 eval FP rate is the early in-repo proxy.

| Version | Date | Primary (observed) | Delta | Notes |
|---------|------|--------------------|-------|-------|
| 1.0 | 2026-06-20 | not-measured | — | initial; accrues as consumers publish |
