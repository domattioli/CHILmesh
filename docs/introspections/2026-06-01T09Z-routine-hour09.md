# CHILmesh Introspection — 2026-06-01T09Z (hour-09 routine)

session_id: session_018MwZyFGcgiSN6dHaxJbfsK
repo: CHILmesh
branch: daily-maintenance
session_sha_start: c7be8bc
session_sha_end: f0997ae
commits: 4
files_changed: 4
insertions: 70
deletions: 4
issues_closed: 3
wall_clock_min: ~25

## Work done

- #185 domi-sync: .domi-pin 89a93bc → 9a2ec97 (pin-only, no local skills)
- #154 docs: created CONTEXT.md — CHILmesh vocabulary glossary + .chil design conflict notes
- #155 feat: added lifecycle stages table to BENCHMARK.md; confirmed stages locked
- fix: pyproject.toml TOML syntax `{"": "src"}` → `{"" = "src"}` (unblocked pip install on fresh containers)

## Pre-flight

- branch_policy_conflict: NO (daily-maintenance correct, harness injected claude/lucid-fermi-PyiIM — ignored per CLAUDE.md)
- mcp_scope_gap: NO
- label_scheme_mismatch: NO (priority:now correctly filtered)

## Pain points

```yaml
pains:
  - pain: >
      git commit signing-server returned status 400 "missing source" — 
      all 4 commits this session needed commit.gpgsign=false workaround
    frequency: recurring-across-sessions
    severity: medium
    evidence: "cc7b78d c1cf5b0 a3ce406 f0997ae — all signed with gpgsign=false"
    existing_skill_should_have_caught_it: git-push-fallback
    missing_skill_would_have_prevented_it: none — workaround documented in CLAUDE.md; root cause is infra
    domi_issue: "#119"
    saved_time_estimate_min: 2

  - pain: >
      ensure-test-venv.sh failed with TOML parse error on pyproject.toml
      {"": "src"} — JSON colon syntax not valid TOML; pip 24 vendored tomli
      rejects it. Root cause: bug in CHILmesh pyproject.toml, not the skill.
    frequency: once
    severity: medium
    evidence: "tomllib.TOMLDecodeError: Expected '=' after key at line 66 col 18"
    existing_skill_should_have_caught_it: none
    missing_skill_would_have_prevented_it: none — pyproject lint would have caught it
    domi_issue: null
    saved_time_estimate_min: 3

  - pain: >
      sync-from-domi plugin not installed at container start — manual pin
      update required; sync-from-domi automated drift detection skipped
    frequency: recurring-across-sessions
    severity: medium
    evidence: "instructions_on_start.sh: ⚠ sync-from-domi not installed"
    existing_skill_should_have_caught_it: plugin-install-with-vendored-fallback
    missing_skill_would_have_prevented_it: plugin-install-with-vendored-fallback
    domi_issue: "#114"
    saved_time_estimate_min: 5

  - pain: >
      caveman plugin not installed at container start — ultra mode emulated
      inline from SKILL.md; no slash-command available
    frequency: recurring-across-sessions
    severity: low
    evidence: "ToolSearch: caveman not found; emulated from plugins/caveman/skills/caveman/SKILL.md"
    existing_skill_should_have_caught_it: plugin-install-with-vendored-fallback
    missing_skill_would_have_prevented_it: plugin-install-with-vendored-fallback
    domi_issue: "#114"
    saved_time_estimate_min: 1
```

## Introspection metadata

- pain_points_captured: 4
- recurring_cross_session: 3 (signing #119, sync-from-domi #114, caveman #114)
- new_single_occurrence: 1 (pyproject.toml TOML bug — CHILmesh-local fix, not DomI skill gap)
- votes_cast: 2 (#119, #114)
- new_requests_filed: 0 (frequency gate — all map to existing issues)
- corpus_written: docs/introspections/2026-06-01T09Z-routine-hour09.md
- pre_flight: branch-policy: no conflict, mcp-scope: ok, label-scheme: ok
