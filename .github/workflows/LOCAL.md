# LOCAL.md — repo-local workflow registry (spec-010 v2.3)

Workflows listed here are intentionally repo-local (not DomI-managed copies). Adding a
new local workflow requires a row here in the same PR — unlisted local
workflows fail the workflow-conformance gate.

| Workflow | Justification |
|---|---|
| `python-package.yml` | full cross-OS test matrix incl. macOS lanes (main-push gated) — macOS gating is repo-local by design (spec-010 v2.2 rule 8) |
| `publish-pypi.yml` | PyPI release, tag-triggered — repo-specific release pipeline |
