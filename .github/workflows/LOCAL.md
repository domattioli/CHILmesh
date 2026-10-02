# LOCAL.md — repo-local workflow registry (spec-010 v2.3)

Workflows listed here are intentionally repo-local (not DomI-managed copies). Adding a
new local workflow requires a row here in the same PR — unlisted local
workflows fail the workflow-conformance gate.

| Workflow | Justification |
|---|---|
| `python-package.yml` | full cross-OS test matrix incl. macOS lanes (main-push gated) — macOS gating is repo-local by design (spec-010 v2.2 rule 8) |
| `publish-pypi.yml` | PyPI release, tag-triggered — repo-specific release pipeline |
| `build-cpp-wheels.yml` | build-only wheel validation for the chilmesh_cpp binary backend on Linux/macOS/Windows (workflow_dispatch, artifacts only, no publish) — repo-specific packaging (#229, #256) |
| `publish-cpp-wheels.yml` | release/dispatch-gated PyPI publish of chilmesh_cpp wheels (Linux/macOS/Windows) + sdist, then a post-publish install smoke on 3 OSes (no push trigger; `environment: pypi` protected; Trusted Publishing/OIDC) — repo-specific packaging (#256) |
| `code-smell.yml` | ruff F + bandit MEDIUM gate, style/complexity advisory, weekly cron — repo-local static-analysis lane (#267) |
