# Security Policy

CHILmesh is a Python library for 2D triangular, quadrilateral, and mixed-element
mesh generation and manipulation. It includes a hybrid Python/C++ backend and is
available on PyPI as `chilmesh`. This security policy covers the published package
and source repository.

## Supported versions

| Version | Supported |
|---|---|
| `1.2.x` (current PyPI release) | ✅ — receives fixes |
| `< 1.2` (older tags) | ⚠️ — best-effort only |

The current stable release is `1.2.2`. We prioritize security fixes for the latest
`1.2.x` series; patches for older major/minor versions are made on a best-effort
basis and are not guaranteed.

## Reporting a vulnerability

**Please do not open a public issue for security problems.**

Report privately through GitHub's
[private vulnerability reporting](https://github.com/domattioli/CHILmesh/security/advisories/new)
for this repository. If that is unavailable to you, contact a maintainer via
their GitHub profile ([@domattioli](https://github.com/domattioli)) and request
a private channel.

When reporting, please include:

- A description of the issue and its impact.
- Steps to reproduce (a minimal proof of concept where practical).
- Affected components (library API, C++ extension, test utilities, CI).
- Any suggested remediation.

We aim to acknowledge a report within a reasonable window and will coordinate a
fix and disclosure timeline with you.

## Scope notes

- **C++ extension memory safety.** CHILmesh's hybrid Python/C++ backend is in scope
  for security review. Crashes, memory corruption, or undefined behavior triggered
  by malformed mesh input (e.g. invalid vertex coordinates, degenerate elements,
  topological violations) are treated as security concerns, as they could represent
  memory safety vulnerabilities.
- **Secrets.** No tokens, API keys, or production credentials are committed to the
  repository. All write surfaces (CI, deployment, PyPI publishing) read secrets
  from the environment / GitHub Actions secrets only.
- **Attestation / signing keys.** No signing keys are committed. Releases use
  PyPI's trusted publishing (OpenID Connect) and do not require locally stored
  credentials.

## Out of scope

- Findings that require a compromised maintainer machine or stolen credentials.
- Denial of service against optional external services (e.g. HuggingFace Hub
  connectivity for domain registry lookups).
- Dependency vulnerabilities in transitive (not direct) packages, unless we can
  confirm CHILmesh's usage amplifies the risk.
