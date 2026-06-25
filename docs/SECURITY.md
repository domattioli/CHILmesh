# Security Policy

## Supported versions

CHILmesh is maintained on a rolling basis; only the latest published release
receives security fixes.

| Version | Supported |
|---------|-----------|
| 1.0.x   | ✅ |
| < 1.0   | ❌ |

## Reporting a vulnerability

**Please do not open a public issue or pull request for security
vulnerabilities.** Public disclosure before a fix is available puts every user
at risk.

Report privately through one of:

1. **GitHub private vulnerability reporting** (preferred) — on this repository,
   go to the **Security** tab → **Report a vulnerability**. This opens a private
   advisory visible only to you and the maintainers.
2. Contact the maintainer directly via the GitHub profile
   [@domattioli](https://github.com/domattioli) if private reporting is
   unavailable.

Please include, where possible:

- the affected version / commit;
- a minimal reproduction (e.g. a crafted `fort.14` / `.2dm` input, a code
  snippet, and the observed vs. expected behavior);
- the impact you believe it has (crash / DoS, memory corruption, data
  tampering, etc.).

## What to expect

- **Acknowledgement:** we aim to respond within a few days.
- **Assessment:** we triage severity and confirm the issue.
- **Fix & disclosure:** once a fix is ready we coordinate a release and, where
  appropriate, publish a GitHub Security Advisory (GHSA) crediting the reporter
  (unless you prefer to remain anonymous).

## Scope notes

CHILmesh parses untrusted mesh files (`fort.14`, `.2dm`). Parser-level
robustness (size/record caps, dtype and index validation) is in scope. The
optional C++/Rust backends are **experimental**; to run the pure-Python backend
only, set `CHILMESH_BACKEND=python`.
