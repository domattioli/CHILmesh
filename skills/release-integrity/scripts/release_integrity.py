"""release-integrity CLI: api-gate, release-check, pr-changelog.

Pure decision helpers (run_api_gate, run_release_check) take gathered inputs so
they unit-test without I/O. main() gathers git/griffe/twine state and dispatches.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from typing import List, Optional, Tuple

import semver_core
import changelog as _changelog


def run_api_gate(breaking: bool, old_version: str, new_version: str, mode: str) -> Tuple[int, str]:
    ok, msg = semver_core.decide(breaking, old_version, new_version)
    if ok:
        return 0, msg
    if mode == "enforce":
        return 1, msg
    return 0, f"WARN (warn-first mode): {msg}"


def run_release_check(tag: str, version: str, changelog_text: str, twine_ok: bool) -> Tuple[int, str]:
    problems: List[str] = []
    if tag.lstrip("v") != version:
        problems.append(f"tag {tag!r} != version {version!r}")
    if not _changelog.has_dated_section(changelog_text, version):
        problems.append(f"no dated CHANGELOG section for {version}")
    if not twine_ok:
        problems.append("twine check failed")
    if problems:
        return 1, "Release-consistency FAILED: " + "; ".join(problems)
    return 0, "Release-consistency OK."


def _run(cmd: List[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True)


def _latest_tag() -> Optional[str]:
    p = _run(["git", "describe", "--tags", "--abbrev=0"])
    return p.stdout.strip() if p.returncode == 0 else None


def _version_at(ref: Optional[str]) -> str:
    if ref is None:
        with open("pyproject.toml") as f:
            text = f.read()
    else:
        text = _run(["git", "show", f"{ref}:pyproject.toml"]).stdout
    return semver_core.extract_version(text)


def _cmd_api_gate(args) -> int:
    tag = _latest_tag()
    if tag is None:
        print("No baseline tag; skipping API gate.")
        return 0
    g = _run([sys.executable, "-m", "griffe", "check", args.package, "-s", args.src, "-a", tag])
    if g.stdout:
        print(g.stdout)
    if g.stderr:
        print(g.stderr, file=sys.stderr)
    code, msg = run_api_gate(g.returncode != 0, _version_at(tag), _version_at(None), args.mode)
    print(msg)
    return code


def _cmd_release_check(args) -> int:
    version = _version_at(None)
    tag = _run(["git", "describe", "--tags", "--exact-match"]).stdout.strip() or f"v{version}"
    with open("CHANGELOG.md") as f:
        cl = f.read()
    twine_ok = True
    if args.pypi:
        twine_ok = _run([sys.executable, "-m", "twine", "check"]).returncode == 0
    code, msg = run_release_check(tag, version, cl, twine_ok)
    print(msg)
    return code


def _cmd_pr_changelog(args) -> int:
    diff = _run(["git", "diff", "--name-only", f"{args.base}...HEAD"]).stdout
    touched_src = any(line.startswith(args.src + "/") for line in diff.splitlines())
    if not touched_src:
        return 0
    with open("CHANGELOG.md") as f:
        cl = f.read()
    if not _changelog.has_unreleased_bullet(cl):
        print("WARN: src/ changed but no '## [Unreleased]' bullet added.")
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(prog="release-integrity")
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("api-gate")
    a.add_argument("--package", required=True)
    a.add_argument("--src", default="src")
    a.add_argument("--mode", choices=["warn", "enforce"], default="warn")
    a.set_defaults(fn=_cmd_api_gate)
    r = sub.add_parser("release-check")
    r.add_argument("--pypi", action="store_true")
    r.set_defaults(fn=_cmd_release_check)
    c = sub.add_parser("pr-changelog")
    c.add_argument("--src", default="src")
    c.add_argument("--base", default="origin/main")
    c.set_defaults(fn=_cmd_pr_changelog)
    args = p.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
