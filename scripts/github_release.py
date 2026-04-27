#!/usr/bin/env python3
"""GitHub Release Skill - Create releases on GitHub with automatic credential detection.

Non-interactive skill that:
1. Auto-detects GitHub credentials
2. Auto-detects version from pyproject.toml
3. Auto-detects repo from git remote
4. Extracts release notes from CHANGELOG.md
5. Creates release via gh CLI
6. Reports success with release URL

Usage:
  github_release.py [--version VERSION] [--repo owner/repo]

Credential Options:
  1. Local gh CLI:
     $ gh auth login
     Then run: python scripts/github_release.py

  2. Environment variable (in Claude Code):
     Use: /update-config to set GITHUB_TOKEN securely
     Then run: python scripts/github_release.py

  3. Environment variable (terminal):
     $ GITHUB_TOKEN=your_token python scripts/github_release.py

Exit codes:
  0 = success
  1 = failure
"""

import subprocess
import sys
import re
import os
from pathlib import Path
from typing import Optional, Tuple


def run(cmd: list[str], check: bool = False) -> Tuple[int, str, str]:
    """Run command, return (code, stdout, stderr)."""
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        return result.returncode, result.stdout.strip(), result.stderr.strip()
    except Exception as e:
        return 1, "", str(e)


def detect_github_credentials() -> Tuple[bool, Optional[str]]:
    """Detect GitHub credentials non-interactively."""
    # Check gh auth status
    rc, _, _ = run(["gh", "auth", "status"])
    if rc == 0:
        return True, "gh CLI authenticated"

    # Check GITHUB_TOKEN env var
    if os.getenv("GITHUB_TOKEN"):
        return True, "GITHUB_TOKEN env var"

    # Check gh config files
    for gh_config in ["~/.config/gh/hosts.yml", "~/.gh/hosts.yml"]:
        path = Path(gh_config).expanduser()
        if path.exists() and "oauth_token" in path.read_text():
            return True, f"gh config at {gh_config}"

    # Try to validate via gh API
    rc, _, _ = run(["gh", "api", "user"])
    if rc == 0:
        return True, "gh API accessible"

    return False, None


def extract_version(arg_version: Optional[str]) -> Optional[str]:
    """Auto-detect version from pyproject.toml or use provided."""
    if arg_version:
        return arg_version

    try:
        with open("pyproject.toml") as f:
            match = re.search(r'version\s*=\s*["\']([^"\']+)["\']', f.read())
            return match.group(1) if match else None
    except:
        return None


def detect_repo() -> Optional[str]:
    """Get repo from git remote origin."""
    rc, stdout, _ = run(["git", "remote", "get-url", "origin"])
    if rc != 0:
        return None

    # Parse owner/repo from git URL
    match = re.search(r'github\.com[:/]([^/]+)/([^/\.]+)(?:\.git)?$', stdout)
    return f"{match.group(1)}/{match.group(2)}" if match else None


def extract_changelog_section(version: str) -> str:
    """Extract release notes for version from CHANGELOG.md."""
    try:
        with open("CHANGELOG.md") as f:
            content = f.read()

        # Pattern: ## VERSION or ## [VERSION] or ## vVERSION
        pattern = rf"^##\s+(?:\[?v?)?{re.escape(version)}(?:\])?.*?(?=^##\s|$)"
        match = re.search(pattern, content, re.MULTILINE | re.DOTALL)

        if match:
            notes = match.group(0).strip()
            # Remove header, keep content
            lines = notes.split('\n', 1)
            return lines[1].strip() if len(lines) > 1 else f"Release {version}"

        return f"Release {version}"
    except:
        return f"Release {version}"


def create_release(version: str, repo: str, notes: str) -> Tuple[bool, str]:
    """Create GitHub release."""
    cmd = ["gh", "release", "create", f"v{version}", "--title", version, "--notes", notes]
    if repo:
        cmd.extend(["--repo", repo])

    rc, stdout, stderr = run(cmd)

    if rc == 0:
        # Extract URL from output
        for line in stdout.split('\n'):
            if 'github.com' in line and 'releases' in line:
                return True, line.strip()
        return True, stdout

    return False, stderr


def main() -> int:
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--version")
    parser.add_argument("--repo")
    args = parser.parse_args()

    # 1. Check credentials
    has_creds, creds_src = detect_github_credentials()
    if not has_creds:
        print("No GitHub credentials found. Run: gh auth login")
        return 1

    # 2. Get version
    version = extract_version(args.version)
    if not version:
        print("Could not determine version. Set in pyproject.toml or use --version")
        return 1

    # 3. Get repo
    repo = args.repo or detect_repo()
    if not repo:
        print("Could not determine repo. Use git remote or --repo owner/repo")
        return 1

    # 4. Extract release notes
    notes = extract_changelog_section(version)

    # 5. Create release
    success, url = create_release(version, repo, notes)
    if not success:
        print(f"Release creation failed: {url}")
        return 1

    print(f"✅ GitHub Release: {url}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
