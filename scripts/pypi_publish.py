#!/usr/bin/env python3
"""PyPI Publish Skill - Publish packages to PyPI with automatic credential detection.

Non-interactive skill that:
1. Auto-detects PyPI credentials
2. Auto-detects package name/version from pyproject.toml
3. Auto-builds if dist/ missing
4. Validates distribution files
5. Publishes via twine with retry logic
6. Verifies on PyPI
7. Reports success with package URL

Usage:
  pypi_publish.py [--version VERSION] [--repo /path/to/project]

Credential Options:
  1. PyPI config file (~/.pypirc):
     [pypi]
     repository = https://upload.pypi.org/legacy/
     username = __token__
     password = pypi-your_token_here
     Then run: python scripts/pypi_publish.py

  2. Environment variable (in Claude Code):
     Use: /update-config to set PYPI_TOKEN securely
     Then run: python scripts/pypi_publish.py

  3. Environment variable (terminal):
     $ PYPI_TOKEN=pypi-your_token_here python scripts/pypi_publish.py

Exit codes:
  0 = success
  1 = failure
"""

import subprocess
import sys
import re
import os
import time
import urllib.request
from pathlib import Path
from typing import Optional, Tuple


def run(cmd: list[str]) -> Tuple[int, str, str]:
    """Run command, return (code, stdout, stderr)."""
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        return result.returncode, result.stdout.strip(), result.stderr.strip()
    except Exception as e:
        return 1, "", str(e)


def detect_pypi_credentials() -> Tuple[bool, Optional[str]]:
    """Detect PyPI credentials non-interactively."""
    # Check PYPI_TOKEN env var
    if os.getenv("PYPI_TOKEN"):
        return True, "PYPI_TOKEN env var"

    # Check ~/.pypirc
    pypirc = Path("~/.pypirc").expanduser()
    if pypirc.exists():
        try:
            content = pypirc.read_text()
            if "password" in content and "pypi-" in content:
                return True, "~/.pypirc"
        except:
            pass

    # Try twine check (validates tokens)
    if Path("dist").exists():
        rc, _, _ = run(["twine", "check", "dist/*"])
        if rc == 0:
            return True, "twine check passed"

    return False, None


def extract_pyproject_data(project_root: str) -> Tuple[Optional[str], Optional[str]]:
    """Extract package name and version from pyproject.toml."""
    try:
        pyproject = Path(project_root) / "pyproject.toml"
        content = pyproject.read_text()

        name_match = re.search(r'name\s*=\s*["\']([^"\']+)["\']', content)
        version_match = re.search(r'version\s*=\s*["\']([^"\']+)["\']', content)

        name = name_match.group(1) if name_match else None
        version = version_match.group(1) if version_match else None

        return name, version
    except:
        return None, None


def ensure_build(project_root: str) -> Tuple[bool, str]:
    """Build distributions if missing."""
    dist_dir = Path(project_root) / "dist"

    if not dist_dir.exists() or not list(dist_dir.glob("*.whl")) or not list(dist_dir.glob("*.tar.gz")):
        cwd = os.getcwd()
        try:
            os.chdir(project_root)
            rc, stdout, stderr = run(["python", "-m", "build"])
            os.chdir(cwd)
            if rc != 0:
                return False, f"Build failed: {stderr}"
        except Exception as e:
            os.chdir(cwd)
            return False, str(e)

    return True, "Build complete"


def validate_distributions(project_root: str, pkg_name: str, version: str) -> Tuple[bool, str]:
    """Validate distribution files exist."""
    dist_dir = Path(project_root) / "dist"

    # Normalize package name (replace - with _)
    normalized = pkg_name.replace("-", "_")

    # Check for .whl and .tar.gz
    wheels = list(dist_dir.glob(f"{normalized}-{version}-py3-none-any.whl"))
    tarballs = list(dist_dir.glob(f"{pkg_name}-{version}.tar.gz"))

    if not wheels and not tarballs:
        return False, f"No distributions found for {pkg_name} {version}"

    return True, f"Found {len(wheels)} wheels, {len(tarballs)} tarballs"


def upload_with_retry(project_root: str, pkg_name: str, version: str) -> Tuple[bool, str]:
    """Upload to PyPI with retry logic (2s, 4s, 8s backoff)."""
    cmd = ["twine", "upload", f"dist/{pkg_name}-{version}*", "--skip-existing"]

    backoffs = [2, 4, 8]
    for attempt in range(len(backoffs) + 1):
        rc, stdout, stderr = run(cmd)

        if rc == 0:
            return True, stdout

        if attempt < len(backoffs):
            time.sleep(backoffs[attempt])
            continue

        return False, f"Upload failed after retries: {stderr}"


def verify_on_pypi(pkg_name: str, version: str) -> Tuple[bool, str]:
    """Verify package is on PyPI."""
    url = f"https://pypi.org/project/{pkg_name}/{version}/"

    try:
        with urllib.request.urlopen(url) as response:
            if response.status == 200:
                return True, url
            return False, f"PyPI returned {response.status}"
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return False, "Package not found on PyPI yet (may take a moment)"
        return False, f"PyPI error: {e}"
    except Exception as e:
        return False, str(e)


def main() -> int:
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--version")
    parser.add_argument("--repo", default=".")
    args = parser.parse_args()

    project_root = args.repo

    # 1. Check credentials
    has_creds, creds_src = detect_pypi_credentials()
    if not has_creds:
        print("No PyPI credentials found. Create ~/.pypirc or set PYPI_TOKEN")
        return 1

    # 2. Extract package info
    pkg_name, version = extract_pyproject_data(project_root)
    if not pkg_name or not (version or args.version):
        print("Could not determine package name/version. Check pyproject.toml")
        return 1

    version = args.version or version

    # 3. Ensure build
    success, msg = ensure_build(project_root)
    if not success:
        print(f"Build failed: {msg}")
        return 1

    # 4. Validate distributions
    success, msg = validate_distributions(project_root, pkg_name, version)
    if not success:
        print(f"Validation failed: {msg}")
        return 1

    # 5. Upload with retry
    success, msg = upload_with_retry(project_root, pkg_name, version)
    if not success:
        print(f"Upload failed: {msg}")
        return 1

    # 6. Verify on PyPI
    success, url = verify_on_pypi(pkg_name, version)
    if not success:
        print(f"Verification failed: {url}")
        return 1

    print(f"✅ PyPI Published: {url}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
