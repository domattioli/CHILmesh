#!/usr/bin/env python3
"""Verify WNAT_Hagen SHA-256 pin matches expected value.

Exits non-zero if mismatch. Used by benchmark scripts to ensure
reproducible input across runs."""

import hashlib
import json
import sys
from pathlib import Path


def verify_mesh_pin(
    mesh_path: str = "/tmp/admesh-domains/registry_data/meshes/WNAT_Hagen.14",
    pin_file: str = "output/benchmark.json"
) -> bool:
    """Verify mesh file matches pinned SHA-256.

    Args:
        mesh_path: Path to the WNAT_Hagen mesh file
        pin_file: Path to benchmark.json containing pinned hash

    Returns:
        True if hash matches, False otherwise
    """
    # Compute actual hash
    sha256 = hashlib.sha256()
    try:
        with open(mesh_path, "rb") as f:
            for chunk in iter(lambda: f.read(65536), b""):
                sha256.update(chunk)
    except FileNotFoundError:
        print(f"ERROR: Mesh file not found: {mesh_path}", file=sys.stderr)
        return False

    actual_hash = sha256.hexdigest()

    # Load pinned hash
    try:
        with open(pin_file, "r") as f:
            data = json.load(f)
            expected_hash = data["metadata"]["mesh_sha256"]
    except (FileNotFoundError, KeyError, json.JSONDecodeError) as e:
        print(f"ERROR: Could not read pinned hash from {pin_file}: {e}", file=sys.stderr)
        return False

    # Compare
    if actual_hash == expected_hash:
        print(f"✓ Mesh pin verified: {actual_hash[:12]}...", file=sys.stderr)
        return True
    else:
        print(f"ERROR: Mesh pin mismatch!", file=sys.stderr)
        print(f"  Expected: {expected_hash}", file=sys.stderr)
        print(f"  Actual:   {actual_hash}", file=sys.stderr)
        return False


if __name__ == "__main__":
    success = verify_mesh_pin()
    sys.exit(0 if success else 1)
