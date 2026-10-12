#!/usr/bin/env python3
"""Download SHA-256-pinned real ADCIRC/SMS files and parse them with chilmesh readers.

Usage: python scripts/check_adcirc_testsuite.py [--stage B1|B2|B3|all] [--cache DIR]

Exit 0 only when every selected manifest entry downloads, is not a Git LFS
pointer, matches its sha256, and parses with the reader registered for its
format. Downloaded files are data only and are never executed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MANIFEST = Path(__file__).resolve().parent / "adcirc_testsuite_manifest.json"
LFS_HEAD = b"version https://git-lfs.github.com/spec/"


def _read_2dm(path: Path) -> str:
    from chilmesh import CHILmesh
    m = CHILmesh.read_from_2dm(path, compute_layers=False)
    return f"{m.n_verts} nodes, {m.n_elems} elems, {len(m.boundary_segments)} nodestrings"


def _read_output(path: Path) -> str:
    from chilmesh import read_adcirc_output
    s = read_adcirc_output(path)
    if s.dropped_records:
        raise ValueError(f"{s.dropped_records} record(s) dropped (truncated file)")
    return f"{s.kind} {s.location} {s.values.shape}, t=[{s.times[0]:g}..{s.times[-1]:g}]"


def _read_output_nc(path: Path) -> str:
    from chilmesh import read_adcirc_output_nc
    s = read_adcirc_output_nc(path)
    return f"{s.kind} {s.location} {s.values.shape}, t=[{s.times[0]:g}..{s.times[-1]:g}]"


_B2_ASCII = ["fort.61", "fort.62", "fort.63", "fort.64", "maxele.63", "maxvel.63"]

# format -> reader(path) -> description; raises on parse failure.
READERS = {"2dm": _read_2dm}
READERS.update({f: _read_output for f in _B2_ASCII})
READERS.update({f + ".nc": _read_output_nc for f in _B2_ASCII})


def _fetch(url: str, dest: Path) -> None:
    req = urllib.request.Request(url, headers={"User-Agent": "chilmesh-adcirc-gate/1"})
    with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as f:
        while chunk := r.read(1 << 20):
            f.write(chunk)


def check_entry(e: dict, cache: Path) -> str | None:
    """Return None on success or a failure message."""
    dest = cache / e["name"]
    try:
        if not dest.exists():
            _fetch(e["url"], dest)
        data = dest.read_bytes()
    except Exception as exc:  # network / IO
        return f"download failed: {exc}"
    if data.startswith(LFS_HEAD):
        dest.unlink()
        return "file is a Git LFS pointer (use the media.githubusercontent.com URL)"
    digest = hashlib.sha256(data).hexdigest()
    if digest != e["sha256"]:
        dest.unlink()
        return f"sha256 mismatch: expected {e['sha256']} got {digest}"
    if e.get("size") is not None and len(data) != e["size"]:
        return f"size mismatch: expected {e['size']} got {len(data)}"
    reader = READERS.get(e["format"])
    if reader is None:
        return f"no reader registered for format {e['format']!r}"
    try:
        info = reader(dest)
    except Exception as exc:
        return f"parse failed: {type(exc).__name__}: {exc}"
    print(f"  parsed: {info}")
    return None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--stage", default="all", choices=["B1", "B2", "B3", "all"])
    ap.add_argument("--cache", type=Path, default=ROOT / ".adcirc_testsuite")
    args = ap.parse_args(argv)
    entries = json.loads(MANIFEST.read_text())["entries"]
    if args.stage != "all":
        entries = [e for e in entries if e["stage"] == args.stage]
    args.cache.mkdir(parents=True, exist_ok=True)
    failures = 0
    for e in entries:
        print(f"[{e['stage']}] {e['name']} ({e['format']})")
        err = check_entry(e, args.cache)
        if err:
            failures += 1
            print(f"  FAIL: {err}")
        else:
            print("  PASS")
    print(f"{len(entries) - failures}/{len(entries)} passed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
