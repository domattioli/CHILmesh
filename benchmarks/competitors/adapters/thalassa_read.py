"""Standalone adapter: open a mesh with Thalassa and print one JSON line.

Run ONLY in a subprocess inside the Thalassa environment (see envs/). The call
shape and the input format Thalassa documents as its default path are unverified
until research.md (T002/T005); failure is recorded as ``failed``.
"""
import argparse
import json
import sys
import time


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    args = ap.parse_args()

    import thalassa
    from thalassa import api

    t0 = time.perf_counter()
    ds = api.open_dataset(args.mesh)
    parse_s = time.perf_counter() - t0
    dims = dict(getattr(ds, "sizes", {}) or {})
    print(json.dumps({"tool": "thalassa", "version": getattr(thalassa, "__version__", "unknown"),
                      "parse_s": parse_s, "n_nodes": dims.get("node"), "n_elems": dims.get("element"),
                      "input_format": "as accepted by thalassa.api.open_dataset"}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
