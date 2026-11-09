"""Standalone adapter: read a fort.14 with adcircpy and print one JSON line.

Run ONLY in a subprocess inside the adcircpy environment (see envs/). adcircpy is
GPL-licensed; this script is the only place it is imported, and CHILmesh never
imports this file. API call shape is unverified until research.md (T002/T005)
confirms it; a mismatch surfaces as a non-zero exit, recorded as ``failed``.
"""
import argparse
import json
import sys
import time


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    args = ap.parse_args()

    import adcircpy

    t0 = time.perf_counter()
    mesh = adcircpy.AdcircMesh.open(args.mesh, crs=None)
    parse_s = time.perf_counter() - t0
    n_nodes = len(mesh.coords) if hasattr(mesh, "coords") else None
    n_elems = len(mesh.elements.elements) if hasattr(mesh, "elements") and hasattr(mesh.elements, "elements") else None
    print(json.dumps({"tool": "adcircpy", "version": getattr(adcircpy, "__version__", "unknown"),
                      "parse_s": parse_s, "n_nodes": n_nodes, "n_elems": n_elems, "input_format": "fort.14"}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
