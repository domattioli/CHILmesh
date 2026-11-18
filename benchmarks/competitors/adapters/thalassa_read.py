"""Standalone adapter: open a mesh with Thalassa and print one JSON line.

Run ONLY in a subprocess inside the Thalassa environment (see envs/). The call
shape and the input format Thalassa documents as its default path are unverified
until research.md (T002/T005); failure is recorded as ``failed``.
"""
import argparse
import json
import importlib.metadata
import sys
import time


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--fair", action="store_true", help="emit materialized fair-campaign evidence")
    args = ap.parse_args()

    import thalassa
    from thalassa import api
    import xarray  # noqa: F401  backend imports land before READY
    import netCDF4  # noqa: F401
    import thalassa.api  # noqa: F401

    if args.fair:
        from benchmarks.competitors.protocol import run_adapter, eager_dataset, thalassa_output, warm_netcdf_reader
        warm_netcdf_reader(netCDF4, api.open_dataset, element_variable="element", vertex_dimension="nvertex")
        return run_adapter("thalassa", lambda: api.open_dataset(args.mesh), materialize=eager_dataset,
                           extract=thalassa_output)

    t0 = time.perf_counter()
    ds = api.open_dataset(args.mesh)
    ds.load()
    parse_s = time.perf_counter() - t0
    dims = dict(getattr(ds, "sizes", {}) or {})
    print(json.dumps({"tool": "thalassa", "version": importlib.metadata.version("thalassa"),
                      "parse_s": parse_s, "n_nodes": dims.get("node"), "n_elems": dims.get("face", dims.get("nele")),
                      "input_format": "as accepted by thalassa.api.open_dataset"}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
