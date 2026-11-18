"""Standalone adapter: open a NetCDF mesh/output with xarray (netCDF4 engine) and print one JSON line.

xarray cannot read fort.14. Pass a NetCDF file (for example ADCIRC ``fort.63.nc``).
Whether the W2 raw-xarray row reads fort.14 or NetCDF is spec open question 7;
for a fort.14 input the runner records ``not-run`` instead of calling this script.
"""
import argparse
import json
import importlib.metadata
import sys
import time


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True, help="path to a NetCDF file")
    ap.add_argument("--fair", action="store_true", help="emit materialized fair-campaign evidence")
    args = ap.parse_args()

    import xarray as xr
    import netCDF4  # noqa: F401  backend imported before READY, outside the timed region

    if args.fair:
        from benchmarks.competitors.protocol import run_adapter, eager_dataset, dataset_output, warm_netcdf_reader
        warm_netcdf_reader(netCDF4, lambda path: xr.open_dataset(path, engine="netcdf4"))
        return run_adapter("xarray", lambda: xr.open_dataset(args.mesh, engine="netcdf4"), materialize=eager_dataset,
                           extract=lambda ds: dataset_output(ds, "nele"))

    t0 = time.perf_counter()
    ds = xr.open_dataset(args.mesh, engine="netcdf4")
    ds.load()  # force the read; open_dataset alone is lazy
    parse_s = time.perf_counter() - t0
    print(json.dumps({"tool": "xarray", "version": importlib.metadata.version("xarray"), "parse_s": parse_s,
                      "n_nodes": int(ds.sizes.get("node", 0)) or None,
                      "n_elems": int(ds.sizes.get("nele", 0)) or None, "input_format": "NetCDF"}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
