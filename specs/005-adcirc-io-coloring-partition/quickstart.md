<!-- provenance: author=domattioli model=claude-opus-5-5 effort=medium date=2026-10-12 skill=speckit-pipeline repo=CHILmesh session=4e33ac09-e474-4ebc-a2dd-027a9e0a0442 -->
# Quickstart: #289 B1-B3

```bash
pip install -e .[dev]            # add [netcdf] for NetCDF readers
pytest tests/test_2dm_nodestrings.py tests/test_coloring.py tests/test_partition_report.py   # B1
pytest tests/test_adcirc_output_ascii.py tests/test_adcirc_output_netcdf.py tests/test_plot_field.py  # B2
pytest tests/test_fort22.py tests/test_owi.py tests/test_boundary_forcing.py                   # B3
python scripts/check_adcirc_testsuite.py --stage all   # local real-file gate, needs network, not CI
```
