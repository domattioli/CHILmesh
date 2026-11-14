# Plan 008: v1.7.0 release prep

Spec: `spec.md`. Branch `release/v1.7.0`.

## Approach

- #309 (FR-1..3): `_read_fort14_open_segments` reads the second token of the NVDLL line into `ibtype` (int, or `None`). `_write_fort14_boundaries` writes `NVDLL IBTYPEE` only when `ibtype` is not `None`. Flow boundaries unchanged.
- #306 (FR-4..5): module-level `_quad_corners_ok(px, py)` in `CHILmesh.py` runs the all-corner scaled-Jacobian test (`_SJ_MIN = 1e-6`). `_local_min_quality_fast` in `angle_based_smoother` calls it for quads after the existing diagonal checks and returns -1 on failure. Triangles untouched.
- C++ backend (FR-6): `_local_min_quality_fast` is a pure-Python closure with no compiled mirror, so there is no parity work.
- FR-7: CHANGELOG entry. No version bump.

## Tests

- `tests/test_fort14_roundtrip.py`: IBTYPEE 7 and absent, both round-trip.
- `tests/test_angle_based_smoother.py`: dart quad rejected; jittered 7x7 grid (seeds 4, 14, 41) gains no inverted quads.

## Out of scope

`.chil`, #301 harness, #307 decision record, version bump, tagging, licence files.
