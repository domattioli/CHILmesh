# Size-weighted smoothing benchmark (#197)

Base vs size-weighted `direct_smoother` (`fem`) and `angle_based_smoother` (`angle-based`). `variant=base` is the unchanged size-blind algorithm (`size_fn=None`); `variant=size` passes the target field as `size_fn`; `variant=input` is the mesh before smoothing, for reference.

**Method.** Meshes: the five bundled fixtures plus the Valence meshes found on the benchmark machine. Size field `own` is the input mesh's own local size (per-vertex mean incident edge length, linearly interpolated over the vertex Delaunay triangulation), i.e. the original sizes the smoother should keep; `graded-x` is an analytic field growing linearly with x from h0 to 2.5 h0 (annulus only). `size err` is `|L - h(mid)| / h(mid)` over unique edges at the smoothed edge midpoints, with `h` the target field (mean, 95th percentile, max). `AR` is `chilmesh.element_quality` aspect ratio (1 = equilateral); `% AR<0.5` is the share of elements below 0.5. `min/max angle` are in degrees over all elements; `skew` is `equiangle_skewness` (0 = ideal, 1 = worst; a concave quad scores 1.0 by the documented rule, #275). `inverted` counts elements with non-positive signed area. `time s` is the median of 3 runs of the smoother call alone. Smoother defaults are used for every other parameter (`angle-based`: 100 passes). Regenerate with `PYTHONPATH=src python scripts/benchmark_size_smoothing.py`.

| mesh | size field | method | variant | size err mean | size err p95 | size err max | AR mean | AR min | % AR<0.5 | min angle | max angle | skew mean | skew worst | inverted | time s |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| annulus | own | - | input | 0.279 | 0.651 | 0.955 | 0.638 | 0.003 | 30.000 | 0.7 | 172.0 | 0.506 | 0.988 | 0 | - |
| annulus | own | fem | base | 0.203 | 0.520 | 1.161 | 0.772 | 0.010 | 14.483 | 3.7 | 171.9 | 0.365 | 0.938 | 0 | 0.01 |
| annulus | own | fem | size | 0.191 | 0.506 | 0.817 | 0.775 | 0.072 | 13.103 | 5.8 | 155.3 | 0.365 | 0.903 | 0 | 0.02 |
| annulus | own | angle-based | base | 0.257 | 0.620 | 1.004 | 0.714 | 0.094 | 17.931 | 8.2 | 154.7 | 0.423 | 0.863 | 0 | 2.90 |
| annulus | own | angle-based | size | 0.207 | 0.526 | 0.716 | 0.755 | 0.095 | 15.172 | 8.6 | 154.5 | 0.384 | 0.856 | 0 | 1.19 |
| annulus | graded-x | - | input | 0.487 | 0.847 | 1.255 | 0.638 | 0.003 | 30.000 | 0.7 | 172.0 | 0.506 | 0.988 | 0 | - |
| annulus | graded-x | fem | base | 0.479 | 0.836 | 0.957 | 0.772 | 0.010 | 14.483 | 3.7 | 171.9 | 0.365 | 0.938 | 0 | 0.01 |
| annulus | graded-x | fem | size | 0.451 | 0.820 | 0.885 | 0.763 | 0.079 | 15.862 | 6.8 | 155.3 | 0.383 | 0.886 | 0 | 0.02 |
| annulus | graded-x | angle-based | base | 0.519 | 0.846 | 1.358 | 0.714 | 0.094 | 17.931 | 8.2 | 154.7 | 0.423 | 0.863 | 0 | 2.92 |
| annulus | graded-x | angle-based | size | 0.451 | 0.821 | 0.938 | 0.743 | 0.106 | 16.207 | 8.6 | 152.7 | 0.401 | 0.856 | 0 | 2.87 |
| donut | own | - | input | 0.159 | 0.356 | 0.554 | 0.826 | 0.405 | 1.449 | 22.0 | 124.4 | 0.326 | 0.633 | 0 | - |
| donut | own | fem | base | 0.145 | 0.314 | 0.494 | 0.864 | 0.591 | 0.000 | 30.0 | 110.2 | 0.280 | 0.500 | 0 | 0.01 |
| donut | own | fem | size | 0.130 | 0.312 | 0.463 | 0.865 | 0.595 | 0.000 | 28.1 | 107.1 | 0.284 | 0.532 | 0 | 0.02 |
| donut | own | angle-based | base | 0.186 | 0.400 | 0.534 | 0.834 | 0.577 | 0.000 | 31.0 | 111.1 | 0.311 | 0.483 | 0 | 1.18 |
| donut | own | angle-based | size | 0.166 | 0.377 | 0.491 | 0.848 | 0.566 | 0.000 | 28.9 | 110.8 | 0.299 | 0.519 | 0 | 0.43 |
| block_o | own | - | input | 0.042 | 0.115 | 0.353 | 0.977 | 0.592 | 0.000 | 31.8 | 110.0 | 0.110 | 0.470 | 0 | - |
| block_o | own | fem | base | 0.045 | 0.142 | 0.452 | 0.984 | 0.550 | 0.000 | 30.5 | 112.9 | 0.087 | 0.492 | 0 | 0.10 |
| block_o | own | fem | size | 0.039 | 0.120 | 0.264 | 0.983 | 0.583 | 0.000 | 33.8 | 110.7 | 0.089 | 0.437 | 0 | 0.18 |
| block_o | own | angle-based | base | 0.052 | 0.157 | 0.361 | 0.979 | 0.571 | 0.000 | 33.8 | 111.7 | 0.104 | 0.436 | 0 | 33.15 |
| block_o | own | angle-based | size | 0.044 | 0.125 | 0.262 | 0.980 | 0.568 | 0.000 | 34.0 | 111.9 | 0.102 | 0.433 | 0 | 35.26 |
| structured | own | - | input | 0.235 | 0.345 | 0.351 | 0.747 | 0.747 | 0.000 | 31.2 | 90.0 | 0.480 | 0.480 | 0 | - |
| structured | own | fem | base | 0.235 | 0.345 | 0.351 | 0.747 | 0.747 | 0.000 | 31.2 | 90.0 | 0.480 | 0.480 | 0 | 0.01 |
| structured | own | fem | size | 0.235 | 0.345 | 0.353 | 0.747 | 0.742 | 0.000 | 31.1 | 90.8 | 0.480 | 0.482 | 0 | 0.02 |
| structured | own | angle-based | base | 0.235 | 0.345 | 0.351 | 0.747 | 0.747 | 0.000 | 31.2 | 90.0 | 0.480 | 0.480 | 0 | 0.02 |
| structured | own | angle-based | size | 0.235 | 0.345 | 0.351 | 0.747 | 0.747 | 0.000 | 31.2 | 90.0 | 0.480 | 0.480 | 0 | 0.03 |
| quad_2x2 | own | - | input | 0.000 | 0.000 | 0.000 | 0.828 | 0.828 | 0.000 | 90.0 | 90.0 | 0.000 | 0.000 | 0 | - |
| quad_2x2 | own | fem | base | 0.000 | 0.000 | 0.000 | 0.828 | 0.828 | 0.000 | 90.0 | 90.0 | 0.000 | 0.000 | 0 | 0.00 |
| quad_2x2 | own | fem | size | 0.000 | 0.000 | 0.000 | 0.828 | 0.828 | 0.000 | 90.0 | 90.0 | 0.000 | 0.000 | 0 | 0.00 |
| quad_2x2 | own | angle-based | base | 0.000 | 0.000 | 0.000 | 0.828 | 0.828 | 0.000 | 90.0 | 90.0 | 0.000 | 0.000 | 0 | 0.00 |
| quad_2x2 | own | angle-based | size | 0.000 | 0.000 | 0.000 | 0.828 | 0.828 | 0.000 | 90.0 | 90.0 | 0.000 | 0.000 | 0 | 0.00 |
| Lake_Erie_mesh_refined | own | - | input | 0.117 | 0.260 | 0.471 | 0.910 | 0.500 | 0.041 | 30.8 | 117.2 | 0.235 | 0.486 | 0 | - |
| Lake_Erie_mesh_refined | own | fem | base | 0.118 | 0.259 | 0.583 | 0.932 | 0.488 | 0.021 | 19.8 | 114.7 | 0.209 | 0.670 | 0 | 0.18 |
| Lake_Erie_mesh_refined | own | fem | size | 0.113 | 0.235 | 0.509 | 0.931 | 0.539 | 0.000 | 23.9 | 112.7 | 0.210 | 0.602 | 0 | 0.32 |
| Lake_Erie_mesh_refined | own | angle-based | base | 0.124 | 0.264 | 0.648 | 0.923 | 0.581 | 0.000 | 30.8 | 110.7 | 0.220 | 0.486 | 0 | 64.16 |
| Lake_Erie_mesh_refined | own | angle-based | size | 0.115 | 0.237 | 0.456 | 0.925 | 0.575 | 0.000 | 30.8 | 111.4 | 0.215 | 0.486 | 0 | 59.31 |
| Baranja_Hill | own | - | input | 0.049 | 0.126 | 0.591 | 0.976 | 0.420 | 0.084 | 26.3 | 123.3 | 0.105 | 0.562 | 0 | - |
| Baranja_Hill | own | fem | base | 0.049 | 0.143 | 0.492 | 0.986 | 0.705 | 0.000 | 36.1 | 101.0 | 0.078 | 0.398 | 0 | 0.03 |
| Baranja_Hill | own | fem | size | 0.044 | 0.126 | 0.492 | 0.985 | 0.652 | 0.000 | 34.1 | 105.0 | 0.080 | 0.432 | 0 | 0.05 |
| Baranja_Hill | own | angle-based | base | 0.061 | 0.175 | 0.591 | 0.980 | 0.710 | 0.000 | 36.5 | 100.1 | 0.101 | 0.391 | 0 | 8.70 |
| Baranja_Hill | own | angle-based | size | 0.049 | 0.134 | 0.591 | 0.982 | 0.703 | 0.000 | 37.1 | 101.0 | 0.093 | 0.381 | 0 | 7.01 |
| Test_Case_2 | own | - | input | 0.044 | 0.121 | 0.527 | 0.975 | 0.280 | 0.029 | 21.4 | 135.1 | 0.116 | 0.644 | 0 | - |
| Test_Case_2 | own | fem | base | 0.048 | 0.141 | 0.527 | 0.983 | 0.614 | 0.000 | 30.1 | 106.7 | 0.100 | 0.498 | 0 | 0.07 |
| Test_Case_2 | own | fem | size | 0.041 | 0.121 | 0.527 | 0.983 | 0.453 | 0.029 | 25.8 | 120.2 | 0.101 | 0.570 | 0 | 0.12 |
| Test_Case_2 | own | angle-based | base | 0.059 | 0.167 | 0.527 | 0.976 | 0.583 | 0.000 | 33.6 | 110.7 | 0.119 | 0.439 | 0 | 23.64 |
| Test_Case_2 | own | angle-based | size | 0.048 | 0.137 | 0.527 | 0.979 | 0.662 | 0.000 | 31.8 | 102.7 | 0.113 | 0.470 | 0 | 16.71 |
