#!/usr/bin/env python
"""Benchmark size-weighted vs base smoothing (#197, option A).

Compares ``direct_smoother`` (FEM) and ``angle_based_smoother`` with and without
the optional ``size_fn`` on the five bundled fixtures plus the Valence corpus
meshes that exist on this machine, and writes ``docs/benchmarks/size_smoothing.md``.

Run from the repository root::

    PYTHONPATH=src python scripts/benchmark_size_smoothing.py

Size fields
-----------
``own``     the input mesh's own local size: per-vertex mean incident edge
            length, linearly interpolated over the vertex Delaunay
            triangulation (nearest vertex outside the hull). This is the
            "keep the original truss-solved sizes" case.
``graded``  analytic field growing linearly with x (h0 at the left edge,
            2.5 h0 at the right edge), annulus only.

Metrics are computed on the smoothed points of the unchanged connectivity.
The size error is ``|L_e - h(mid_e)| / h(mid_e)`` over unique edges, with
``h`` the target field evaluated at the smoothed edge midpoints.
"""
from __future__ import annotations

import argparse
import statistics
import sys
import time
import warnings
from pathlib import Path

import numpy as np
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

import chilmesh
from chilmesh import CHILmesh, element_quality
from chilmesh.CHILmesh import _shoelace, _unique_in_order

REPO = Path(__file__).resolve().parents[1]
VALENCE = Path("/Users/domattioli/Projects/Valence/registry_data/meshes")
VALENCE_MESHES = ["Lake_Erie_mesh_refined.14", "Baranja_Hill.14", "Test_Case_2.14"]
FIXTURES = ["annulus", "donut", "block_o", "structured", "quad_2x2"]

COLUMNS = [
    "mesh", "size field", "method", "variant",
    "size err mean", "size err p95", "size err max",
    "AR mean", "AR min", "% AR<0.5",
    "min angle", "max angle", "skew mean", "skew worst",
    "inverted", "time s",
]


def edge_lengths(mesh: CHILmesh, pts: np.ndarray):
    e = mesh.adjacencies["Edge2Vert"]
    p1, p2 = pts[e[:, 0], :2], pts[e[:, 1], :2]
    return np.linalg.norm(p1 - p2, axis=1), 0.5 * (p1 + p2)


def own_size_fn(mesh: CHILmesh):
    """Per-vertex mean incident edge length, linearly interpolated (nearest outside hull)."""
    xy = mesh.points[:, :2]
    e = mesh.adjacencies["Edge2Vert"]
    length, _ = edge_lengths(mesh, mesh.points)
    total = np.zeros(mesh.n_verts)
    count = np.zeros(mesh.n_verts)
    for col in (0, 1):
        np.add.at(total, e[:, col], length)
        np.add.at(count, e[:, col], 1.0)
    h_vert = total / np.maximum(count, 1.0)
    lin = LinearNDInterpolator(xy, h_vert)
    near = NearestNDInterpolator(xy, h_vert)

    def fn(q: np.ndarray) -> np.ndarray:
        q = np.asarray(q, dtype=float)
        h = lin(q)
        bad = ~np.isfinite(h)
        if bad.any():
            h[bad] = near(q[bad])
        return h

    return fn


def graded_size_fn(mesh: CHILmesh):
    xy = mesh.points[:, :2]
    lo, hi = xy[:, 0].min(), xy[:, 0].max()
    length, _ = edge_lengths(mesh, mesh.points)
    h0 = float(length.mean())
    return lambda q: h0 * (1.0 + 1.5 * (np.asarray(q)[:, 0] - lo) / (hi - lo))


def metrics(mesh: CHILmesh, pts: np.ndarray, size_fn) -> dict:
    length, mid = edge_lengths(mesh, pts)
    target = size_fn(mid)
    err = np.abs(length - target) / target
    xy = pts[:, :2]
    conn = mesh.connectivity_list
    ar = element_quality(xy, conn, "aspect_ratio")
    ang_min = np.degrees(element_quality(xy, conn, "min_angle"))
    ang_max = np.degrees(element_quality(xy, conn, "max_angle"))
    skew = element_quality(xy, conn, "equiangle_skewness")
    inverted = sum(_shoelace(xy[_unique_in_order(r)]) <= 0.0 for r in conn)
    return {
        "size err mean": err.mean(), "size err p95": np.percentile(err, 95),
        "size err max": err.max(),
        "AR mean": ar.mean(), "AR min": ar.min(), "% AR<0.5": 100.0 * np.mean(ar < 0.5),
        "min angle": ang_min.min(), "max angle": ang_max.max(),
        "skew mean": skew.mean(), "skew worst": skew.max(),
        "inverted": int(inverted),
    }


def timed(fn, repeats: int):
    times, out = [], None
    for _ in range(repeats):
        t0 = time.perf_counter()
        out = fn()
        times.append(time.perf_counter() - t0)
    return out, statistics.median(times)


def load_meshes() -> list[tuple[str, CHILmesh]]:
    meshes = [(n, getattr(chilmesh.examples, n)()) for n in FIXTURES]
    for fname in VALENCE_MESHES:
        path = VALENCE / fname
        if path.exists():
            meshes.append((fname.removesuffix(".14"), CHILmesh.read_from_fort14(path)))
        else:
            print(f"skip (missing): {path}", file=sys.stderr)
    return meshes


def run(repeats: int) -> list[dict]:
    rows = []
    for name, mesh in load_meshes():
        fields = [("own", own_size_fn(mesh))]
        if name == "annulus":
            fields.append(("graded-x", graded_size_fn(mesh)))
        for field_name, fn in fields:
            target = fn  # target for the error metric is always this field
            row0 = {"mesh": name, "size field": field_name, "method": "-", "variant": "input",
                    "time s": float("nan")}
            row0.update(metrics(mesh, mesh.points, target))
            rows.append(row0)
            methods = [
                ("fem", lambda **kw: mesh.direct_smoother(**kw)),
                ("angle-based", lambda **kw: mesh.angle_based_smoother(**kw)),
            ]
            for method, call in methods:
                for variant, kw in (("base", {}), ("size", {"size_fn": fn})):
                    print(f"  {name} / {field_name} / {method} / {variant}", flush=True)
                    pts, sec = timed(lambda: call(**kw), repeats)
                    row = {"mesh": name, "size field": field_name, "method": method,
                           "variant": variant, "time s": sec}
                    row.update(metrics(mesh, pts, target))
                    rows.append(row)
    return rows


def fmt(col: str, v) -> str:
    if isinstance(v, str):
        return v
    if col == "inverted":
        return str(int(v))
    if col == "time s":
        return "-" if not np.isfinite(v) else f"{v:.2f}"
    if col in ("min angle", "max angle"):
        return f"{v:.1f}"
    return f"{v:.3f}"


def check_cells(rows: list[dict]) -> None:
    for r in rows:
        for c in COLUMNS[4:]:
            v = r[c]
            if r["variant"] == "input" and c == "time s":
                continue
            if not np.isfinite(v):
                raise SystemExit(f"non-finite cell: {r['mesh']}/{r['method']}/{r['variant']}/{c}")


def write_markdown(rows: list[dict], out: Path, repeats: int) -> None:
    lines = [
        "# Size-weighted smoothing benchmark (#197)",
        "",
        "Base vs size-weighted `direct_smoother` (`fem`) and `angle_based_smoother` "
        "(`angle-based`). `variant=base` is the unchanged size-blind algorithm "
        "(`size_fn=None`); `variant=size` passes the target field as `size_fn`; "
        "`variant=input` is the mesh before smoothing, for reference.",
        "",
        "**Method.** Meshes: the five bundled fixtures plus the Valence meshes found "
        "on the benchmark machine. Size field `own` is the input mesh's own local size "
        "(per-vertex mean incident edge length, linearly interpolated over the vertex "
        "Delaunay triangulation), i.e. the original sizes the smoother should keep; "
        "`graded-x` is an analytic field growing linearly with x from h0 to 2.5 h0 "
        "(annulus only). `size err` is `|L - h(mid)| / h(mid)` over unique edges at the "
        "smoothed edge midpoints, with `h` the target field (mean, 95th percentile, max). "
        "`AR` is `chilmesh.element_quality` aspect ratio (1 = equilateral); `% AR<0.5` is "
        "the share of elements below 0.5. `min/max angle` are in degrees over all "
        "elements; `skew` is `equiangle_skewness` (0 = ideal, 1 = worst; a concave quad "
        "scores 1.0 by the documented rule, #275). `inverted` counts elements with "
        f"non-positive signed area. `time s` is the median of {repeats} runs of the "
        "smoother call alone. Smoother defaults are used for every other parameter "
        "(`angle-based`: 100 passes). Regenerate with "
        "`PYTHONPATH=src python scripts/benchmark_size_smoothing.py`.",
        "",
        "| " + " | ".join(COLUMNS) + " |",
        "|" + "|".join(["---"] * 4 + ["---:"] * (len(COLUMNS) - 4)) + "|",
    ]
    for r in rows:
        lines.append("| " + " | ".join(fmt(c, r[c]) for c in COLUMNS) + " |")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--repeats", type=int, default=3, help="timing repeats (median)")
    ap.add_argument("--out", type=Path, default=REPO / "docs/benchmarks/size_smoothing.md")
    args = ap.parse_args()
    warnings.simplefilter("ignore")
    rows = run(args.repeats)
    check_cells(rows)
    write_markdown(rows, args.out, args.repeats)
    print(args.out.read_text())
    return 0


if __name__ == "__main__":
    sys.exit(main())
