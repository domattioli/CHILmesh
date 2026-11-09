"""Command line: ``python -m benchmarks.competitors w2|w1|report``. Real runs need CHILMESH_RUN_BENCH=1."""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

COMMANDS = ("w2", "w1", "report")
GATE_MESSAGE = ("refusing to run: real benchmark runs download meshes and take a long time. "
                "Set CHILMESH_RUN_BENCH=1 to opt in.")

W1_METRICS = ("q_min_after", "q_mean_minus_3sigma_after", "pct_angles_lt_10_after", "max_courant_after",
              "dt_max_after", "boundary_displacement", "wall_s")


def _parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="python -m benchmarks.competitors")
    sub = ap.add_subparsers(dest="command", required=True)
    for name in ("w2", "w1"):
        p = sub.add_parser(name, help=f"run workflow {name} (needs CHILMESH_RUN_BENCH=1)")
        p.add_argument("--mesh", required=True, action="append", help="Valence mesh name; repeatable")
        p.add_argument("--out", required=True, help="directory for run records")
        p.add_argument("--timeout", type=float, default=3600.0)
        if name == "w2":
            p.add_argument("--repeats", type=int, default=3)
            p.add_argument("--cache", choices=["cold", "warm", "unknown"], default="unknown")
            p.add_argument("--aux-dir", help="directory with the fetched CO-OPS station list and "
                                             "points.cwl.nc; absent files record station sampling as not run")
    r = sub.add_parser("report", help="render Markdown from run records")
    r.add_argument("--records", required=True)
    r.add_argument("--workflow", choices=["w1", "w2"], required=True)
    r.add_argument("--published", action="store_true", help="exclude non-cold runs")
    r.add_argument("--rights", help="rights.json; meshes without a publishable entry are excluded")
    return ap


def _report(args) -> int:
    from . import records, valence
    from .report import render
    from .w2 import tools as w2tools

    directory = Path(args.records)
    recs = [r for r in (records.read_records(directory) if directory.is_dir() else []) if r["workflow"] == args.workflow]
    if not recs:
        print(f"no records for {args.workflow} in {directory}")
        return 0
    exclusions = []
    if args.rights:
        names = sorted({r["mesh"]["name"] for r in recs})
        allowed, exclusions = valence.apply_rights(names, valence.load_rights(args.rights))
        recs = [r for r in recs if r["mesh"]["name"] in allowed]
    if args.workflow == "w2":
        tools, metrics = w2tools.all_tools(), list(w2tools.METRICS)
    else:
        tools = sorted({r["tool"] for r in recs})
        metrics = list(W1_METRICS)
    print(render.render_sections(recs, tools, metrics, published=args.published, exclusions=exclusions))
    return 0


def main(argv=None) -> int:
    args = _parser().parse_args(argv)
    if args.command == "report":
        return _report(args)
    if os.environ.get("CHILMESH_RUN_BENCH") != "1":
        print(GATE_MESSAGE, file=sys.stderr)
        return 2
    if args.command == "w2":
        from .w2 import run
        all_recs = []
        for name in args.mesh:
            all_recs += run.run_mesh(name, args.out, repeats=args.repeats, timeout_s=args.timeout,
                                     cache=args.cache, aux_dir=args.aux_dir)
        scaling = run.collect_scaling(all_recs, args.out)
        print(f"scaling: sizes={scaling['sizes']} summary={scaling['summary'].get('best') or scaling['summary']['note']}"
              f" not_run={len(scaling['not_run'])}")
    else:
        from .w1 import run
        for name in args.mesh:
            run.run_mesh(name, args.out, timeout_s=args.timeout)
    return 0


if __name__ == "__main__":
    sys.exit(main())
