"""Command line: ``python -m benchmarks.competitors w2|w1|report``. Real runs need CHILMESH_RUN_BENCH=1."""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

COMMANDS = ("w2", "w1", "report", "fair-preflight", "fair-run", "fair-report")
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
        p.add_argument("--cache", choices=["cold", "warm", "unknown"], default="unknown")
        if name == "w2":
            p.add_argument("--repeats", type=int, default=3)
            p.add_argument("--aux-dir", help="directory with the fetched CO-OPS station list and "
                                             "points.cwl.nc; absent files record station sampling as not run")
    r = sub.add_parser("report", help="render Markdown from run records")
    r.add_argument("--records", required=True)
    r.add_argument("--workflow", choices=["w1", "w2"], required=True)
    r.add_argument("--published", action="store_true", help="exclude non-cold runs")
    r.add_argument("--rights", help="rights.json; meshes without a publishable entry are excluded")
    for name in ('fair-preflight', 'fair-run'):
        p = sub.add_parser(name, help='validate a frozen fair campaign' if name == 'fair-preflight' else 'execute a prepared fair campaign (opt-in)')
        p.add_argument('--manifest', required=True)
        p.add_argument('--out', required=True)
    p = sub.add_parser('fair-report', help='recompute saved fair evidence offline')
    p.add_argument('--campaign', required=True)
    p.add_argument('--published', action='store_true')
    p.add_argument('--rights')
    p.add_argument('--claims')
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


def _fair(args):
    import json
    from . import campaign, records, valence
    if args.command == 'fair-run' and os.environ.get('CHILMESH_RUN_BENCH') != '1':
        print(GATE_MESSAGE, file=sys.stderr)
        return 2
    try:
        if args.command == 'fair-report':
            from .report.render import fair_report
            rights = valence.load_rights(args.rights) if args.rights else {}
            proposals = json.loads(Path(args.claims).read_text()) if args.claims else []
            if not isinstance(proposals, list):
                raise ValueError('claims must be a list of structured proposals')
            result = fair_report(args.campaign, published=args.published, rights=rights, proposals=proposals)
            print(result['markdown'])
            return 0 if result['publication_passed'] else 1
        manifest = campaign.validate_manifest(json.loads(Path(args.manifest).read_text()))
        preflight = campaign.execution_preflight(manifest)
        out = Path(args.out)
        path = out / 'evidence' / 'preflight.json'
        # Preserve earlier preflight evidence, never replace it silently.
        if not path.exists():
            records._write_new((json.dumps(preflight, sort_keys=True, allow_nan=False)+'\n').encode(), path)
        elif json.loads(path.read_text()) != preflight:
            raise ValueError('preflight changed; use a new output directory/campaign')
        if args.command == 'fair-preflight':
            print(json.dumps(preflight, sort_keys=True))
            return 0 if preflight['ready'] else 1
        prepared_path = Path(args.manifest).parent / 'evidence' / 'prepared.json'
        reason = '; '.join(preflight['reasons'])
        if not prepared_path.exists():
            reason = (reason+'; missing evidence/prepared.json').lstrip('; ')
        if reason:
            campaign.finalize_unstarted_campaign(manifest, out, 'preflight: '+reason)
            print('execution blocked: '+reason, file=sys.stderr)
            return 1
        prepared = json.loads(prepared_path.read_text())
        from .w2.run import run_fair_campaign
        try:
            rows = run_fair_campaign(manifest, prepared['identities'], prepared['references'], prepared['commands'], out,
                                     conversions=prepared.get('conversions'))
        except (campaign.CampaignError, KeyError, OSError, ValueError) as exc:
            if not (out/'evidence'/'campaign.json').exists():
                campaign.finalize_unstarted_campaign(manifest,out,'preparation/preflight: '+str(exc))
            raise
        return 1 if any('controller abort' in r['reason'] for r in rows) else 0
    except (OSError, ValueError, TypeError, KeyError, records.MissingCellError) as exc:
        # Public integrity errors cannot echo private record values or diagnostics.
        print('fair command failed: '+ (type(exc).__name__ if getattr(args,'published',False) else str(exc)), file=sys.stderr)
        return 1


def main(argv=None) -> int:
    args = _parser().parse_args(argv)
    if args.command.startswith('fair-'):
        return _fair(args)
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
            run.run_mesh(name, args.out, timeout_s=args.timeout, cache=args.cache)
    return 0


if __name__ == "__main__":
    sys.exit(main())
