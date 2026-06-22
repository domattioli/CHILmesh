"""Command-line interface for chilmesh.

Wraps the public API in a small argparse entry point so non-Python consumers
can inspect, convert, smooth, and plot meshes without writing scripts. See
``chilmesh --help`` for the full subcommand list.

The entry point is registered in ``pyproject.toml`` as
``chilmesh = "chilmesh.cli:main"`` and is also runnable via
``python -m chilmesh``.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence

from . import __version__
from .CHILmesh import CHILmesh
from .summary_io import summary

_READERS = {
    "fort14": "read_from_fort14",
    "2dm": "read_from_2dm",
}

_WRITERS = {
    "fort14": "write_to_fort14",
}

_SMOOTH_METHODS = ("angle-based", "fem", "laplacian")


def _detect_format(path: Path) -> str:
    """Infer mesh format from a filename suffix.

    ``.14`` and ``.fort.14`` → ``fort14``; ``.2dm`` → ``2dm``. Raises
    :class:`SystemExit` with a descriptive message on unknown suffixes so the
    caller sees a clean CLI error rather than a stack trace.
    """
    name = path.name.lower()
    if name.endswith(".2dm"):
        return "2dm"
    if name.endswith(".fort.14") or name.endswith(".14") or name.endswith(".grd"):
        return "fort14"
    raise SystemExit(
        f"error: cannot infer mesh format from '{path.name}'. "
        "Supported suffixes: .14, .fort.14, .grd, .2dm"
    )


def _load_mesh(path: Path, compute_layers: bool = True) -> CHILmesh:
    fmt = _detect_format(path)
    reader_name = _READERS.get(fmt)
    if reader_name is None:
        raise SystemExit(f"error: no reader available for format '{fmt}'")
    reader = getattr(CHILmesh, reader_name)
    return reader(path, compute_layers=compute_layers)


def _write_mesh(mesh: CHILmesh, path: Path) -> None:
    fmt = _detect_format(path)
    writer_name = _WRITERS.get(fmt)
    if writer_name is None:
        raise SystemExit(
            f"error: writing format '{fmt}' is not yet supported. "
            "Supported output formats: fort.14 (.14, .fort.14)"
        )
    writer = getattr(mesh, writer_name)
    writer(str(path), grid_name=mesh.grid_name)


# ----------------------------------------------------------------------
# Subcommand implementations
# ----------------------------------------------------------------------


def cmd_info(args: argparse.Namespace) -> int:
    """Print mesh statistics: vertex/element/edge/layer counts + quality stats."""
    mesh = _load_mesh(Path(args.input), compute_layers=not args.no_layers)

    print(f"Mesh: {mesh.grid_name}")
    print(f"  Source:        {args.input}")
    print(f"  Vertices:      {mesh.n_verts}")
    print(f"  Elements:      {mesh.n_elems}")
    if args.no_layers:
        # `compute_layers=False` also skips adjacency building, so n_edges
        # and n_layers are not populated. Be explicit rather than report 0.
        print("  Edges:         (skipped, --no-layers)")
        print("  Layers:        (skipped, --no-layers)")
    else:
        print(f"  Edges:         {mesh.n_edges}")
        print(f"  Layers:        {mesh.n_layers}")

    if not args.no_quality:
        _, _, stats = mesh.elem_quality()
        print("  Quality (skew):")
        print(f"    min:    {stats['min']:.4f}")
        print(f"    median: {stats['median']:.4f}")
        print(f"    mean:   {stats['mean']:.4f}")
        print(f"    max:    {stats['max']:.4f}")
        print(f"    std:    {stats['std']:.4f}")
    return 0


def cmd_convert(args: argparse.Namespace) -> int:
    """Convert a mesh between supported formats.

    Currently reads fort.14 / .2dm and writes fort.14. The output format is
    inferred from the output suffix.
    """
    in_path = Path(args.input)
    out_path = Path(args.output)
    mesh = _load_mesh(in_path, compute_layers=False)
    _write_mesh(mesh, out_path)
    print(f"Wrote {out_path} ({mesh.n_verts} verts, {mesh.n_elems} elems)")
    return 0


def cmd_smooth(args: argparse.Namespace) -> int:
    """Apply a smoothing algorithm and write the result to a new file."""
    method = args.method
    if method not in _SMOOTH_METHODS:
        # argparse `choices` already covers this, but keep the explicit guard
        # for callers that bypass argparse.
        raise SystemExit(
            f"error: unknown smoothing method '{method}'. "
            f"Choose one of: {', '.join(_SMOOTH_METHODS)}"
        )

    mesh = _load_mesh(Path(args.input), compute_layers=True)
    _, _, pre = mesh.elem_quality()

    if method == "laplacian":
        # The public dispatcher uses 'fem' for the FEM/Laplacian path and
        # 'angle-based' for the angle-based smoother. Map laplacian → fem so
        # the CLI surface matches what the issue spec promises.
        method_internal = "fem"
    else:
        method_internal = method

    for _ in range(max(1, args.iter)):
        mesh.smooth_mesh(method=method_internal, acknowledge_change=True)

    _, _, post = mesh.elem_quality()
    _write_mesh(mesh, Path(args.output))
    print(
        f"Smoothed {args.input} with method={method}, iter={args.iter}. "
        f"median quality {pre['median']:.4f} -> {post['median']:.4f}; "
        f"min {pre['min']:.4f} -> {post['min']:.4f}"
    )
    print(f"Wrote {args.output}")
    return 0


def cmd_plot(args: argparse.Namespace) -> int:
    """Render the mesh to a static figure (PNG/PDF/SVG inferred from suffix)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Default plot needs Elem2Edge adjacency, which is built alongside layers
    # under the current compute_layers flag. Always build them for plotting.
    mesh = _load_mesh(Path(args.input), compute_layers=True)

    if args.quality:
        fig, _ = mesh.plot_quality()
    elif args.layers:
        fig, _ = mesh.plot_layer()
    else:
        fig, _ = mesh.plot()

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")
    return 0


def cmd_summary(args: argparse.Namespace) -> int:
    """Print lazy header-only mesh metadata (counts, no layers/quality)."""
    data = summary(Path(args.input), deep=args.deep)

    # Pretty-print with human-friendly format
    grid_name = data.get("grid_name", "Mesh")
    print(f"Summary: {grid_name}")
    print(f"  Format:     {data.get('format', 'N/A')}")

    file_bytes = data.get("file_bytes")
    if file_bytes is not None:
        print(f"  File bytes: {file_bytes}")
    else:
        print(f"  File bytes: N/A")

    print(f"  Nodes:      {data.get('n_nodes', 'N/A')}")
    print(f"  Elements:   {data.get('n_elems', 'N/A')}")

    # Optional fields (present when deep=True or for mesh objects)
    if "element_type" in data:
        print(f"  Type:       {data['element_type']}")

    if "bbox" in data:
        bbox = data["bbox"]
        min_x = bbox.get("min_x", "N/A")
        max_x = bbox.get("max_x", "N/A")
        min_y = bbox.get("min_y", "N/A")
        max_y = bbox.get("max_y", "N/A")
        print(f"  Bbox:       x[{min_x}..{max_x}] y[{min_y}..{max_y}]")

    return 0


# ----------------------------------------------------------------------
# Parser
# ----------------------------------------------------------------------


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="chilmesh",
        description=(
            "Shell-friendly mesh inspection, conversion, smoothing and "
            "plotting for chilmesh."
        ),
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"chilmesh {__version__}",
    )
    sub = parser.add_subparsers(dest="command", metavar="COMMAND", required=True)

    p_info = sub.add_parser(
        "info",
        help="Print mesh statistics (verts, elems, edges, layers, quality).",
        description=(
            "Load a mesh and print summary statistics.\n\n"
            "Example:\n"
            "  chilmesh info path/to/mesh.fort.14"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_info.add_argument("input", help="Input mesh path (.14, .fort.14, .2dm).")
    p_info.add_argument(
        "--no-layers",
        action="store_true",
        help="Skip skeletonization (faster, no layer count).",
    )
    p_info.add_argument(
        "--no-quality",
        action="store_true",
        help="Skip quality statistics (faster).",
    )
    p_info.set_defaults(func=cmd_info)

    p_summary = sub.add_parser(
        "summary",
        help="Print lazy header-only mesh metadata (counts, no layers/quality).",
        description=(
            "Read only the file header to report counts without loading the full "
            "mesh. Use --deep to also load points for bbox + element type.\n\n"
            "Example:\n  chilmesh summary path/to/mesh.fort.14"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_summary.add_argument("input", help="Input mesh path (.14, .fort.14, .grd, .2dm).")
    p_summary.add_argument(
        "--deep",
        action="store_true",
        help="Also load points for bbox + element type.",
    )
    p_summary.set_defaults(func=cmd_summary)

    p_conv = sub.add_parser(
        "convert",
        help="Convert a mesh between supported formats.",
        description=(
            "Convert a mesh between formats. Output format is inferred from "
            "the output suffix.\n\n"
            "Example:\n"
            "  chilmesh convert mesh.2dm mesh.fort.14"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_conv.add_argument("input", help="Input mesh path (.14, .fort.14, .2dm).")
    p_conv.add_argument("output", help="Output mesh path (.14, .fort.14).")
    p_conv.set_defaults(func=cmd_convert)

    p_smooth = sub.add_parser(
        "smooth",
        help="Smooth a mesh and write the result.",
        description=(
            "Apply mesh smoothing in place and write to a new file.\n\n"
            "Example:\n"
            "  chilmesh smooth mesh.fort.14 -o smoothed.fort.14 "
            "--method angle-based --iter 50"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_smooth.add_argument("input", help="Input mesh path.")
    p_smooth.add_argument(
        "-o", "--output", required=True, help="Output mesh path."
    )
    p_smooth.add_argument(
        "--method",
        choices=_SMOOTH_METHODS,
        default="angle-based",
        help="Smoothing algorithm to apply (default: angle-based).",
    )
    p_smooth.add_argument(
        "--iter",
        type=int,
        default=1,
        help="Number of full smoothing passes (default: 1).",
    )
    p_smooth.set_defaults(func=cmd_smooth)

    p_plot = sub.add_parser(
        "plot",
        help="Render a static figure of the mesh.",
        description=(
            "Render the mesh to a figure (PNG/PDF/SVG by suffix).\n\n"
            "Example:\n"
            "  chilmesh plot mesh.fort.14 -o mesh.png --quality"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_plot.add_argument("input", help="Input mesh path.")
    p_plot.add_argument(
        "-o", "--output", required=True, help="Output figure path."
    )
    p_plot.add_argument(
        "--layers",
        action="store_true",
        help="Color elements by skeletonization layer.",
    )
    p_plot.add_argument(
        "--quality",
        action="store_true",
        help="Color elements by skew quality (overrides --layers).",
    )
    p_plot.add_argument(
        "--dpi", type=int, default=150, help="Figure resolution (default: 150)."
    )
    p_plot.set_defaults(func=cmd_plot)

    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Entry point used both by the installed script and ``python -m chilmesh``."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except SystemExit:
        raise
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except Exception as exc:  # pragma: no cover - defensive top-level guard
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
