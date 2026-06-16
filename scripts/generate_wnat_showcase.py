"""Regenerate ``docs/gallery/wnat_hagen_showcase.png``.

Top panel: ``plot_quality()`` colormap on the full WNAT_Hagen mesh.
Bottom panel: ``plot_quality_histogram()`` (100 bins, cool_r-coded
bars matching the top panel's colormap).

WNAT_Hagen is not bundled with the wheel. Configure its path via:

- ``--mesh /path/to/WNAT_Hagen.14``, or
- environment variable ``WNAT_HAGEN_PATH``, or
- Valence catalog at ``/tmp/valence-domains/registry_data/meshes/WNAT_Hagen.14``.

If none of those resolve, falls back to the largest bundled fixture
(``block_o``) so the script always produces an image; caption in the
README should be updated to reflect the actual mesh used.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib.pyplot as plt

from chilmesh import CHILmesh, examples


def _candidate_paths() -> list[Path]:
    out: list[Path] = []
    env = os.environ.get("WNAT_HAGEN_PATH")
    if env:
        out.append(Path(env))
    out.append(Path("/tmp/valence-domains/registry_data/meshes/WNAT_Hagen.14"))
    out.append(Path.home() / "valence-domains/registry_data/meshes/WNAT_Hagen.14")
    return out


def _resolve_mesh(cli_path: str | None) -> CHILmesh:
    if cli_path:
        return CHILmesh.read_from_fort14(Path(cli_path))
    for candidate in _candidate_paths():
        if candidate.is_file():
            print(f"[mesh] using WNAT_Hagen at {candidate}")
            return CHILmesh.read_from_fort14(candidate)
    print("[mesh] WNAT_Hagen not found; falling back to block_o fixture")
    return examples.block_o()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", help="Path to fort.14 file (overrides discovery)")
    parser.add_argument(
        "--output", default="docs/gallery/wnat_hagen_showcase.png",
        help="Output image path (default: docs/gallery/wnat_hagen_showcase.png)",
    )
    parser.add_argument("--bins", type=int, default=100,
                        help="Histogram bin count (default: 100)")
    parser.add_argument("--dpi", type=int, default=120,
                        help="Figure DPI (default: 120)")
    args = parser.parse_args()

    mesh = _resolve_mesh(args.mesh)
    print(f"[mesh] {mesh.n_verts:,} vertices · {mesh.n_elems:,} elements")

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(12, 12),
        gridspec_kw={"height_ratios": [3, 1]},
    )

    mesh.plot_quality(ax=ax_top)
    # auto_norm=True compresses the cmap to the actual quality range so
    # bar colours separate visually even on narrow distributions.
    mesh.plot_quality_histogram(bins=args.bins, auto_norm=True, ax=ax_bot)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=args.dpi)
    print(f"[wrote] {out_path}")


if __name__ == "__main__":
    main()
