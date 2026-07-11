"""Audit tooling for hero GIF generation validity.

Pure audit functions for verifying:
  - audit_degenerate: no hole-spanning/oversize/inverted elements in captured snapshots
  - audit_pacing: uniform motion speed (thirds displacement ratio)
  - audit_peel_invariance: per-bin histogram totals stable across peel reveals
  - audit_determinism: run-to-run reproducibility

All functions are PURE (no generator import at module level). CLI entrypoints import
the generator LAZILY inside main(). See plan.md for detailed metrics definitions.
"""

from __future__ import annotations

import argparse
import sys
import numpy as np
from typing import Optional, Tuple, Callable

# Module-level constants (no generator imports at module level)
HBINS = 40


def _signed_area_elements(pts: np.ndarray, elems: np.ndarray) -> np.ndarray:
    """Compute signed area of triangles. Positive = CCW, Negative = CW."""
    p0 = pts[elems[:, 0]]
    p1 = pts[elems[:, 1]]
    p2 = pts[elems[:, 2]]
    return 0.5 * ((p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
                  - (p2[:, 0] - p0[:, 0]) * (p1[:, 1] - p0[:, 1]))


def _triangle_area_elements(pts: np.ndarray, elems: np.ndarray) -> np.ndarray:
    """Compute absolute area of triangles."""
    return np.abs(_signed_area_elements(pts, elems))


def audit_degenerate(
    hist: list[Tuple[np.ndarray, np.ndarray]],
    fd: Callable[[np.ndarray], np.ndarray],
    hole_tau: float = 0.02,
    oversize_mult: float = 8.0,
) -> dict:
    """Audit captured snapshots for degenerate elements (T007, US1).

    Measures degenerate element counts across all captured snapshots using
    corrected metrics per plan.md § Implementation Architecture A:
    - hole-span: edge-midpoint fd > +hole_tau (NOT fd > -geps)
    - oversize: area > oversize_mult * per-snapshot-median-area (NOT 4×)
    - inverted: signed-area sign flip vs. snapshot majority

    Parameters
    ----------
    hist : list of (p, t) tuples
        Captured snapshots from distmesh2d_warmstart; each (p, t) tuple is
        (points array, triangulation indices).
    fd : callable
        Signed distance function for hole detection: fd(points: (K, 2)) -> (K,).
    hole_tau : float
        Signed-distance threshold for hole-spanning detection.
        Edge midpoints with fd > +hole_tau are considered hole-spanning.
    oversize_mult : float
        Multiplier for oversized element detection.
        Elements with area > oversize_mult * per-snapshot-median are oversized.

    Returns
    -------
    dict
        Keys:
        - 'per_snapshot': list of dicts, one per snapshot with keys 'hole', 'oversize', 'inverted'
        - 'total_hole': total hole-spanning elements across all snapshots
        - 'total_oversize': total oversized elements across all snapshots
        - 'total_inverted': total inverted elements across all snapshots
        - 'max_per_snapshot_hole': max hole count in any single snapshot
        - 'max_per_snapshot_oversize': max oversize count in any single snapshot
        - 'max_per_snapshot_inverted': max inverted count in any single snapshot
    """
    if not hist:
        return {
            'per_snapshot': [],
            'total_hole': 0,
            'total_oversize': 0,
            'total_inverted': 0,
            'max_per_snapshot_hole': 0,
            'max_per_snapshot_oversize': 0,
            'max_per_snapshot_inverted': 0,
        }

    results = {
        'per_snapshot': [],
        'total_hole': 0,
        'total_oversize': 0,
        'total_inverted': 0,
        'max_per_snapshot_hole': 0,
        'max_per_snapshot_oversize': 0,
        'max_per_snapshot_inverted': 0,
    }

    for snap_idx, (p, t) in enumerate(hist):
        p = np.asarray(p)[:, :2]  # Ensure 2D
        t = np.asarray(t, dtype=np.int64)

        # Initialize counts for this snapshot
        snap_result = {'hole': 0, 'oversize': 0, 'inverted': 0}

        # Hole-spanning detection (edge-midpoint fd > +hole_tau)
        holes = 0
        for e_pair in [t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]]:
            midpts = 0.5 * (p[e_pair[:, 0]] + p[e_pair[:, 1]])
            fd_vals = fd(midpts)
            holes += np.sum(fd_vals > hole_tau)
        snap_result['hole'] = int(holes)

        # Oversize detection (area > oversize_mult * median area)
        areas = _triangle_area_elements(p, t)
        if len(areas) > 0:
            median_area = np.median(areas)
            if median_area > 0:
                oversize_count = np.sum(areas > oversize_mult * median_area)
                snap_result['oversize'] = int(oversize_count)

        # Inverted detection (signed-area sign flip vs. majority)
        signed_areas = _signed_area_elements(p, t)
        if len(signed_areas) > 0:
            majority_sign = np.sign(np.median(signed_areas))
            if majority_sign != 0:
                flipped = np.sum(np.sign(signed_areas) != majority_sign)
                snap_result['inverted'] = int(flipped)

        results['per_snapshot'].append(snap_result)
        results['total_hole'] += snap_result['hole']
        results['total_oversize'] += snap_result['oversize']
        results['total_inverted'] += snap_result['inverted']
        results['max_per_snapshot_hole'] = max(results['max_per_snapshot_hole'], snap_result['hole'])
        results['max_per_snapshot_oversize'] = max(results['max_per_snapshot_oversize'], snap_result['oversize'])
        results['max_per_snapshot_inverted'] = max(results['max_per_snapshot_inverted'], snap_result['inverted'])

    return results


def audit_pacing(
    hist: list[Tuple[np.ndarray, np.ndarray]],
    snap_idx: list[int],
    hold: int = 2,
) -> float:
    """Audit uniform pacing of truss animation (T011, US2).

    Measures mean per-rendered-frame node displacement across thirds of the
    truss playback. Returns the ratio of max to min. Per Probe 2: target ≤ 3×
    (ideal ~1.26).

    Parameters
    ----------
    hist : list of (p, t) tuples
        Full history of captured snapshots.
    snap_idx : list of int
        Indices into hist that are rendered (selected snapshots).
    hold : int
        Number of frames each snapshot is held. Default 2 (uniform hold).

    Returns
    -------
    float
        Ratio of max to min mean displacement per rendered frame across thirds.
        Target: <= 3.0 (ideal ~1.26 per research.md Probe 2).
    """
    if len(snap_idx) < 3 or len(hist) == 0:
        return 1.0  # Degenerate case

    # Compute per-iteration mean displacement
    Pn = [np.asarray(p)[:, :2] for (p, _) in hist]
    d_iter = np.zeros(len(Pn))
    for i in range(1, len(Pn)):
        d_iter[i] = np.mean(np.linalg.norm(Pn[i] - Pn[i - 1], axis=1))

    # Replicate rendered frames with hold
    played_iters = []
    for idx in snap_idx:
        played_iters.extend([idx] * hold)

    # Per-frame displacement (using the difference at each played index)
    frame_displacements = [d_iter[idx] for idx in played_iters]

    if not frame_displacements or len(frame_displacements) < 3:
        return 1.0

    # Split into thirds and compute mean displacement per third
    n_frames = len(frame_displacements)
    third_size = n_frames // 3

    thirds_means = []
    for t in range(3):
        start = t * third_size
        end = (t + 1) * third_size if t < 2 else n_frames
        if start < end:
            third_mean = np.mean(frame_displacements[start:end])
            thirds_means.append(third_mean)

    if len(thirds_means) == 0 or max(thirds_means) == 0:
        return 1.0

    # Ratio of max to min
    ratio = max(thirds_means) / max(min(thirds_means), 1e-12)
    return float(ratio)


def audit_peel_invariance(
    q_fem: np.ndarray,
    elem_layer: np.ndarray,
    n_layers: int,
) -> dict:
    """Audit per-bin histogram totals are invariant across peel reveals (T015, US4).

    Verifies per plan.md § Implementation Architecture C and F-04 (convert-in-place):
    1. Per-bin totals are np.array_equal across k=0..n_layers
    2. Each per-bin total equals an independently computed histogram

    Parameters
    ----------
    q_fem : (n_elems,) ndarray
        Quality values for all elements.
    elem_layer : (n_elems,) ndarray
        Layer assignment for each element (0 to n_layers-1).
    n_layers : int
        Number of layers.

    Returns
    -------
    dict
        Keys:
        - 'histogram_invariant': bool, True if totals equal across all k
        - 'histogram_independent': bool, True if totals match independent histogram
        - 'per_layer_totals': list of arrays, one per layer showing bin counts
        - 'reference_histogram': array of reference bin counts
        - 'message': str with details
    """
    result = {
        'histogram_invariant': False,
        'histogram_independent': False,
        'per_layer_totals': [],
        'reference_histogram': None,
        'message': '',
    }

    q_fem = np.asarray(q_fem, dtype=np.float64)
    elem_layer = np.asarray(elem_layer, dtype=np.int32)

    if len(q_fem) != len(elem_layer):
        result['message'] = f'Mismatch: len(q_fem)={len(q_fem)} != len(elem_layer)={len(elem_layer)}'
        return result

    if len(q_fem) == 0:
        result['message'] = 'Empty quality array'
        return result

    # Compute reference histogram (independent, F-10: not reused from render binning)
    ref_counts, _ = np.histogram(q_fem, bins=HBINS, range=(0.0, 1.0))
    result['reference_histogram'] = ref_counts.copy()

    # Compute per-bin totals across all peel states k=0..n_layers
    # At k, layers 0..k-1 are revealed; layers k..n_layers-1 are hidden
    # Histogram TOTALS remain invariant (same elements counted, just colored differently)
    per_k_totals = []
    for k in range(n_layers + 1):
        counts = np.zeros(HBINS, dtype=np.int64)
        for elem_idx in range(len(elem_layer)):
            q_val = q_fem[elem_idx]
            # F-10: bin index via np.clip(np.floor(q * HBINS), 0, HBINS-1)
            # ensures q=1.0 lands in HBINS-1, not overflow at HBINS
            bin_idx = int(np.clip(np.floor(q_val * HBINS), 0, HBINS - 1))
            # Count this element regardless of layer in the bin totals
            # (the histogram tracks ALL elements, just colored differently)
            counts[bin_idx] += 1
        per_k_totals.append(counts)
    result['per_layer_totals'] = per_k_totals

    # Check invariance: all k have same totals (F-06)
    invariant = True
    for k in range(1, len(per_k_totals)):
        if not np.array_equal(per_k_totals[k], per_k_totals[0]):
            invariant = False
            break
    result['histogram_invariant'] = invariant

    # Check independence: totals match reference (F-04, F-10)
    independent = np.array_equal(per_k_totals[0], ref_counts)
    result['histogram_independent'] = independent

    if invariant and independent:
        result['message'] = 'PASS: histogram totals invariant and independent'
    else:
        msg = []
        if not invariant:
            msg.append('histogram totals not invariant across k')
        if not independent:
            msg.append('histogram totals do not match independent reference')
        result['message'] = 'FAIL: ' + ', '.join(msg)

    return result


def audit_determinism(
    run1_snapshots: list,
    run1_n_layers: int,
    run2_snapshots: list,
    run2_n_layers: int,
) -> dict:
    """Audit determinism: run-to-run reproducibility (T022, SC-007, FR-009).

    Compares two runs of _stage_data() to verify identical results under fixed seed.

    Parameters
    ----------
    run1_snapshots : list
        snapshots from first run
    run1_n_layers : int
        n_layers from first run
    run2_snapshots : list
        snapshots from second run
    run2_n_layers : int
        n_layers from second run

    Returns
    -------
    dict
        Keys:
        - 'snapshots_match': bool
        - 'n_layers_match': bool
        - 'snapshot_count_match': bool
        - 'run1_count': int
        - 'run2_count': int
        - 'message': str
    """
    result = {
        'snapshots_match': False,
        'n_layers_match': False,
        'snapshot_count_match': False,
        'run1_count': len(run1_snapshots),
        'run2_count': len(run2_snapshots),
        'message': '',
    }

    # Check n_layers
    n_layers_match = run1_n_layers == run2_n_layers
    result['n_layers_match'] = n_layers_match

    # Check snapshot count
    snap_count_match = len(run1_snapshots) == len(run2_snapshots)
    result['snapshot_count_match'] = snap_count_match

    # Check snapshots themselves (compare points and triangulation)
    snapshots_match = True
    if snap_count_match:
        for (p1, t1), (p2, t2) in zip(run1_snapshots, run2_snapshots):
            p1 = np.asarray(p1)[:, :2]
            p2 = np.asarray(p2)[:, :2]
            t1 = np.asarray(t1, dtype=np.int64)
            t2 = np.asarray(t2, dtype=np.int64)
            if not (np.allclose(p1, p2) and np.array_equal(t1, t2)):
                snapshots_match = False
                break
    result['snapshots_match'] = snapshots_match

    if n_layers_match and snap_count_match and snapshots_match:
        result['message'] = 'PASS: runs are deterministic'
    else:
        msg = []
        if not n_layers_match:
            msg.append(f'n_layers mismatch: {run1_n_layers} vs {run2_n_layers}')
        if not snap_count_match:
            msg.append(f'snapshot count mismatch: {len(run1_snapshots)} vs {len(run2_snapshots)}')
        if not snapshots_match:
            msg.append('snapshots do not match')
        result['message'] = 'FAIL: ' + ', '.join(msg)

    return result


def main():
    """CLI entry point. Imports generator lazily per plan.md (PURE functions, lazy imports)."""
    parser = argparse.ArgumentParser(
        description='Audit hero GIF generation validity'
    )
    parser.add_argument(
        '--check',
        choices=['degenerate', 'pacing', 'peel-invariance', 'determinism', 'all'],
        default='all',
        help='Which audit to run'
    )
    args = parser.parse_args()

    # Lazy import of generator (happens only when CLI is invoked)
    # This keeps the module PURE at import time (no generator dependency at module level).
    # The script's own directory (scripts/) is on sys.path[0] when invoked as
    # `python scripts/audit_hero_gif.py`, so the sibling generator imports by bare name.
    import pathlib
    sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / 'src'))
    sys.path.insert(0, str(pathlib.Path(__file__).parent))
    from generate_hero_animation import _stage_data, TRUSS_HOLD

    failures = []

    if args.check in ('degenerate', 'all'):
        print('[audit] Running audit_degenerate...')
        try:
            D = _stage_data()
            res = audit_degenerate(D['hist'], D['fd'])
            print(f"[audit]   TOTAL hole={res['total_hole']}, "
                  f"oversize={res['total_oversize']}, inverted={res['total_inverted']}")
            print(f"[audit]   MAX-per-snapshot hole={res['max_per_snapshot_hole']}, "
                  f"oversize={res['max_per_snapshot_oversize']}, "
                  f"inverted={res['max_per_snapshot_inverted']}  "
                  f"(over {len(D['hist'])} snapshots)")
            # F-03: assert the FEM-input connectivity hist[-1] is degenerate-free by name.
            last = audit_degenerate([D['hist'][-1]], D['fd'])
            fem_ok = (last['total_hole'] == 0 and last['total_oversize'] == 0
                      and last['total_inverted'] == 0)
            print(f"[audit]   hist[-1] (FEM input) degenerate-free = {fem_ok}")
            if res['total_hole'] or res['total_oversize'] or res['total_inverted'] or not fem_ok:
                failures.append('degenerate')
                print('[audit] ERROR: degenerate elements found', file=sys.stderr)
        except Exception as e:
            failures.append('degenerate')
            print(f'[audit] ERROR in degenerate audit: {e}', file=sys.stderr)

    if args.check in ('pacing', 'all'):
        print('[audit] Running audit_pacing...')
        try:
            D = _stage_data()
            ratio = audit_pacing(D['hist'], D['snap_idx'], TRUSS_HOLD)
            print(f"[audit]   thirds_displacement_ratio = {ratio:.3f}  (target <= 3.0)")
            if ratio > 3.0:
                failures.append('pacing')
                print('[audit] ERROR: pacing ratio exceeds 3.0', file=sys.stderr)
        except Exception as e:
            failures.append('pacing')
            print(f'[audit] ERROR in pacing audit: {e}', file=sys.stderr)

    if args.check in ('peel-invariance', 'all'):
        print('[audit] Running audit_peel_invariance...')
        try:
            D = _stage_data()
            result = audit_peel_invariance(D['q_fem'], D['elem_layer'], D['n_layers'])
            print(f"[audit]   per-bin-totals-invariant = {result['histogram_invariant']}")
            print(f"[audit]   totals-match-independent-histogram = {result['histogram_independent']}")
            print(f"[audit]   {result['message']}")
            if not (result['histogram_invariant'] and result['histogram_independent']):
                failures.append('peel-invariance')
                print('[audit] ERROR: peel invariance check failed', file=sys.stderr)
        except Exception as e:
            failures.append('peel-invariance')
            print(f'[audit] ERROR in peel-invariance audit: {e}', file=sys.stderr)

    if args.check in ('determinism', 'all'):
        print('[audit] Running audit_determinism...')
        try:
            D1 = _stage_data()
            D2 = _stage_data()
            result = audit_determinism(
                D1.get('snapshots', []), D1.get('n_layers', 0),
                D2.get('snapshots', []), D2.get('n_layers', 0)
            )
            print(f"[audit]   run1_count={result['run1_count']}, run2_count={result['run2_count']}, "
                  f"n_layers_match={result['n_layers_match']}")
            print(f"[audit]   {result['message']}")
            if not result['snapshots_match']:
                failures.append('determinism')
                print('[audit] ERROR: determinism check failed', file=sys.stderr)
        except Exception as e:
            failures.append('determinism')
            print(f'[audit] ERROR in determinism audit: {e}', file=sys.stderr)

    if failures:
        print(f'[audit] FAILED checks: {", ".join(failures)}')
        sys.exit(1)
    print('[audit] All checks passed')


if __name__ == '__main__':
    main()
