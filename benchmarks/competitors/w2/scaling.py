"""Test the "O(log n) point location" claim: fit query time against log n and against n (FR-013, SC-004)."""
from __future__ import annotations

from typing import Sequence

import numpy as np

CONSTANT_SPREAD = 0.05  # (max-min)/mean below this counts as flat

NOTE = ("Least-squares fits of time = a + b*log(n) and time = a + b*n. O(log n) is supported only "
        "when the log fit has the higher R^2. Cache effects (tree and element data outgrowing CPU "
        "cache) can bend the curve at large n; interpret with the per-size timings.")


def _r2(x: np.ndarray, y: np.ndarray) -> float:
    slope, intercept = np.polyfit(x, y, 1)
    resid = y - (slope * x + intercept)
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return 1.0 if ss_tot == 0 else 1.0 - ss_res / ss_tot


def fit_scaling(sizes: Sequence[float], times: Sequence[float]) -> dict:
    n = np.asarray(sizes, dtype=float)
    t = np.asarray(times, dtype=float)
    if len(n) != len(t):
        raise ValueError("sizes and times must have the same length")
    if len(n) < 3:
        raise ValueError("need at least 3 mesh sizes (SC-004)")
    if np.any(n <= 0):
        raise ValueError("sizes must be positive")
    base = {"sizes": [int(s) if float(s).is_integer() else float(s) for s in sizes],
            "times": [float(x) for x in t], "note": NOTE}
    if t.mean() > 0 and (t.max() - t.min()) / t.mean() < CONSTANT_SPREAD:
        return {**base, "best": "constant", "r2_log": None, "r2_linear": None, "o_log_n_supported": False}
    r2_log, r2_lin = _r2(np.log(n), t), _r2(n, t)
    best = "log-n" if r2_log > r2_lin else "linear-n"
    return {**base, "best": best, "r2_log": r2_log, "r2_linear": r2_lin,
            "o_log_n_supported": best == "log-n"}
