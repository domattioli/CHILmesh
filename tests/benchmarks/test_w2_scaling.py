"""Scaling fit: log-n versus linear-n on synthetic timings (FR-013, SC-004)."""
from __future__ import annotations

import math

import pytest

from benchmarks.competitors.w2 import scaling

SIZES = [5_780, 98_365, 531_680, 24_875_336]


def test_log_growth_is_supported():
    times = [2e-5 + 1e-5 * math.log(n) for n in SIZES]
    fit = scaling.fit_scaling(SIZES, times)
    assert fit["best"] == "log-n"
    assert fit["o_log_n_supported"] is True
    assert fit["r2_log"] > fit["r2_linear"]


def test_linear_growth_refutes_log():
    times = [1e-9 * n + 1e-5 for n in SIZES]
    fit = scaling.fit_scaling(SIZES, times)
    assert fit["best"] == "linear-n"
    assert fit["o_log_n_supported"] is False


def test_needs_three_sizes():
    with pytest.raises(ValueError):
        scaling.fit_scaling(SIZES[:2], [1.0, 2.0])


def test_mismatched_lengths_rejected():
    with pytest.raises(ValueError):
        scaling.fit_scaling(SIZES, [1.0, 2.0, 3.0])


def test_nonpositive_size_rejected():
    with pytest.raises(ValueError):
        scaling.fit_scaling([0, 10, 100], [1.0, 2.0, 3.0])


def test_constant_times_not_called_log():
    fit = scaling.fit_scaling(SIZES, [1e-4] * 4)
    assert fit["best"] == "constant"
    assert fit["o_log_n_supported"] is False


def test_result_lists_inputs_and_cache_note():
    fit = scaling.fit_scaling(SIZES, [2e-5 + 1e-5 * math.log(n) for n in SIZES])
    assert fit["sizes"] == SIZES
    assert "cache" in fit["note"].lower()
