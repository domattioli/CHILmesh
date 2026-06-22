from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import run_eval  # noqa: E402


def test_corpus_full_accuracy():
    s = run_eval.score()
    assert s["correct"] == s["total"], f"eval regressions: {s}"
