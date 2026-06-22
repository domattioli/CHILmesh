"""Run the release-integrity eval corpus and emit a scorecard."""
from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "scripts"))
import release_integrity as ri  # noqa: E402

from cases import CASES  # noqa: E402


def _verdict(case) -> int:
    if case["kind"] == "api":
        code, _ = ri.run_api_gate(case["breaking"], case["old"], case["new"], case["mode"])
        return code
    code, _ = ri.run_release_check(case["tag"], case["version"], case["changelog"], case["twine_ok"])
    return code


def score(cases=CASES) -> dict:
    total = len(cases)
    correct = sum(1 for c in cases if _verdict(c) == c["expect_code"])
    return {"total": total, "correct": correct, "accuracy": correct / total}


if __name__ == "__main__":
    s = score()
    print(f"release-integrity eval: {s['correct']}/{s['total']} correct "
          f"(accuracy {s['accuracy']:.0%})")
    sys.exit(0 if s["correct"] == s["total"] else 1)
