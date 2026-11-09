"""time_adjacency must time construction with layers.

Adjacency building is pure Python on every backend. Only layerization differs between
the C++ and Python backends, so adj_cpp_s / adj_python_s only differ if layers are built.
No chilmesh_cpp needed.
"""
from __future__ import annotations

import inspect
from pathlib import Path

import numpy as np
import pytest

import chilmesh
from benchmarks.competitors.w2 import measures, tools

REPO = Path(__file__).resolve().parents[2]
README = REPO / "benchmarks" / "competitors" / "README.md"


def _tiny():
    pts = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [2.0, 0.5, 0.0]])
    conn = np.array([[0, 1, 2], [0, 2, 3], [1, 4, 2]])
    return conn, pts


def test_time_adjacency_builds_layers(monkeypatch):
    seen = {}
    from chilmesh import CHILmesh

    real_init = CHILmesh.__init__

    def spy(self, *args, **kwargs):
        seen.update(kwargs)
        return real_init(self, *args, **kwargs)

    monkeypatch.setattr(CHILmesh, "__init__", spy)
    measures.time_adjacency(*_tiny())
    assert seen.get("compute_layers") is True
    assert seen.get("compute_adjacencies") is True


def test_time_adjacency_keeps_result_keys():
    out = measures.time_adjacency(*_tiny())
    assert {"adjacency_s", "backend_selected", "backend_info"} <= set(out)
    assert out["adjacency_s"] >= 0


def test_time_adjacency_docstring_says_construction_with_layers():
    doc = inspect.getdoc(measures.time_adjacency).lower()
    assert "layer" in doc
    assert "construction" in doc


def test_time_adjacency_cpp_backend_info_when_available():
    if "cpp" not in str(chilmesh.backend_info()).lower():
        pytest.skip("chilmesh_cpp not installed")
    assert measures.time_adjacency(*_tiny())["backend_selected"] is not None


def test_tools_describe_adj_metrics_as_construction_with_layers():
    desc = tools.METRIC_DESCRIPTIONS
    for name in ("adj_cpp_s", "adj_python_s"):
        text = desc[name].lower()
        assert "construction" in text and "layers" in text
        assert "adjacency build" not in text
    assert "c++" in desc["adj_cpp_s"].lower()
    assert "python" in desc["adj_python_s"].lower()


def test_metric_names_unchanged():
    assert "adj_cpp_s" in tools.SECONDARY and "adj_python_s" in tools.SECONDARY


def test_readme_describes_construction_with_layers():
    text = README.read_text(encoding="utf-8")
    assert "adjacency build" not in text.lower()
    assert "adj_cpp_s" in text and "adj_python_s" in text
    low = text.lower()
    assert "construction with adjacencies and layers" in low
    assert "c++" in low
