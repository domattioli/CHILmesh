"""Spec 007 — mesh element validity test suite.

One test per fixture asserting ``report.ok`` (FR-017). Failure messages
aggregate all violation categories via ``format_failures``.
"""
from __future__ import annotations

import pytest

import chilmesh
from chilmesh.tri2quad import tri_to_quad
from tests._validity import validate_mesh_elements
from tests._validity.fixtures import (
    bowtie_quad_mesh,
    edge_crossing_mesh,
    interior_triangle_mesh,
    overlapping_quads_mesh,
    pentagon_mesh,
)
from tests._validity.validator import format_failures

BUILTIN_FIXTURES = ["annulus", "donut", "structured"]

RUNTIME_BUDGET_S = {
    "annulus": 10.0,
    "donut": 15.0,
    "structured": 15.0,
    "block_o": 90.0,
}


# annulus & donut: dual-graph of input triangulation has Tutte-Berge
# obstructions — a few interior tris cannot be paired by any matching
# without non-conforming subdivision. Tracked separately; not required
# by the spec-007 acceptance for block_o.
TOPOLOGY_OBSTRUCTED = {"annulus", "donut"}


@pytest.fixture(scope="module", params=BUILTIN_FIXTURES)
def quadified_builtin(request):
    raw = getattr(chilmesh.examples, request.param)()
    try:
        return tri_to_quad(raw), request.param, None
    except RuntimeError as exc:
        return None, request.param, exc


def test_builtin_fixtures_pass(quadified_builtin):
    mesh, name, exc = quadified_builtin
    if exc is not None:
        if name in TOPOLOGY_OBSTRUCTED:
            pytest.xfail(f"{name}: max matching leaves interior tris ({exc})")
        raise exc
    report = validate_mesh_elements(mesh)
    assert report.ok, format_failures(report)
    assert report.runtime_s < RUNTIME_BUDGET_S[name], (
        f"{name}: validator took {report.runtime_s:.2f}s (budget {RUNTIME_BUDGET_S[name]}s)"
    )


@pytest.fixture(scope="module")
def quadified_block_o():
    raw = chilmesh.examples.block_o()
    return tri_to_quad(raw)


@pytest.mark.slow
def test_block_o_passes(quadified_block_o):
    report = validate_mesh_elements(quadified_block_o)
    assert report.ok, format_failures(report)
    assert report.runtime_s < RUNTIME_BUDGET_S["block_o"], (
        f"block_o: validator took {report.runtime_s:.2f}s (budget {RUNTIME_BUDGET_S['block_o']}s)"
    )


@pytest.mark.parametrize(
    "factory,expected_category",
    [
        (bowtie_quad_mesh, "SELF_INTERSECTING_QUAD"),
        (interior_triangle_mesh, "INTERIOR_TRIANGLE_FORBIDDEN"),
        (pentagon_mesh, "UNSUPPORTED_ELEMENT_ARITY"),
        (overlapping_quads_mesh, "INTERIOR_OVERLAP"),
        (edge_crossing_mesh, "EDGE_CROSSING"),
    ],
)
def test_synthetic_negatives(factory, expected_category):
    mesh = factory()
    report = validate_mesh_elements(mesh)
    cats = {v.category for v in report.violations}
    assert expected_category in cats, (
        f"expected {expected_category} in report; got {cats}. "
        + format_failures(report)
    )
    assert not report.ok


def test_failure_message_format():
    mesh = bowtie_quad_mesh()
    report = validate_mesh_elements(mesh)
    msg = format_failures(report)
    assert "SELF_INTERSECTING_QUAD" in msg
    assert "elems=" in msg


def test_degenerate_quad_accepted():
    """Padded-triangle and duplicate-vertex quads are notes, not violations."""
    import numpy as np

    pts = np.array(
        [
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [2.0, 1.0, 0.0],
        ]
    )
    conn = np.array(
        [
            [0, 1, 4, 3],
            [1, 2, 5, 4],
            [3, 4, 1, -1],
        ]
    )
    mesh = chilmesh.CHILmesh(connectivity=conn, points=pts, compute_layers=True)
    report = validate_mesh_elements(mesh)
    assert report.ok, format_failures(report)
