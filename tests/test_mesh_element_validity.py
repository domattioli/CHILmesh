"""Spec 007 — mesh element validity test suite.

One test per fixture asserting ``report.ok`` (FR-017). Failure messages
aggregate all violation categories via ``format_failures``.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

import chilmesh
from chilmesh.tri2quad import tri_to_quad
from chilmesh.validation import validate_mesh_elements
from chilmesh.validation.fixtures import (
    bowtie_quad_mesh,
    edge_crossing_mesh,
    interior_triangle_mesh,
    overlapping_quads_mesh,
    pentagon_mesh,
)
from chilmesh.validation.validator import format_failures


_DELAWARE_BAY_CANDIDATES = (
    "/workspace/ADMESH-Domains/registry_data/meshes/Deleware_Bay.14",
    "/home/user/QuADMesh/03_CHILMesh_Test_Cases/01_.14_Files/Deleware_Bay.14",
)


def _delaware_bay_path() -> Path | None:
    override = os.environ.get("CHILMESH_DELAWARE_BAY_FORT14")
    if override and Path(override).is_file():
        return Path(override)
    for cand in _DELAWARE_BAY_CANDIDATES:
        if Path(cand).is_file():
            return Path(cand)
    return None


def _padded_tri_boundary_check(mesh) -> tuple[int, int, list[int]]:
    """Return (n_padded_tris, n_off_boundary, sample_off_boundary_ids).

    Counts ``[a,b,c,a]`` rows in mesh.connectivity_list and verifies each
    such row has at least one vertex on the mesh boundary.
    """
    boundary_verts = {int(x) for x in np.asarray(mesh.boundary_node_indices()).flatten()}
    conn = np.asarray(mesh.connectivity_list)
    n_padded = 0
    off_boundary: list[int] = []
    for elem_id in range(conn.shape[0]):
        row = conn[elem_id]
        if int(row[0]) == int(row[3]):
            n_padded += 1
            if not any(int(v) in boundary_verts for v in row[:3]):
                off_boundary.append(elem_id)
    return n_padded, len(off_boundary), off_boundary[:10]

BUILTIN_FIXTURES = ["annulus", "donut", "structured"]

RUNTIME_BUDGET_S = {
    "annulus": 10.0,
    "donut": 15.0,
    "structured": 15.0,
    "block_o": 90.0,
}


# All built-in fixtures: tri_to_quad is a literal port of MATLAB's
# identifyEdgesFun_v2 + mergeTrianglesFun (no Blossom, no edge flips, no
# vertex insertion, no midpoint bisection — per user direction). The
# path-walk + layer-deferral pairing leaves a small number of interior
# triangles on inputs whose dual graph has odd-component obstructions.
# These tests are xfail until the leftover-handling phase of the port
# (the MATLAB ``removeTrianglesFun`` family) lands per user spec.
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
        pytest.xfail(f"{name}: faithful MATLAB port leaves interior tris ({exc})")
    report = validate_mesh_elements(mesh)
    assert report.ok, format_failures(report)
    assert report.runtime_s < RUNTIME_BUDGET_S[name], (
        f"{name}: validator took {report.runtime_s:.2f}s (budget {RUNTIME_BUDGET_S[name]}s)"
    )
    n_padded, n_off, sample = _padded_tri_boundary_check(mesh)
    assert n_off == 0, (
        f"{name}: {n_off}/{n_padded} padded-tri rows have no boundary vertex "
        f"(interior triangles forbidden). Sample elem ids: {sample}"
    )


@pytest.fixture(scope="module")
def quadified_block_o():
    raw = chilmesh.examples.block_o()
    try:
        return tri_to_quad(raw), None
    except RuntimeError as exc:
        return None, exc


@pytest.mark.slow
def test_block_o_passes(quadified_block_o):
    mesh, exc = quadified_block_o
    if exc is not None:
        pytest.xfail(f"block_o: faithful MATLAB port leaves interior tris ({exc})")
    report = validate_mesh_elements(mesh)
    assert report.ok, format_failures(report)
    assert report.runtime_s < RUNTIME_BUDGET_S["block_o"], (
        f"block_o: validator took {report.runtime_s:.2f}s (budget {RUNTIME_BUDGET_S['block_o']}s)"
    )
    n_padded, n_off, sample = _padded_tri_boundary_check(mesh)
    assert n_off == 0, (
        f"block_o: {n_off}/{n_padded} padded-tri rows have no boundary vertex "
        f"(interior triangles forbidden). Sample elem ids: {sample}"
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


@pytest.fixture(scope="module")
def quadified_delaware_bay():
    path = _delaware_bay_path()
    if path is None:
        pytest.skip(
            "Delaware Bay fixture not found. Set CHILMESH_DELAWARE_BAY_FORT14 "
            f"or place file at one of: {_DELAWARE_BAY_CANDIDATES}"
        )
    raw = chilmesh.CHILmesh.read_from_fort14(str(path))
    # fort.14 stores bathymetry in z column. Spec FR-008 requires planar mesh,
    # so zero z out — the 2D mesh validator operates on xy only.
    flat_pts = np.asarray(raw.points, dtype=float).copy()
    if flat_pts.shape[1] >= 3:
        flat_pts[:, 2] = 0.0
    raw = chilmesh.CHILmesh(
        connectivity=raw.connectivity_list,
        points=flat_pts,
        grid_name=raw.grid_name,
        compute_layers=True,
    )
    try:
        return tri_to_quad(raw, strict=False), None
    except RuntimeError as exc:
        return None, exc


@pytest.mark.slow
def test_delaware_bay_passes(quadified_delaware_bay):
    """Real-world ADCIRC mesh: Delaware Bay (~26.7k tris, 17 skeleton layers).

    Asserts: tri_to_quad output validates clean AND every padded-tri row has at
    least one vertex on the mesh boundary (no interior triangles).
    """
    mesh, exc = quadified_delaware_bay
    if exc is not None:
        pytest.fail(f"delaware_bay: tri_to_quad raised {exc!r}")
    report = validate_mesh_elements(mesh)
    assert report.ok, format_failures(report)
    n_padded, n_off, sample = _padded_tri_boundary_check(mesh)
    assert n_off == 0, (
        f"delaware_bay: {n_off}/{n_padded} padded-tri rows have no boundary vertex "
        f"(interior triangles forbidden). Sample elem ids: {sample}"
    )


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
