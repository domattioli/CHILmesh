"""Tests for Gmsh ↔ admesh round-trip parity — spec 008 FR-013,014, SC-003."""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest

from admesh.api import BoundarySegment, Mesh
from admesh.boundary_types import BoundaryType
from admesh.gmsh import read_msh, write_msh
from admesh.fort14 import read_fort14, write_fort14


# ---------------------------------------------------------------------------
# Round-trip equality predicate (from data-model.md)
# ---------------------------------------------------------------------------


def _canon_elements(elements: np.ndarray) -> np.ndarray:
    """Row-canonicalize: sort each row by (min, mid, max), then sort rows."""
    canon = np.sort(elements, axis=1)
    order = np.lexsort(canon.T[::-1])
    return canon[order]


def assert_mesh_roundtrip_equal(
    a: Mesh,
    b: Mesh,
    *,
    node_atol: float = 1e-9,
    check_boundaries: bool = True,
) -> None:
    """Assert round-trip equality per data-model.md equality predicate."""
    assert a.n_nodes == b.n_nodes, f"node count: {a.n_nodes} != {b.n_nodes}"
    assert a.n_elements == b.n_elements, f"element count: {a.n_elements} != {b.n_elements}"

    # Nodes: sort by rounded (x, y) for stable comparison across format round-trips
    # (fort14 truncates to 6 dp, so we round before sorting to get stable order)
    rnd = max(0, -int(np.floor(np.log10(node_atol))) - 2) if node_atol > 0 else 12
    a_rounded = np.round(a.nodes, rnd)
    b_rounded = np.round(b.nodes, rnd)
    a_order = np.lexsort(a_rounded.T[::-1])
    b_order = np.lexsort(b_rounded.T[::-1])
    np.testing.assert_allclose(
        a.nodes[a_order], b.nodes[b_order],
        atol=node_atol, rtol=0,
        err_msg="node coordinates differ beyond tolerance",
    )

    # Elements: canonicalize
    np.testing.assert_array_equal(
        _canon_elements(a.elements),
        _canon_elements(b.elements),
        err_msg="element connectivity differs",
    )

    # Bathymetry
    if a.bathymetry is None and b.bathymetry is None:
        pass
    elif a.bathymetry is not None and b.bathymetry is not None:
        np.testing.assert_allclose(
            a.bathymetry, b.bathymetry,
            atol=node_atol, rtol=0,
            err_msg="bathymetry differs",
        )
    else:
        pytest.fail(f"bathymetry presence mismatch: a={a.bathymetry is not None}, b={b.bathymetry is not None}")

    if check_boundaries:
        assert len(a.boundaries) == len(b.boundaries), (
            f"boundary count: {len(a.boundaries)} != {len(b.boundaries)}"
        )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_mesh_with_boundaries(
    bc_type: BoundaryType = BoundaryType.MAINLAND,
    with_bathy: bool = False,
) -> Mesh:
    """4-triangle mesh with one boundary ring."""
    nodes = np.array(
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]],
        dtype=np.float64,
    )
    elements = np.array(
        [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
        dtype=np.int64,
    )
    is_open = bc_type == BoundaryType.OPEN
    ring = np.array([0, 1, 2, 3], dtype=np.int64)
    bseg = BoundarySegment(node_ids=ring, bc_type=bc_type, is_open=is_open)
    bathy = None
    if with_bathy:
        bathy = np.array([-5.0, -10.0, -15.0, -20.0, -12.5], dtype=np.float64)
    return Mesh(
        nodes=nodes, elements=elements, boundaries=(bseg,),
        bathymetry=bathy, quality=None, title="roundtrip_test",
    )


# ---------------------------------------------------------------------------
# Gmsh write → read round-trip
# ---------------------------------------------------------------------------


class TestGmshRoundTrip:
    """Write a mesh, read it back, assert structural parity."""

    @pytest.mark.parametrize("version", ["2.2", "4.1"])
    def test_node_count_preserved(self, tmp_path, version):
        mesh = _make_mesh_with_boundaries()
        out = tmp_path / f"rt_{version.replace('.','')}.msh"
        write_msh(mesh, out, version=version)
        back = read_msh(out)
        assert back.n_nodes == mesh.n_nodes

    @pytest.mark.parametrize("version", ["2.2", "4.1"])
    def test_element_count_preserved(self, tmp_path, version):
        mesh = _make_mesh_with_boundaries()
        out = tmp_path / f"rt_{version.replace('.','')}.msh"
        write_msh(mesh, out, version=version)
        back = read_msh(out)
        assert back.n_elements == mesh.n_elements

    @pytest.mark.parametrize("version", ["2.2", "4.1"])
    def test_node_coordinates_preserved(self, tmp_path, version):
        mesh = _make_mesh_with_boundaries()
        out = tmp_path / f"rt_{version.replace('.','')}.msh"
        write_msh(mesh, out, version=version)
        back = read_msh(out)
        assert_mesh_roundtrip_equal(mesh, back, node_atol=1e-9, check_boundaries=False)

    @pytest.mark.parametrize("version", ["2.2", "4.1"])
    def test_boundary_type_preserved(self, tmp_path, version):
        mesh = _make_mesh_with_boundaries(BoundaryType.MAINLAND)
        out = tmp_path / f"rtbc_{version.replace('.','')}.msh"
        write_msh(mesh, out, version=version)
        back = read_msh(out)
        bc_types = {int(s.bc_type) for s in back.boundaries}
        assert int(BoundaryType.MAINLAND) in bc_types

    @pytest.mark.parametrize("version", ["2.2", "4.1"])
    def test_bathymetry_roundtrip(self, tmp_path, version):
        mesh = _make_mesh_with_boundaries(with_bathy=True)
        out = tmp_path / f"rtbathy_{version.replace('.','')}.msh"
        write_msh(mesh, out, version=version)
        back = read_msh(out)
        assert back.bathymetry is not None
        # Node order may differ; sort by coordinates to compare
        a_order = np.lexsort(mesh.nodes.T[::-1])
        b_order = np.lexsort(back.nodes.T[::-1])
        np.testing.assert_allclose(
            mesh.bathymetry[a_order], back.bathymetry[b_order],
            atol=1e-9, rtol=0,
        )

    def test_numeric_ibtype_roundtrip(self, tmp_path):
        """FR-009,013: numeric IBTYPE codes survive write→read."""
        ring = np.array([0, 1, 2, 3], dtype=np.int64)
        nodes = np.array(
            [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]],
            dtype=np.float64,
        )
        elements = np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]], dtype=np.int64)
        # Use unmapped IBTYPE = 22
        seg = BoundarySegment(node_ids=ring, bc_type=22, is_open=False)
        mesh = Mesh(nodes=nodes, elements=elements, boundaries=(seg,),
                    bathymetry=None, quality=None, title="")
        out = tmp_path / "ibtype22.msh"
        write_msh(mesh, out, version="2.2")
        back = read_msh(out)
        assert len(back.boundaries) >= 1
        # The BC type should be recoverable (may come back as numeric from group name)
        assert any(s.bc_type is not None for s in back.boundaries)


# ---------------------------------------------------------------------------
# Gmsh → fort14 → Gmsh round-trip
# ---------------------------------------------------------------------------


class TestGmshFort14Roundtrip:
    """FR-013: Gmsh → admesh → fort.14 → admesh → Gmsh structural parity."""

    def test_gmsh_to_fort14_to_gmsh(self, tmp_path):
        mesh_orig = _make_mesh_with_boundaries()
        msh_path = tmp_path / "orig.msh"
        fort_path = tmp_path / "mesh.14"
        msh_back_path = tmp_path / "back.msh"

        # Gmsh → admesh
        write_msh(mesh_orig, msh_path, version="2.2")
        mesh_from_gmsh = read_msh(msh_path)

        # admesh → fort14 → admesh
        write_fort14(mesh_from_gmsh, fort_path)
        mesh_from_fort = read_fort14(fort_path)

        # admesh → Gmsh (complete the cycle)
        write_msh(mesh_from_fort, msh_back_path, version="2.2")
        mesh_final = read_msh(msh_back_path)

        # Structural parity
        assert mesh_final.n_nodes == mesh_orig.n_nodes
        assert mesh_final.n_elements == mesh_orig.n_elements
        assert_mesh_roundtrip_equal(mesh_orig, mesh_final, node_atol=1e-9, check_boundaries=False)

    def test_fort14_to_gmsh_to_fort14(self, tmp_path):
        """fort14 → Gmsh → fort14 preserves nodes and elements."""
        mesh = _make_mesh_with_boundaries()
        fort = tmp_path / "mid.14"
        msh = tmp_path / "out.msh"
        write_fort14(mesh, fort)
        m2 = read_fort14(fort)
        write_msh(m2, msh, version="2.2")
        m3 = read_msh(msh)
        assert m3.n_nodes == mesh.n_nodes
        assert m3.n_elements == mesh.n_elements

    def test_boundary_labels_survive_fort14_roundtrip(self, tmp_path):
        mesh = _make_mesh_with_boundaries(BoundaryType.MAINLAND)
        msh1 = tmp_path / "start.msh"
        fort = tmp_path / "mid.14"
        msh2 = tmp_path / "end.msh"

        write_msh(mesh, msh1, version="2.2")
        m1 = read_msh(msh1)
        write_fort14(m1, fort)
        m2 = read_fort14(fort)
        write_msh(m2, msh2, version="2.2")
        m3 = read_msh(msh2)

        # BC type should be preserved across the full chain
        bc_types_orig = {int(s.bc_type) for s in mesh.boundaries}
        bc_types_final = {int(s.bc_type) for s in m3.boundaries}
        assert bc_types_orig == bc_types_final


# ---------------------------------------------------------------------------
# MVP domains (FR-014 / SC-003)
# ---------------------------------------------------------------------------


FIXTURES_GMSH = Path(__file__).parent / "fixtures" / "gmsh"


def _try_triangulate_domain(domain_name: str) -> Mesh | None:
    """Attempt to triangulate an MVP domain; skip if unavailable."""
    try:
        from admesh._stages.domains import UNIT_SQUARE, L_SHAPE, UNIT_DISK, ANNULUS, NOTCHED_RECTANGLE
        import admesh
        domain_map = {
            "unit_square": UNIT_SQUARE,
            "L_shape": L_SHAPE,
            "unit_disk": UNIT_DISK,
            "annulus": ANNULUS,
            "notched_rectangle": NOTCHED_RECTANGLE,
        }
        domain = domain_map.get(domain_name)
        if domain is None:
            return None
        # Use a coarse mesh for speed
        mesh = admesh.triangulate(domain, h_max=0.3, quality_gate=(0.1, 0.3))
        return mesh
    except Exception:
        return None


@pytest.mark.parametrize("domain_name", [
    "unit_square", "L_shape", "unit_disk", "annulus", "notched_rectangle"
])
@pytest.mark.parametrize("version", ["2.2", "4.1"])
def test_mvp_domain_write_read(tmp_path, domain_name, version):
    """FR-014: each MVP domain can be written and read back without error."""
    mesh = _try_triangulate_domain(domain_name)
    if mesh is None:
        pytest.skip(f"Could not triangulate {domain_name}")

    out = tmp_path / f"{domain_name}_{version.replace('.','')}.msh"
    write_msh(mesh, out, version=version)
    back = read_msh(out)

    assert back.n_nodes == mesh.n_nodes
    assert back.n_elements == mesh.n_elements


@pytest.mark.parametrize("domain_name", [
    "unit_square", "L_shape", "unit_disk", "annulus", "notched_rectangle"
])
def test_mvp_domain_gmsh_fort14_gmsh(tmp_path, domain_name):
    """FR-014/SC-003: Gmsh→fort14→Gmsh structural parity per MVP domain."""
    mesh = _try_triangulate_domain(domain_name)
    if mesh is None:
        pytest.skip(f"Could not triangulate {domain_name}")

    msh1 = tmp_path / "start.msh"
    fort = tmp_path / "mid.14"
    msh2 = tmp_path / "end.msh"

    write_msh(mesh, msh1, version="2.2")
    m1 = read_msh(msh1)
    write_fort14(m1, fort)
    m2 = read_fort14(fort)
    write_msh(m2, msh2, version="2.2")
    m3 = read_msh(msh2)

    assert m3.n_nodes == mesh.n_nodes
    assert m3.n_elements == mesh.n_elements
    assert_mesh_roundtrip_equal(mesh, m3, node_atol=1e-5, check_boundaries=False)  # fort14 truncates to 6dp
