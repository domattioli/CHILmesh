"""Pin CHILmesh public unification API for MADMESHing/QuADMesh consumers (spec-048 T1).

CHILmesh serves as a bridge library downstream of ADMESH, providing mesh data
structures and skeletonization for MADMESHR and QuADMesh. This contract test
ensures that critical functions remain stable and importable across versions
for callers like MADMESHing (https://github.com/domattioli/MADMESHing).

Reference: domattioli/CHILmesh#196
"""

from __future__ import annotations

import dataclasses
import inspect
import tempfile
from pathlib import Path

import numpy as np


def test_unification_surface_importable():
    """CHILmesh, element_quality, paths_on_outer_vertices are importable."""
    from chilmesh import CHILmesh, element_quality, paths_on_outer_vertices

    assert CHILmesh is not None
    assert element_quality is not None
    assert paths_on_outer_vertices is not None


def test_unification_surface_in_all():
    """CHILmesh, element_quality, paths_on_outer_vertices are in chilmesh.__all__."""
    import chilmesh

    expected = {"CHILmesh", "element_quality", "paths_on_outer_vertices"}
    assert expected.issubset(set(chilmesh.__all__)), (
        f"Missing from chilmesh.__all__: {expected - set(chilmesh.__all__)}"
    )


def test_canonical_modules():
    """Symbols live in pinned canonical modules."""
    from chilmesh import CHILmesh, element_quality, paths_on_outer_vertices

    assert CHILmesh.__module__ == "chilmesh.CHILmesh", (
        f"CHILmesh must live in chilmesh.CHILmesh, found in {CHILmesh.__module__}"
    )
    assert element_quality.__module__ == "chilmesh.quality", (
        f"element_quality must live in chilmesh.quality, "
        f"found in {element_quality.__module__}"
    )
    assert paths_on_outer_vertices.__module__ == "chilmesh.layer_paths", (
        f"paths_on_outer_vertices must live in chilmesh.layer_paths, "
        f"found in {paths_on_outer_vertices.__module__}"
    )


def test_element_quality_signature():
    """element_quality() has required signature for consumers.

    Consumers call: element_quality(verts, conn, metric='aspect_ratio')
    Params 'verts' and 'conn' must be named exactly; 'metric' has default.
    """
    from chilmesh import element_quality

    sig = inspect.signature(element_quality)
    params = sig.parameters

    param_list = list(params.keys())
    assert len(param_list) >= 2, (
        f"element_quality must have at least 2 parameters, got {len(param_list)}"
    )

    assert param_list[0] == "verts", (
        f"First param must be 'verts', got {param_list[0]}"
    )
    assert param_list[1] == "conn", (
        f"Second param must be 'conn', got {param_list[1]}"
    )

    assert "metric" in params, "Missing 'metric' parameter"
    assert params["metric"].default != inspect.Parameter.empty, (
        "metric must have a default value"
    )


def test_paths_on_outer_vertices_signature():
    """paths_on_outer_vertices(mesh, layer_idx) has required signature."""
    from chilmesh import paths_on_outer_vertices

    sig = inspect.signature(paths_on_outer_vertices)
    params = sig.parameters

    param_list = list(params.keys())
    assert len(param_list) == 2, (
        f"paths_on_outer_vertices must have exactly 2 parameters, "
        f"got {len(param_list)}: {param_list}"
    )

    assert param_list[0] == "mesh", (
        f"First param must be 'mesh', got {param_list[0]}"
    )
    assert param_list[1] == "layer_idx", (
        f"Second param must be 'layer_idx', got {param_list[1]}"
    )


def test_chilmesh_submesh_signature():
    """CHILmesh.submesh() has required signature.

    Consumers call: mesh.submesh(elem_ids, compute_layers=..., compute_adjacencies=...)
    """
    from chilmesh import CHILmesh

    sig = inspect.signature(CHILmesh.submesh)
    params = sig.parameters

    # Must have 'elem_ids' as first positional (after self)
    param_list = list(params.keys())
    assert "elem_ids" in param_list, (
        f"submesh must have 'elem_ids' parameter, got {param_list}"
    )

    # compute_layers and compute_adjacencies should be present
    assert "compute_layers" in params, "Missing 'compute_layers' parameter"
    assert "compute_adjacencies" in params, "Missing 'compute_adjacencies' parameter"

    # Both should have defaults
    assert params["compute_layers"].default != inspect.Parameter.empty, (
        "compute_layers must have a default value"
    )
    assert params["compute_adjacencies"].default != inspect.Parameter.empty, (
        "compute_adjacencies must have a default value"
    )


def test_chilmesh_ccw_edges_around_vert_signature():
    """CHILmesh.ccw_edges_around_vert(self, vert_id) has required signature."""
    from chilmesh import CHILmesh

    sig = inspect.signature(CHILmesh.ccw_edges_around_vert)
    params = sig.parameters

    param_list = list(params.keys())
    assert "vert_id" in param_list, (
        f"ccw_edges_around_vert must have 'vert_id' parameter, got {param_list}"
    )


def test_chilmesh_init_compute_flags():
    """CHILmesh.__init__() accepts compute_layers and compute_adjacencies kwargs."""
    from chilmesh import CHILmesh

    sig = inspect.signature(CHILmesh.__init__)
    params = sig.parameters

    assert "compute_layers" in params, (
        "CHILmesh.__init__ missing 'compute_layers' parameter"
    )
    assert "compute_adjacencies" in params, (
        "CHILmesh.__init__ missing 'compute_adjacencies' parameter"
    )

    # Both should have defaults
    assert params["compute_layers"].default != inspect.Parameter.empty, (
        "compute_layers must have a default value"
    )
    assert params["compute_adjacencies"].default != inspect.Parameter.empty, (
        "compute_adjacencies must have a default value"
    )


def test_chilmesh_save_load_exist():
    """CHILmesh.save() and CHILmesh.load() exist with required signatures."""
    from chilmesh import CHILmesh

    assert hasattr(CHILmesh, "save"), "CHILmesh must have save method"
    assert hasattr(CHILmesh, "load"), "CHILmesh must have load method"
    assert callable(CHILmesh.save), "save must be callable"
    assert callable(CHILmesh.load), "load must be callable"

    # Verify save signature: (self, filename)
    save_sig = inspect.signature(CHILmesh.save)
    save_params = list(save_sig.parameters.keys())
    assert "filename" in save_params, (
        f"save must have 'filename' parameter, got {save_params}"
    )

    # Verify load signature: (cls, filename, compute_layers=..., compute_adjacencies=...)
    load_sig = inspect.signature(CHILmesh.load)
    load_params = list(load_sig.parameters.keys())
    assert "filename" in load_params, (
        f"load must have 'filename' parameter, got {load_params}"
    )
    assert "compute_layers" in load_params, (
        f"load must have 'compute_layers' parameter, got {load_params}"
    )
    assert "compute_adjacencies" in load_params, (
        f"load must have 'compute_adjacencies' parameter, got {load_params}"
    )


def _simple_mesh():
    """Create a simple two-triangle mesh for functional tests."""
    pts = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float)
    conn = np.array([[0, 1, 2], [0, 2, 3]], dtype=int)
    return conn, pts


def test_adjacencies_structure():
    """mesh.adjacencies dict contains required keys."""
    from chilmesh import CHILmesh

    conn, pts = _simple_mesh()
    mesh = CHILmesh(
        connectivity=conn,
        points=pts,
        compute_layers=False,
        compute_adjacencies=True,
    )

    # Required adjacency keys (EdgeMap is internal optimization, also present)
    required_keys = {"Elem2Vert", "Edge2Vert", "Elem2Edge", "Vert2Edge", "Vert2Elem", "Edge2Elem"}
    assert "adjacencies" in dir(mesh), "mesh must have 'adjacencies' attribute"
    assert isinstance(mesh.adjacencies, dict), "adjacencies must be a dict"

    missing = required_keys - set(mesh.adjacencies.keys())
    assert not missing, (
        f"adjacencies dict missing required keys: {missing}. "
        f"Found: {set(mesh.adjacencies.keys())} (#196 issue: CHILmesh#196)"
    )


def test_layers_structure():
    """mesh.layers dict contains required per-layer keys."""
    from chilmesh import CHILmesh

    conn, pts = _simple_mesh()
    mesh = CHILmesh(
        connectivity=conn,
        points=pts,
        compute_layers=True,
        compute_adjacencies=True,
    )

    assert "layers" in dir(mesh), "mesh must have 'layers' attribute"
    assert isinstance(mesh.layers, dict), "layers must be a dict"

    # Per-layer structure: each key maps to a list/dict of arrays
    required_layer_keys = {"OE", "IE", "OV", "IV"}
    assert "OE" in mesh.layers, "layers must have 'OE' (outer elements)"
    assert "IE" in mesh.layers, "layers must have 'IE' (inner elements)"
    assert "OV" in mesh.layers, "layers must have 'OV' (outer vertices)"
    assert "IV" in mesh.layers, "layers must have 'IV' (inner vertices)"


def test_element_quality_functional():
    """element_quality() returns correct shape on simple mesh."""
    from chilmesh import element_quality

    verts = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], dtype=float)
    conn = np.array([[0, 1, 2], [0, 2, 3]], dtype=int)

    result = element_quality(verts, conn, metric="aspect_ratio")

    assert isinstance(result, np.ndarray), "element_quality must return ndarray"
    assert result.shape == (2,), (
        f"element_quality shape must be (2,) for 2 elements, got {result.shape}"
    )
    assert np.all(result >= 0), "quality values must be non-negative"
    assert np.all(result <= 1), "aspect_ratio must be in [0, 1]"


def test_submesh_functional():
    """CHILmesh.submesh() returns valid mesh for valid element selection."""
    from chilmesh import CHILmesh

    conn, pts = _simple_mesh()
    mesh = CHILmesh(
        connectivity=conn,
        points=pts,
        compute_layers=False,
        compute_adjacencies=True,
    )

    # Extract just the first element
    sub = mesh.submesh([0], compute_layers=False, compute_adjacencies=True)

    assert isinstance(sub, CHILmesh), "submesh must return CHILmesh instance"
    assert sub.n_elems == 1, f"submesh should have 1 element, got {sub.n_elems}"
    assert sub.n_verts >= 3, (
        f"submesh with 1 triangle should have at least 3 vertices, got {sub.n_verts}"
    )


def test_ccw_edges_around_vert_functional():
    """ccw_edges_around_vert() returns list of edge IDs."""
    from chilmesh import CHILmesh

    conn, pts = _simple_mesh()
    mesh = CHILmesh(
        connectivity=conn,
        points=pts,
        compute_layers=False,
        compute_adjacencies=True,
    )

    # Vertex 0 is incident to edges in both triangles
    edges = mesh.ccw_edges_around_vert(0)

    assert isinstance(edges, list), (
        f"ccw_edges_around_vert must return list, got {type(edges)}"
    )
    assert all(isinstance(e, (int, np.integer)) for e in edges), (
        "ccw_edges_around_vert must return list of integers"
    )


def test_save_load_roundtrip():
    """CHILmesh.save() + CHILmesh.load() roundtrip preserves connectivity and points."""
    from chilmesh import CHILmesh

    conn, pts = _simple_mesh()
    mesh = CHILmesh(
        connectivity=conn,
        points=pts,
        compute_layers=False,
        compute_adjacencies=False,
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        filepath = Path(tmpdir) / "test_mesh.14"
        mesh.save(str(filepath))

        loaded = CHILmesh.load(str(filepath), compute_layers=False, compute_adjacencies=False)

        # Verify connectivity preserved
        np.testing.assert_array_equal(
            loaded.connectivity_list, mesh.connectivity_list,
            err_msg="save/load roundtrip corrupted connectivity"
        )

        # Verify points preserved (allow small floating-point deviation)
        np.testing.assert_allclose(
            loaded.points, mesh.points,
            rtol=1e-10, atol=1e-15,
            err_msg="save/load roundtrip corrupted points"
        )


def test_save_load_dispatch_rejects_unknown_extension():
    """CHILmesh.save() and load() reject unknown file extensions."""
    import pytest
    from chilmesh import CHILmesh

    conn, pts = _simple_mesh()
    mesh = CHILmesh(
        connectivity=conn,
        points=pts,
        compute_layers=False,
        compute_adjacencies=False,
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)

        # Test save with unknown extension
        unknown_save_path = str(tmp / "test_mesh.xyz")
        with pytest.raises(ValueError):
            mesh.save(unknown_save_path)

        # Test load with unknown extension on an existing file
        dummy_xyz_path = tmp / "dummy.xyz"
        dummy_xyz_path.write_text("dummy content")
        with pytest.raises(ValueError):
            CHILmesh.load(str(dummy_xyz_path))
