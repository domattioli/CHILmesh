"""Round-trip tests for fort.13 nodal attribute I/O."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pytest

from chilmesh import read_fort13, write_fort13, Fort13, NodalAttribute


def test_read_fort13_parses_fixture():
    """Test reading the sample.13 fixture file."""
    fixture_path = Path(__file__).parent / "fixtures" / "fort13" / "sample.13"

    f13 = read_fort13(fixture_path)

    # Check header
    assert f13.grid_name == "sample_grid"
    assert f13.num_nodes == 4
    assert len(f13.attributes) == 2

    # Check attribute 1 (primitive_weighting_in_continuity_equation)
    attr1 = f13.attributes[0]
    assert attr1.name == "primitive_weighting_in_continuity_equation"
    assert attr1.units == "unitless"
    assert attr1.values_per_node == 1
    np.testing.assert_array_almost_equal(attr1.default_values, [0.03])
    assert 1 in attr1.nondefault  # 0-based index (node 2 in 1-based)
    np.testing.assert_array_almost_equal(attr1.nondefault[1], [0.02])
    assert len(attr1.nondefault) == 1

    # Check attribute 2 (surface_directional_effective_roughness_length)
    attr2 = f13.attributes[1]
    assert attr2.name == "surface_directional_effective_roughness_length"
    assert attr2.units == "m"
    assert attr2.values_per_node == 3
    np.testing.assert_array_almost_equal(attr2.default_values, [0.0, 0.0, 0.0])
    assert 0 in attr2.nondefault  # 0-based index (node 1 in 1-based)
    assert 3 in attr2.nondefault  # 0-based index (node 4 in 1-based)
    np.testing.assert_array_almost_equal(attr2.nondefault[0], [0.1, 0.2, 0.3])
    np.testing.assert_array_almost_equal(attr2.nondefault[3], [0.4, 0.5, 0.6])
    assert len(attr2.nondefault) == 2


def test_dense_overlay():
    """Test the dense() method overlays nondefault on default values."""
    fixture_path = Path(__file__).parent / "fixtures" / "fort13" / "sample.13"
    f13 = read_fort13(fixture_path)

    # Check dense for attribute 1
    dense1 = f13.dense("primitive_weighting_in_continuity_equation")
    assert dense1.shape == (4, 1)
    np.testing.assert_array_almost_equal(dense1[0, 0], 0.03)
    np.testing.assert_array_almost_equal(dense1[1, 0], 0.02)  # nondefault
    np.testing.assert_array_almost_equal(dense1[2, 0], 0.03)
    np.testing.assert_array_almost_equal(dense1[3, 0], 0.03)

    # Check dense for attribute 2
    dense2 = f13.dense("surface_directional_effective_roughness_length")
    assert dense2.shape == (4, 3)
    np.testing.assert_array_almost_equal(dense2[0], [0.1, 0.2, 0.3])  # nondefault
    np.testing.assert_array_almost_equal(dense2[1], [0.0, 0.0, 0.0])
    np.testing.assert_array_almost_equal(dense2[2], [0.0, 0.0, 0.0])
    np.testing.assert_array_almost_equal(dense2[3], [0.4, 0.5, 0.6])  # nondefault


def test_fort13_roundtrip_identity(tmp_path):
    """Test that reading -> writing -> reading preserves all data."""
    fixture_path = Path(__file__).parent / "fixtures" / "fort13" / "sample.13"
    f13_original = read_fort13(fixture_path)

    # Write to temp file
    temp_file = tmp_path / "roundtrip.13"
    write_fort13(f13_original, temp_file)

    # Read back
    f13_roundtrip = read_fort13(temp_file)

    # Compare structure
    assert f13_roundtrip.grid_name == f13_original.grid_name
    assert f13_roundtrip.num_nodes == f13_original.num_nodes
    assert len(f13_roundtrip.attributes) == len(f13_original.attributes)

    # Compare each attribute
    for orig_attr, round_attr in zip(f13_original.attributes, f13_roundtrip.attributes):
        assert orig_attr.name == round_attr.name
        assert orig_attr.units == round_attr.units
        assert orig_attr.values_per_node == round_attr.values_per_node
        np.testing.assert_array_equal(orig_attr.default_values, round_attr.default_values)

        # Compare nondefault dicts
        assert set(orig_attr.nondefault.keys()) == set(round_attr.nondefault.keys())
        for node_id in orig_attr.nondefault:
            np.testing.assert_array_almost_equal(
                orig_attr.nondefault[node_id],
                round_attr.nondefault[node_id]
            )


def test_public_exports():
    """Smoke test that all public names are exported from chilmesh."""
    from chilmesh import (
        Fort13,
        NodalAttribute,
        read_fort13 as rf13,
        write_fort13 as wf13,
        Fort13ParseError,
    )
    assert Fort13 is not None
    assert NodalAttribute is not None
    assert rf13 is not None
    assert wf13 is not None
    assert Fort13ParseError is not None


def test_attribute_method():
    """Test the Fort13.attribute(name) method."""
    fixture_path = Path(__file__).parent / "fixtures" / "fort13" / "sample.13"
    f13 = read_fort13(fixture_path)

    # Found attribute
    attr = f13.attribute("primitive_weighting_in_continuity_equation")
    assert attr.name == "primitive_weighting_in_continuity_equation"

    # Missing attribute raises KeyError
    with pytest.raises(KeyError, match="not found"):
        f13.attribute("nonexistent_attribute")
