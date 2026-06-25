"""Comprehensive error handling tests for fort.13 I/O.

Tests all Fort13ParseError branches in src/chilmesh/fort13_io.py
(15 error conditions + 2 positive control tests).
"""
import pytest
import numpy as np
from pathlib import Path

from chilmesh import read_fort13, write_fort13, Fort13, NodalAttribute
from chilmesh.fort13_io import Fort13ParseError


def _write(tmp_path, text):
    """Helper: write text to a temporary fort.13 file."""
    p = tmp_path / "test.13"
    p.write_text(text)
    return p


# ============================================================================
# Error Branch Tests (1-15)
# ============================================================================

def test_file_too_short(tmp_path):
    """L77: File with fewer than 3 non-blank lines."""
    text = "Grid Name\n42\n"  # Only 2 non-blank lines
    p = _write(tmp_path, text)
    with pytest.raises(Fort13ParseError, match="too short"):
        read_fort13(p)


def test_header_parse_non_integer_num_nodes(tmp_path):
    """L84-85: num_nodes is not an integer."""
    text = "Grid Name\nnotanint\n2\n"
    p = _write(tmp_path, text)
    with pytest.raises(Fort13ParseError, match="header parse error"):
        read_fort13(p)


def test_header_parse_non_integer_num_attrs(tmp_path):
    """L84-85: num_attrs is not an integer."""
    text = "Grid Name\n10\nnotanint\n"
    p = _write(tmp_path, text)
    with pytest.raises(Fort13ParseError, match="header parse error"):
        read_fort13(p)


def test_metadata_section_incomplete(tmp_path):
    """L94: File declares 1 attribute but ends before its 3 metadata lines."""
    text = "Grid Name\n5\n1\nManningN\n"  # Only attr name, missing units & vpn
    p = _write(tmp_path, text)
    with pytest.raises(Fort13ParseError, match="metadata section incomplete"):
        read_fort13(p)


def test_values_per_node_not_int(tmp_path):
    """L100-101: values_per_node line is not an integer."""
    text = "Grid Name\n5\n1\nManningN\nm/s\nnotanint\n"
    p = _write(tmp_path, text)
    with pytest.raises(Fort13ParseError, match="values_per_node parse error"):
        read_fort13(p)


def test_default_values_line_missing(tmp_path):
    """L107: File ends right after vpn line, before default values line."""
    text = "Grid Name\n5\n1\nManningN\nm/s\n1"  # Missing newline & default values
    p = _write(tmp_path, text)
    with pytest.raises(Fort13ParseError, match="default values line missing"):
        read_fort13(p)


def test_default_token_count_mismatch(tmp_path):
    """L111: Default token count != vpn."""
    text = "Grid Name\n5\n1\nManningN\nm/s\n2\n0.025\n"  # vpn=2 but only 1 default value
    p = _write(tmp_path, text)
    with pytest.raises(Fort13ParseError, match="expects 2 default values, got 1"):
        read_fort13(p)


def test_default_values_not_floats(tmp_path):
    """L117-118: Default values cannot be parsed as floats."""
    text = "Grid Name\n5\n1\nManningN\nm/s\n1\nnotafloat\n"
    p = _write(tmp_path, text)
    with pytest.raises(Fort13ParseError, match="Default values parse error"):
        read_fort13(p)


def test_data_section_incomplete(tmp_path):
    """L134: Valid metadata but data section name line missing."""
    text = "Grid Name\n5\n1\nManningN\nm/s\n1\n0.025\n"  # Metadata OK but no data section
    p = _write(tmp_path, text)
    with pytest.raises(Fort13ParseError, match="data section incomplete"):
        read_fort13(p)


def test_unknown_attribute_in_data_section(tmp_path):
    """L138: Data section references an attribute not in metadata."""
    text = (
        "Grid Name\n5\n1\nManningN\nm/s\n1\n0.025\n"
        "UnknownAttr\n0\n"  # Data section with unknown attr name
    )
    p = _write(tmp_path, text)
    with pytest.raises(Fort13ParseError, match="Unknown attribute"):
        read_fort13(p)


def test_num_nondefault_line_missing(tmp_path):
    """L145: Data section has attr name but file ends before num_nondefault line."""
    text = "Grid Name\n5\n1\nManningN\nm/s\n1\n0.025\nManningN"  # Missing num_nondefault
    p = _write(tmp_path, text)
    with pytest.raises(Fort13ParseError, match="num_nondefault line missing"):
        read_fort13(p)


def test_num_nondefault_not_int(tmp_path):
    """L149-150: num_nondefault line is not an integer."""
    text = (
        "Grid Name\n5\n1\nManningN\nm/s\n1\n0.025\n"
        "ManningN\nnotanint\n"
    )
    p = _write(tmp_path, text)
    with pytest.raises(Fort13ParseError, match="num_nondefault parse error"):
        read_fort13(p)


def test_data_row_missing(tmp_path):
    """L157: Data section declares N nondefault rows but file ends early."""
    text = (
        "Grid Name\n5\n1\nManningN\nm/s\n1\n0.025\n"
        "ManningN\n2\n"  # Declare 2 nondefault rows but provide none
    )
    p = _write(tmp_path, text)
    with pytest.raises(Fort13ParseError, match="data row missing"):
        read_fort13(p)


def test_data_row_token_count_mismatch(tmp_path):
    """L161: Data row has incorrect token count (expected 1 + vpn)."""
    text = (
        "Grid Name\n5\n1\nManningN\nm/s\n2\n0.025 0.030\n"
        "ManningN\n1\n1 0.04\n"  # vpn=2 expects 1 + 2 = 3 tokens, but got 2
    )
    p = _write(tmp_path, text)
    with pytest.raises(Fort13ParseError, match="expects 1 \\+ 2 tokens"):
        read_fort13(p)


def test_node_id_out_of_range(tmp_path):
    """L172: Node ID is outside valid range [1, num_nodes]."""
    text = (
        "Grid Name\n4\n1\nManningN\nm/s\n1\n0.025\n"
        "ManningN\n1\n99 0.04\n"  # num_nodes=4, but node_id=99 is out of range
    )
    p = _write(tmp_path, text)
    with pytest.raises(Fort13ParseError, match="out of range"):
        read_fort13(p)


def test_data_row_values_not_floats(tmp_path):
    """L176-177: Data row values cannot be parsed as floats."""
    text = (
        "Grid Name\n5\n1\nManningN\nm/s\n1\n0.025\n"
        "ManningN\n1\n2 notafloat\n"  # Second token is not a float
    )
    p = _write(tmp_path, text)
    with pytest.raises(Fort13ParseError, match="Data row parse error"):
        read_fort13(p)


# ============================================================================
# Positive Control Tests (16-17)
# ============================================================================

def test_fort13_attribute_missing_raises_keyerror(tmp_path):
    """Positive: Fort13.attribute() raises KeyError for missing attr name."""
    # Build a minimal valid Fort13 via dataclass constructor
    attr = NodalAttribute(
        name="ManningN",
        units="m/s",
        values_per_node=1,
        default_values=np.array([0.025], dtype=np.float64),
        nondefault={}
    )
    f13 = Fort13(grid_name="Test Grid", num_nodes=5, attributes=[attr])

    with pytest.raises(KeyError, match="Attribute 'Missing' not found"):
        f13.attribute("Missing")


def test_valid_fort13_roundtrip(tmp_path):
    """Positive: Valid fort.13 round-trips through read → write → read."""
    # Create a minimal valid fort.13 file
    text = (
        "Test Grid\n"
        "3\n"  # 3 nodes
        "1\n"  # 1 attribute
        "ManningN\n"
        "m/s\n"
        "1\n"  # vpn=1
        "0.025\n"  # default value
        "ManningN\n"
        "1\n"  # 1 nondefault row
        "2 0.030\n"  # node_id=2, value=0.030
    )
    input_path = _write(tmp_path, text)

    # Read fort.13
    f13_read = read_fort13(input_path)
    assert f13_read.grid_name == "Test Grid"
    assert f13_read.num_nodes == 3
    assert len(f13_read.attributes) == 1
    assert f13_read.attributes[0].name == "ManningN"
    assert f13_read.attributes[0].values_per_node == 1
    assert np.allclose(f13_read.attributes[0].default_values, [0.025])
    assert 1 in f13_read.attributes[0].nondefault  # node_id 1 (0-based; 2 from file - 1)
    assert np.allclose(f13_read.attributes[0].nondefault[1], [0.030])

    # Write to a new file
    output_path = tmp_path / "roundtrip.13"
    write_fort13(f13_read, output_path)

    # Read back and verify equality
    f13_reread = read_fort13(output_path)
    assert f13_reread.grid_name == f13_read.grid_name
    assert f13_reread.num_nodes == f13_read.num_nodes
    assert len(f13_reread.attributes) == len(f13_read.attributes)
    for attr_orig, attr_new in zip(f13_read.attributes, f13_reread.attributes):
        assert attr_orig.name == attr_new.name
        assert attr_orig.units == attr_new.units
        assert attr_orig.values_per_node == attr_new.values_per_node
        assert np.allclose(attr_orig.default_values, attr_new.default_values)
        assert set(attr_orig.nondefault.keys()) == set(attr_new.nondefault.keys())
        for nid in attr_orig.nondefault:
            assert np.allclose(attr_orig.nondefault[nid], attr_new.nondefault[nid])
