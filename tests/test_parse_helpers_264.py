"""Direct tests for the parse/write helpers extracted for issue #264.

The public readers are covered elsewhere; these tests pin each helper's
contract, including the exact exception type and message on malformed input.
"""
from __future__ import annotations

import sys

import numpy as np
import pytest

from chilmesh import examples, fort13_io, fort14_io, gmsh_io, summary_io
from chilmesh.CHILmesh import (
    _read_fort14_boundary_segments,
    _read_fort14_connectivity,
    _read_fort14_points,
)
from chilmesh.fort13_io import Fort13ParseError, NodalAttribute
from chilmesh.fort14_io import Fort14ParseError, _Cursor
from chilmesh.gmsh_io import GmshParseError
from chilmesh.mutations import MutableMesh
from chilmesh.summary_io import SummaryError


# --------------------------------------------------------------------------
# gmsh_io
# --------------------------------------------------------------------------

def test_gmsh_find_sections_last_wins_and_missing():
    assert gmsh_io._find_node_and_element_sections(
        ["$Nodes", "$Elements", "$Nodes"]) == (2, 1)
    with pytest.raises(GmshParseError, match=r"Missing \$Nodes section"):
        gmsh_io._find_node_and_element_sections(["$Elements"])
    with pytest.raises(GmshParseError, match=r"Missing \$Elements section"):
        gmsh_io._find_node_and_element_sections(["$Nodes"])


def test_gmsh_v2_count():
    assert gmsh_io._parse_v2_count(["$Nodes", "7"], 0, "$Nodes") == 7
    with pytest.raises(GmshParseError, match=r"\$Nodes section incomplete"):
        gmsh_io._parse_v2_count(["$Nodes"], 0, "$Nodes")
    with pytest.raises(GmshParseError, match=r"\$Elements section count is not an integer"):
        gmsh_io._parse_v2_count(["$Elements", "x"], 0, "$Elements")


def test_gmsh_v2_node_line():
    assert gmsh_io._parse_v2_node_line("3 1 2 3") == (3, [1.0, 2.0, 3.0])
    with pytest.raises(GmshParseError, match="Node line malformed: 1 2"):
        gmsh_io._parse_v2_node_line("1 2")
    with pytest.raises(GmshParseError, match="non-numeric values: 1 a 2 3"):
        gmsh_io._parse_v2_node_line("1 a 2 3")


def test_gmsh_v2_nodes_truncated():
    lines = ["$Nodes", "2", "1 0 0 0"]
    with pytest.raises(GmshParseError, match="expected 2 nodes"):
        gmsh_io._parse_v2_nodes(lines, 0)
    assert gmsh_io._parse_v2_nodes(lines + ["2 1 0 0"], 0) == {
        1: [0.0, 0.0, 0.0], 2: [1.0, 0.0, 0.0]}


def test_gmsh_v2_element_line():
    assert gmsh_io._parse_v2_element_line("1 2 2 0 0 1 2 3") == ("tri", [1, 2, 3])
    assert gmsh_io._parse_v2_element_line("1 3 0 1 2 3 4") == ("quad", [1, 2, 3, 4])
    assert gmsh_io._parse_v2_element_line("1 15 0 1") is None
    with pytest.raises(GmshParseError, match="Element line malformed: 1 2"):
        gmsh_io._parse_v2_element_line("1 2")
    with pytest.raises(GmshParseError, match="Element line too short"):
        gmsh_io._parse_v2_element_line("1 2 0 1 2")
    with pytest.raises(GmshParseError, match="malformed or non-numeric"):
        gmsh_io._parse_v2_element_line("1 2 0 1 2 x")
    # _require_int raises its own GmshParseError (not a ValueError), so it escapes
    with pytest.raises(GmshParseError, match="malformed element ID: non-numeric 'q'"):
        gmsh_io._parse_v2_element_line("q 2 0 1 2 3")


def test_gmsh_v2_elements_truncated():
    with pytest.raises(GmshParseError, match="expected 1 elements"):
        gmsh_io._parse_v2_elements(["$Elements", "1"], 0)


def test_gmsh_v41_section_header():
    assert gmsh_io._parse_v41_section_header(["$Nodes", "2 4 1 4"], 0, "$Nodes", "number of nodes") == 2
    with pytest.raises(GmshParseError, match=r"\$Nodes section incomplete"):
        gmsh_io._parse_v41_section_header(["$Nodes"], 0, "$Nodes", "number of nodes")
    with pytest.raises(GmshParseError, match=r"\$Elements header malformed"):
        gmsh_io._parse_v41_section_header(["$Elements", "1 2"], 0, "$Elements", "x")
    with pytest.raises(GmshParseError, match="malformed number of nodes: non-numeric 'z'"):
        gmsh_io._parse_v41_section_header(["$Nodes", "1 z 1 1"], 0, "$Nodes", "number of nodes")


def test_gmsh_v41_node_block_header():
    assert gmsh_io._parse_v41_node_block_header(["2 1 0 5"], 0, 0) == 5
    with pytest.raises(GmshParseError, match=r"\$Nodes section incomplete"):
        gmsh_io._parse_v41_node_block_header([], 0, 0)
    with pytest.raises(GmshParseError, match="Node block 1 header malformed"):
        gmsh_io._parse_v41_node_block_header(["2 1"], 0, 1)
    with pytest.raises(GmshParseError, match="Node block 0 header has non-numeric values"):
        gmsh_io._parse_v41_node_block_header(["2 1 0 x"], 0, 0)


def test_gmsh_v41_node_tags_and_coords():
    assert gmsh_io._parse_v41_node_tags(["5", "6"], 0, 0, 2) == [5, 6]
    with pytest.raises(GmshParseError, match="Node block 0 node-tags incomplete"):
        gmsh_io._parse_v41_node_tags(["5"], 0, 0, 2)
    with pytest.raises(GmshParseError, match="Node tag line malformed: q"):
        gmsh_io._parse_v41_node_tags(["q"], 0, 0, 1)

    assert gmsh_io._parse_v41_node_coords(["0 1 2"], 0, 0, [9]) == [(9, [0.0, 1.0, 2.0])]
    with pytest.raises(GmshParseError, match="Node block 0 coordinates incomplete"):
        gmsh_io._parse_v41_node_coords([], 0, 0, [9])
    with pytest.raises(GmshParseError, match="Node coordinate line malformed: 0 1"):
        gmsh_io._parse_v41_node_coords(["0 1"], 0, 0, [9])
    with pytest.raises(GmshParseError, match="non-numeric values: 0 1 z"):
        gmsh_io._parse_v41_node_coords(["0 1 z"], 0, 0, [9])


def test_gmsh_v41_nodes_block():
    lines = ["$Nodes", "1 2 1 2", "2 1 0 2", "1", "2", "0 0 0", "1 0 0"]
    assert gmsh_io._parse_v41_nodes(lines, 0) == {1: [0.0, 0.0, 0.0], 2: [1.0, 0.0, 0.0]}


def test_gmsh_v41_element_block_header_and_line():
    assert gmsh_io._parse_v41_element_block_header(["2 1 2 3"], 0, 0) == (2, 3)
    with pytest.raises(GmshParseError, match=r"\$Elements section incomplete"):
        gmsh_io._parse_v41_element_block_header([], 0, 0)
    with pytest.raises(GmshParseError, match="Element block 0 header malformed"):
        gmsh_io._parse_v41_element_block_header(["1 2"], 0, 0)
    with pytest.raises(GmshParseError, match="Element block 2 header has non-numeric values"):
        gmsh_io._parse_v41_element_block_header(["2 1 x 3"], 0, 2)

    assert gmsh_io._parse_v41_element_line("1 1 2 3", 2) == ("tri", [1, 2, 3])
    assert gmsh_io._parse_v41_element_line("1 1 2 3 4", 3) == ("quad", [1, 2, 3, 4])
    assert gmsh_io._parse_v41_element_line("1 1 2", 1) is None
    with pytest.raises(GmshParseError, match="Triangle element line too short"):
        gmsh_io._parse_v41_element_line("1 1 2", 2)
    with pytest.raises(GmshParseError, match="Quad element line too short"):
        gmsh_io._parse_v41_element_line("1 1 2 3", 3)
    with pytest.raises(GmshParseError, match="Element line malformed: 1 1 2 z"):
        gmsh_io._parse_v41_element_line("1 1 2 z", 2)
    with pytest.raises(GmshParseError, match="Element line malformed:"):
        gmsh_io._parse_v41_element_line("", 2)


def test_gmsh_v41_elements_block_incomplete():
    lines = ["$Elements", "1 2 1 2", "2 1 2 2", "1 1 2 3"]
    with pytest.raises(GmshParseError, match="Element block 0 incomplete"):
        gmsh_io._parse_v41_elements(lines, 0)


# --------------------------------------------------------------------------
# fort14_io
# --------------------------------------------------------------------------

def test_fort14_raw_node_rows():
    nid, coords, i = fort14_io._read_node_rows(["1 0 0 5", "2 1 1"], 0, 2, "p")
    assert nid == [1, 2] and coords[2] == (1.0, 1.0, 0.0) and i == 2
    with pytest.raises(Fort14ParseError, match="p: bad node row near line 1"):
        fort14_io._read_node_rows(["1 0 0 5"], 0, 2, "p")
    with pytest.raises(Fort14ParseError, match="bad node row"):
        fort14_io._read_node_rows(["1 a 0 5"], 0, 1, "p")


def test_fort14_raw_element_rows():
    ids, elems, i = fort14_io._read_element_rows(["1 3 1 2 3"], 0, 1, "p")
    assert ids == [1] and elems == {1: (1, 2, 3)} and i == 1
    with pytest.raises(Fort14ParseError, match="element 1 declares 4 nodes, found 3"):
        fort14_io._read_element_rows(["1 4 1 2 3"], 0, 1, "p")
    with pytest.raises(Fort14ParseError, match="bad element row"):
        fort14_io._read_element_rows([], 0, 1, "p")


def test_fort14_raw_open_and_flow_segments():
    lines = ["2 0 extra", "10", "11"]
    seg = fort14_io._parse_open_segment(lines, _Cursor(0))
    assert seg.nodes == [10, 11] and seg.ibtype == 0
    # truncated file: tolerated, returns the nodes seen
    cur = _Cursor(0)
    assert fort14_io._parse_open_segment(["3", "10"], cur).nodes == [10]
    with pytest.raises(ValueError):
        fort14_io._parse_open_segment(["x"], _Cursor(0))

    flow = fort14_io._parse_flow_segment(["1 24", "5 6 1.5 2.0 3.0"], _Cursor(0))
    assert flow.back_nodes == [6] and flow.heights == [1.5] and flow.coeffs == [[2.0, 3.0]]
    one = fort14_io._parse_flow_segment(["1 3", "5 1.5 2.0"], _Cursor(0))
    assert one.heights == [1.5] and one.coeffs == [[2.0]]
    assert fort14_io._parse_flow_boundaries(["x"], _Cursor(1)) == []


def test_fort14_raw_accumulate_flow_node_malformed():
    with pytest.raises(ValueError):
        fort14_io._accumulate_flow_node(["1", "bad"], 4, [], [], [], [], [])


def test_fort14_raw_boundary_block_error_reports_line():
    with pytest.raises(Fort14ParseError, match="malformed boundary block near line 0"):
        fort14_io._read_boundary_block(["zz"], 0, "p")


# --------------------------------------------------------------------------
# fort13_io
# --------------------------------------------------------------------------

def test_fort13_header():
    assert fort13_io._parse_f13_header(["g", "4", "1"]) == ("g", 4, 1)
    with pytest.raises(Fort13ParseError, match="too short"):
        fort13_io._parse_f13_header(["g", "4"])
    with pytest.raises(Fort13ParseError, match="header parse error"):
        fort13_io._parse_f13_header(["g", "x", "1"])


def test_fort13_attribute_meta():
    attr, nxt = fort13_io._parse_f13_attribute_meta(["a", "u", "2", "0.1 0.2"], 0)
    assert attr.name == "a" and attr.values_per_node == 2 and nxt == 4
    with pytest.raises(Fort13ParseError, match="metadata section incomplete"):
        fort13_io._parse_f13_attribute_meta(["a", "u"], 0)
    with pytest.raises(Fort13ParseError, match="values_per_node parse error at line 2"):
        fort13_io._parse_f13_attribute_meta(["a", "u", "z"], 0)
    with pytest.raises(Fort13ParseError, match="default values line missing"):
        fort13_io._parse_f13_attribute_meta(["a", "u", "1"], 0)
    with pytest.raises(Fort13ParseError, match="expects 2 default values, got 1"):
        fort13_io._parse_f13_attribute_meta(["a", "u", "2", "0.1"], 0)
    with pytest.raises(Fort13ParseError, match="Default values parse error"):
        fort13_io._parse_f13_attribute_meta(["a", "u", "1", "q"], 0)


def _attr(vpn=1):
    return NodalAttribute("a", "u", vpn, np.zeros(vpn))


def test_fort13_data_row_and_block():
    nid, vals = fort13_io._parse_f13_data_row("2 0.5", _attr(), 3)
    assert nid == 1 and vals.tolist() == [0.5]
    with pytest.raises(Fort13ParseError, match="expects 1 \\+ 1 tokens, got 3"):
        fort13_io._parse_f13_data_row("1 2 3", _attr(), 3)
    # range error is re-wrapped with the "Data row parse error" prefix
    with pytest.raises(Fort13ParseError, match=r"Data row parse error: Node ID 9 out of range \[1, 3\]"):
        fort13_io._parse_f13_data_row("9 0.5", _attr(), 3)
    with pytest.raises(Fort13ParseError, match="Data row parse error"):
        fort13_io._parse_f13_data_row("1 zz", _attr(), 3)

    attrs = [_attr()]
    end = fort13_io._parse_f13_data_block(["a", "1", "1 0.7"], 0, attrs, 3)
    assert end == 3 and attrs[0].nondefault[0].tolist() == [0.7]
    with pytest.raises(Fort13ParseError, match="data section incomplete"):
        fort13_io._parse_f13_data_block([], 0, attrs, 3)
    with pytest.raises(Fort13ParseError, match="Unknown attribute 'b'"):
        fort13_io._parse_f13_data_block(["b"], 0, attrs, 3)
    with pytest.raises(Fort13ParseError, match="num_nondefault line missing for 'a'"):
        fort13_io._parse_f13_data_block(["a"], 0, attrs, 3)
    with pytest.raises(Fort13ParseError, match="num_nondefault parse error"):
        fort13_io._parse_f13_data_block(["a", "x"], 0, attrs, 3)
    with pytest.raises(Fort13ParseError, match="data row missing for 'a'"):
        fort13_io._parse_f13_data_block(["a", "1"], 0, attrs, 3)


# --------------------------------------------------------------------------
# CHILmesh.read_from_fort14 / write_fort14 helpers
# --------------------------------------------------------------------------

def test_fort14_points_and_connectivity_helpers():
    pts, i = _read_fort14_points(["1 0 0 1", "2 1 0"], 0, 2)
    assert pts.tolist() == [[0, 0, 1], [1, 0, 0]] and i == 2
    with pytest.raises(IndexError):
        _read_fort14_points(["1 0"], 0, 1)
    with pytest.raises(ValueError):
        _read_fort14_points(["1 a 0"], 0, 1)

    conn, i = _read_fort14_connectivity(["1 3 1 2 3", "2 4 1 2 3 4"], 0, 2)
    assert conn.tolist() == [[0, 1, 2, 0], [0, 1, 2, 3]] and i == 2
    tri, _ = _read_fort14_connectivity(["1 3 1 2 3"], 0, 1)
    assert tri.shape == (1, 3)
    with pytest.raises(ValueError):
        _read_fort14_connectivity(["1 x 1 2 3"], 0, 1)


def test_fort14_boundary_segments_helper():
    lines = ["1", "2", "2", "1", "2", "1", "2", "1 0", "3", "4"]
    segs, present = _read_fort14_boundary_segments(lines, 0)
    assert present and [s["kind"] for s in segs] == ["open", "flow"]
    assert segs[0]["nodes"].tolist() == [0, 1]
    # absent block: IndexError swallowed
    assert _read_fort14_boundary_segments([], 0) == ([], False)
    # malformed row: warns, keeps what was read, still reports present
    with pytest.warns(UserWarning, match="Malformed boundary section"):
        segs, present = _read_fort14_boundary_segments(["1", "1", "x"], 0)
    assert present and segs == []


def test_write_fort14_boundaries_helper(tmp_path):
    class _M:
        boundary_segments = [
            {"kind": "open", "ibtype": None, "nodes": np.array([0, 1])},
            {"kind": "flow", "ibtype": 20, "nodes": np.array([2])},
            {"kind": "flow", "ibtype": None, "nodes": np.array([3])},
        ]

    out = tmp_path / "b.txt"
    with open(out, "w") as f:
        CHILmesh_mod._write_fort14_boundaries(f, _M())
    assert out.read_text() == "1\n2\n2\n1\n2\n2\n2\n1 20\n3\n1\n4\n"

    class _E:
        pass

    with open(out, "w") as f:
        CHILmesh_mod._write_fort14_boundaries(f, _E())
    assert out.read_text() == "0\n0\n0\n0\n"


CHILmesh_mod = sys.modules["chilmesh.CHILmesh"]  # the package re-exports the class under the same name


# --------------------------------------------------------------------------
# summary_io
# --------------------------------------------------------------------------

@pytest.mark.parametrize("name,fmt", [
    ("a.14", "fort14"), ("a.GRD", "fort14"), ("a.fort14", "fort14"),
    ("a.2dm", "2dm"), ("a.13", "fort13"), ("a.15", "fort15"),
    ("a.npy", "npy"), ("a.npz", "npz"), ("a.msh", "gmsh"),
])
def test_detect_format_known(name, fmt):
    from pathlib import Path
    assert summary_io._detect_format(Path(name)) == fmt


def test_detect_format_generic_and_unknown():
    from pathlib import Path
    assert summary_io._detect_format(Path("fort.14")) == "fort14"
    with pytest.raises(SummaryError, match="Unknown mesh format: .txt"):
        summary_io._detect_format(Path("mesh.txt"))


def test_apply_deep_summary_wraps_errors(tmp_path):
    bad = tmp_path / "x.14"
    bad.write_text("nonsense\n")
    with pytest.raises(SummaryError, match="Failed to load mesh for deep summary"):
        summary_io._apply_deep_summary(bad, {})


# --------------------------------------------------------------------------
# mutations helpers
# --------------------------------------------------------------------------

@pytest.fixture
def mm():
    return MutableMesh(examples.annulus())


def test_elems_referencing_vertex_skips_excluded_and_tombstones(mm):
    conn = mm.mesh.connectivity_list
    v = int(conn[0, 0])
    full = mm._elems_referencing_vertex(v, exclude=set())
    assert 0 in full
    assert 0 not in mm._elems_referencing_vertex(v, exclude={0})
    saved = conn[0].copy()
    conn[0] = -1
    try:
        assert 0 not in mm._elems_referencing_vertex(v, exclude=set())
    finally:
        conn[0] = saved


def test_check_collapse_no_inversion_raises(mm):
    conn = mm.mesh.connectivity_list
    row = conn[0]
    # Collapse vertex 1 of a triangle onto vertex 0: area drops to zero.
    with pytest.raises(RuntimeError, match=r"collapse_edge\(7\) would invert element 0; aborted"):
        mm._check_collapse_no_inversion(7, [0], int(row[0]), int(row[1]))
    mm._check_collapse_no_inversion(7, [], 0, 1)


def test_first_layer_containing_and_replay(mm):
    layers = mm.mesh.layers
    first_oe = int(layers["OE"][0][0])
    assert mm._first_layer_containing({first_oe}) == 0
    assert mm._first_layer_containing({-5}) is None
    e2v, e2e = mm._replay_consumed_layers(0)
    assert np.array_equal(e2v, mm.mesh.adjacencies["Edge2Vert"])
    assert np.array_equal(e2e, mm.mesh.adjacencies["Edge2Elem"])
    if mm.mesh.n_layers > 1:
        e2v1, e2e1 = mm._replay_consumed_layers(1)
        assert (e2e1 == -1).sum() > (e2e == -1).sum()


def test_peel_one_layer_stops_when_no_boundary(mm):
    e2v = mm.mesh.adjacencies["Edge2Vert"].copy()
    e2e = np.full_like(mm.mesh.adjacencies["Edge2Elem"], -1)
    new = {k: [] for k in ("OE", "IE", "OV", "IV", "bEdgeIDs")}
    assert mm._peel_one_layer(e2v, e2e, new) is False
    assert all(v == [] for v in new.values())


def test_peel_one_layer_matches_full_peel_first_layer(mm):
    e2v = mm.mesh.adjacencies["Edge2Vert"].copy()
    e2e = mm.mesh.adjacencies["Edge2Elem"].copy()
    new = {k: [] for k in ("OE", "IE", "OV", "IV", "bEdgeIDs")}
    assert mm._peel_one_layer(e2v, e2e, new) is True
    for key in new:
        assert np.array_equal(new[key][0], mm.mesh.layers[key][0])
