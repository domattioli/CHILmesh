"""ADCIRC fort.14 boundary-type parsing + round-trip (#129).

CHILmesh historically read only the node + element blocks of a fort.14 file and
discarded the trailing NOPE/NBOU boundary section, so open-ocean vs land
boundary types never survived the I/O pipeline. These tests pin the foundation
for type-aware skeletonization: boundary segments are parsed, exposed on
``mesh.boundary_segments``, and written back losslessly.
"""
from __future__ import annotations

import warnings

import numpy as np

import chilmesh
from chilmesh.CHILmesh import CHILmesh

# A 4-node unit square split into two triangles, with one open-ocean boundary
# segment (bottom edge, nodes 1-2) and one land boundary segment (top edge,
# IBTYPE 0, nodes 3-4).
TYPED_FORT14 = """test mesh
2 4
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 1.0 1.0 0.0
4 0.0 1.0 0.0
1 3 1 2 3
2 3 1 3 4
1
2
2
1
2
1
2
2 0
3
4
"""


def _write(tmp_path, text, name="typed.fort.14"):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def test_typed_boundaries_parsed(tmp_path):
    path = _write(tmp_path, TYPED_FORT14)
    mesh = CHILmesh.read_from_fort14(path, compute_layers=False, compute_adjacencies=False)

    assert len(mesh.boundary_segments) == 2, "one open + one land segment expected"

    open_seg, land_seg = mesh.boundary_segments
    assert open_seg["kind"] == "open"
    assert open_seg["ibtype"] is None
    np.testing.assert_array_equal(open_seg["nodes"], [0, 1])

    assert land_seg["kind"] == "flow"
    assert land_seg["ibtype"] == 0
    np.testing.assert_array_equal(land_seg["nodes"], [2, 3])


def test_legacy_fort14_has_empty_boundary_segments():
    # Bundled fixtures carry no boundary section; load must be unchanged and the
    # attribute must default to an empty list rather than raising or being None.
    mesh = chilmesh.examples.annulus()
    assert mesh.boundary_segments == []


def test_roundtrip_preserves_boundary_segments(tmp_path):
    src = _write(tmp_path, TYPED_FORT14)
    mesh = CHILmesh.read_from_fort14(src, compute_layers=False, compute_adjacencies=False)

    dst = tmp_path / "roundtrip.fort.14"
    assert mesh.write_to_fort14(str(dst)) is True

    reloaded = CHILmesh.read_from_fort14(dst, compute_layers=False, compute_adjacencies=False)

    assert len(reloaded.boundary_segments) == len(mesh.boundary_segments)
    for before, after in zip(mesh.boundary_segments, reloaded.boundary_segments):
        assert before["kind"] == after["kind"]
        assert before["ibtype"] == after["ibtype"]
        np.testing.assert_array_equal(before["nodes"], after["nodes"])


def test_write_without_boundaries_omits_section(tmp_path):
    # A mesh with no boundary segments must produce byte-identical legacy output
    # (no trailing NOPE/NBOU block).
    mesh = CHILmesh.read_from_fort14(_write(tmp_path, TYPED_FORT14),
                                     compute_layers=False, compute_adjacencies=False)
    mesh.boundary_segments = []

    dst = tmp_path / "plain.fort.14"
    mesh.write_to_fort14(str(dst))
    text = dst.read_text(encoding="utf-8")

    # Last non-empty line is the final element row, not a boundary count.
    last = [ln for ln in text.splitlines() if ln.strip()][-1]
    assert last.split()[:2] == ["2", "3"], "expected final element row, no boundary section"


def test_malformed_boundary_section_warns_not_raises(tmp_path):
    bad = TYPED_FORT14.replace("1\n2\n2\n1\n2\n", "not-an-int\n", 1)
    path = _write(tmp_path, bad, name="bad.fort.14")

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        mesh = CHILmesh.read_from_fort14(path, compute_layers=False, compute_adjacencies=False)

    assert mesh.boundary_segments == []
    assert any("boundary section" in str(w.message) for w in caught)
