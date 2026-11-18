"""fort.14 open-boundary IBTYPEE survives read, write, re-read (#309)."""
from __future__ import annotations

import numpy as np

from chilmesh.CHILmesh import CHILmesh

_NODES_ELEMS = """ibtypee mesh
2 4
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 1.0 1.0 0.0
4 0.0 1.0 0.0
1 3 1 2 3
2 3 1 3 4
"""

# Two open segments: IBTYPEE 0 (with a trailing comment) and IBTYPEE 1.
TWO_OPEN = _NODES_ELEMS + """2
4
2 0   !NVDLL(1),IBTYPEE
1
2
2 1
3
4
0
0
"""

NO_OPEN = _NODES_ELEMS + """0
0
1
2
2 0
3
4
"""


def _read(path):
    return CHILmesh.read_from_fort14(path, compute_layers=False, compute_adjacencies=False)


def _roundtrip(tmp_path, text):
    src = tmp_path / "in.fort.14"
    src.write_text(text, encoding="utf-8")
    first = _read(src)
    out = tmp_path / "out.fort.14"
    first.write_to_fort14(str(out))
    return first, _read(out), out.read_text(encoding="utf-8")


def _summary(mesh):
    return [(s["kind"], s["ibtype"], s.get("ibtypee"), s["nodes"].tolist())
            for s in mesh.boundary_segments]


def test_open_ibtypee_roundtrip(tmp_path):
    first, second, written = _roundtrip(tmp_path, TWO_OPEN)
    assert [s.get("ibtypee") for s in first.boundary_segments] == [0, 1]
    assert _summary(second) == _summary(first)
    assert "2 0\n" in written and "2 1\n" in written
    # Flow-style ``ibtype`` stays None on open segments (seed_ibtypes semantics).
    assert all(s["ibtype"] is None for s in first.boundary_segments)


def test_open_header_without_ibtypee_writes_bare_count(tmp_path):
    text = TWO_OPEN.replace("2 0   !NVDLL(1),IBTYPEE", "2").replace("2 1\n3", "2\n3")
    first, second, written = _roundtrip(tmp_path, text)
    assert all("ibtypee" not in s for s in first.boundary_segments)
    assert _summary(second) == _summary(first)
    lines = written.splitlines()
    assert lines[10] == "2" and lines[13] == "2"


def test_mesh_without_open_boundaries_unchanged(tmp_path):
    first, second, written = _roundtrip(tmp_path, NO_OPEN)
    assert _summary(first) == [("flow", 0, None, [2, 3])]
    assert _summary(second) == _summary(first)
    np.testing.assert_allclose(second.points, first.points)
    tail = written.splitlines()[8:]
    assert tail == ["0", "0", "1", "2", "2 0", "3", "4"]
