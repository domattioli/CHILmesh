"""Functional smoke test for a built chilmesh_cpp wheel (CIBW_TEST_COMMAND).

Lives in a file, not an inline ``python -c`` string, because the Windows leg
runs the test command through cmd.exe, where nested quoting of a multi-statement
one-liner is fragile (#256).

Runs ``full_init`` on a 2-triangle unit square and asserts the topology
invariants: 4 verts, 2 elems, 5 unique edges (4 boundary + 1 shared diagonal),
adjacency built, at least one peel layer, edge2vert shape (5, 2).
"""
import numpy as np

import chilmesh_cpp

pts = np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]], dtype=np.float64)
conn = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int32)
m = chilmesh_cpp.full_init(pts, conn)
assert m.n_verts == 4, m.n_verts
assert m.n_elems == 2, m.n_elems
assert m.n_edges == 5, m.n_edges
assert m.adjacency_built, "adjacency not built"
assert m.n_layers >= 1, m.n_layers
assert np.asarray(m.edge2vert).shape == (5, 2), np.asarray(m.edge2vert).shape
print("chilmesh_cpp full_init smoke OK:", m.n_verts, "v", m.n_elems, "e",
      m.n_edges, "edges", m.n_layers, "layers")
