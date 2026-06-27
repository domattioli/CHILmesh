"""fort.14 I/O roundtrip: read an ADCIRC mesh, write it back, verify identity.

Run:

    python examples/02_fort14_roundtrip.py

Reads the bundled annulus fixture, writes it to a temp path, reads it
back, and asserts vertex/element counts match. Demonstrates the I/O API.
"""
from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np

import chilmesh
from chilmesh import CHILmesh, write_fort14


def main() -> None:
    src = chilmesh.examples.annulus()
    print(f"loaded annulus: {src.n_verts} verts, {src.n_elems} elems")

    with tempfile.TemporaryDirectory() as td:
        out_path = Path(td) / "annulus_copy.fort.14"

        src.grid_name = "annulus-roundtrip"
        write_fort14(src, out_path)
        print(f"wrote {out_path.name} ({out_path.stat().st_size:,} bytes)")

        loaded = CHILmesh.read_from_fort14(out_path)
        print(f"reloaded: {loaded.n_verts} verts, {loaded.n_elems} elems")

        assert loaded.n_verts == src.n_verts, "vert count drift"
        assert loaded.n_elems == src.n_elems, "elem count drift"
        np.testing.assert_allclose(
            loaded.points[:, :2], src.points[:, :2], atol=1e-8,
            err_msg="coordinate drift across roundtrip",
        )
        print("OK: roundtrip preserves vertex/element counts and coordinates")


if __name__ == "__main__":
    main()
