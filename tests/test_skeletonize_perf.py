"""Perf-regression guard for the pure-Python skeletonize path (#202).

#202 reported a >200s pure-Python full-init cliff on Block_O (~5k elems).
Re-measurement on current ``development`` shows that cliff is gone: Block_O
full init (read + adjacency + skeletonize + spatial index) is ~0.2-0.3s and
scaling is linear (~1s per 60k elems). This test pins that — a regression
back toward the cliff trips the generous 30s ceiling well before it becomes
a >200s hang. The ceiling is a tripwire, not a target; it is deliberately
loose to avoid flakiness on slow / loaded CI runners.
"""
from __future__ import annotations

import time
from importlib.resources import files

from chilmesh import CHILmesh

# Generous regression tripwire. Observed ~0.2-0.3s; >200s was the original
# #202 cliff. 30s catches a real regression long before a hang.
BLOCK_O_INIT_CEILING_S = 30.0


def test_block_o_pure_python_init_under_ceiling():
    """Block_O full init must stay well under the historical #202 cliff."""
    path = files("chilmesh.data") / "Block_O.14"
    start = time.perf_counter()
    mesh = CHILmesh.read_from_fort14(str(path))
    elapsed = time.perf_counter() - start

    # Sanity: this is the real Block_O mesh and layers were computed.
    assert mesh.n_elems > 5000
    assert mesh.n_layers > 0
    assert elapsed < BLOCK_O_INIT_CEILING_S, (
        f"Block_O pure-Python full init took {elapsed:.1f}s "
        f"(ceiling {BLOCK_O_INIT_CEILING_S}s) — possible #202 perf regression."
    )
