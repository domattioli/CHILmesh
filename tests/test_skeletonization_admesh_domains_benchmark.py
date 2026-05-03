"""Benchmark: MATLAB-parity validation for skeletonization against ADMESH-Domains meshes.

The Python ``_skeletonize()`` method (see ``src/chilmesh/CHILmesh.py``) is a port of
the original MATLAB ``meshLayers`` function from the QuADMesh+ codebase
(``domattioli/QuADMesh-MATLAB/blob/master/00_CHILMesh_Class/@CHILmesh/CHILmesh.m``).

Faithfulness of the port is established by two complementary checks:

1. **Algorithmic correctness** (``test_skeletonization_invariant.py``): the medial-axis
   layer separation invariant — vertices in layer ``k`` cannot appear in elements of
   layer ``k+2`` or beyond — must hold across all bundled fixtures. This is verified
   in CI on every push.

2. **Numerical parity with MATLAB** (this file): for real-world hydrodynamic domains
   from the ADMESH-Domains catalog, the Python implementation must produce the same
   ``n_layers`` value as the original MATLAB ``meshLayers`` function. These
   reference values were captured by the project maintainer from external MATLAB
   runs and are documented below.

This module is **skipped by default** because the meshes are not bundled with
CHILmesh — they live in the external ``ADMESH-Domains`` catalog. To run the
benchmark, install the catalog and set the ``CHILMESH_RUN_DOMAINS_BENCHMARK``
environment variable:

.. code-block:: bash

   pip install admesh-domains
   CHILMESH_RUN_DOMAINS_BENCHMARK=1 python -m pytest tests/test_skeletonization_admesh_domains_benchmark.py -v

Contributing additional reference data
--------------------------------------
If you have a MATLAB QuADMesh+ install and a mesh that is not yet listed below,
run::

    addpath('00_CHILMesh_Class');
    cm = CHILmesh(...your fixture...);
    cm = cm.meshLayers();
    fprintf('n_layers = %d\\n', cm.nLayers);

then add the result to the ``MATLAB_REFERENCE_LAYER_COUNTS`` table below with a
short comment naming the variant and the date the reference was captured.
"""
from __future__ import annotations

import os
from typing import Optional, Tuple

import pytest


# Maintainer-provided reference layer counts from external MATLAB ``meshLayers`` runs.
# Each entry is keyed by ADMESH-Domains identifier (catalog mesh name + variant tag).
# Values:
#   - int: exact n_layers known
#   - (int, int) tuple: range (lo, hi) when source mesh variant is unconfirmed
#   - None: known to be in the catalog but no MATLAB count captured yet
MATLAB_REFERENCE_LAYER_COUNTS: dict[str, Optional[int | Tuple[int, int]]] = {
    "italy@default-v1": 15,
    "lake-erie@5k-nodes-v1": 17,
    "lake-erie@refined-v1": None,
    "delaware-bay@default-v1": 17,
    "delaware-bay@refined-100-20000-v1": None,
    "wnat@hagen-v1": (35, 45),    # WNAT family observed at ~39; exact variant unconfirmed
    "wnat@onur-v1": (35, 45),
    "wnat@test-v1": (35, 45),
    "wnat@nc-inundation-v6c-v1": (35, 45),
    "misc-tests@wetting-drying-v1": 15,
    "chesapeake-bay@default-v1": None,
    "great-lakes@default-v1": None,
    "lake-michigan@default-v1": None,
    "baranja-hill@default-v1": None,
    "baranja-hill@admesh-v2-v1": None,
}


def _admesh_domains_available() -> bool:
    try:
        import admesh_domains  # noqa: F401
        return True
    except ImportError:
        return False


def _load_mesh(catalog_id: str):
    """Resolve an ADMESH-Domains catalog ID to a loaded CHILmesh.

    Catalog ID convention: ``"<mesh-name>@<variant>-v<version>"``, e.g.
    ``"italy@default-v1"``. The exact resolution depends on the
    ``admesh-domains`` Python loader; we use a placeholder that will be wired
    up once the loader API is finalized.
    """
    import admesh_domains  # type: ignore

    mesh_name, variant_with_version = catalog_id.split("@", 1)
    variant, _, version = variant_with_version.rpartition("-v")

    record = admesh_domains.get(mesh_name, variant=variant, version=int(version))
    record.load()  # ADMESH-Domains lazy-load contract

    from chilmesh import CHILmesh
    return CHILmesh.from_admesh_domain(record, compute_layers=True)


_RUN_BENCHMARK = bool(os.environ.get("CHILMESH_RUN_DOMAINS_BENCHMARK"))
_SKIP_REASON_BENCHMARK = (
    "Set CHILMESH_RUN_DOMAINS_BENCHMARK=1 to run the ADMESH-Domains MATLAB-parity "
    "benchmark. Requires `pip install admesh-domains` and network access to fetch meshes."
)


@pytest.mark.skipif(not _RUN_BENCHMARK, reason=_SKIP_REASON_BENCHMARK)
@pytest.mark.skipif(
    not _admesh_domains_available(),
    reason="admesh-domains package not installed; `pip install admesh-domains`",
)
@pytest.mark.parametrize(
    "catalog_id,expected",
    [(cid, exp) for cid, exp in MATLAB_REFERENCE_LAYER_COUNTS.items() if exp is not None],
)
def test_layer_count_matches_matlab(catalog_id: str, expected) -> None:
    """For each ADMESH-Domains mesh with a known MATLAB ``n_layers``,
    the Python port must produce the same value (or fall in the documented range).
    """
    mesh = _load_mesh(catalog_id)

    if isinstance(expected, tuple):
        lo, hi = expected
        assert lo <= mesh.n_layers <= hi, (
            f"{catalog_id}: expected MATLAB n_layers in [{lo}, {hi}], "
            f"got {mesh.n_layers} from Python port"
        )
    else:
        assert mesh.n_layers == expected, (
            f"{catalog_id}: expected MATLAB n_layers == {expected}, "
            f"got {mesh.n_layers} from Python port"
        )


@pytest.mark.skipif(not _RUN_BENCHMARK, reason=_SKIP_REASON_BENCHMARK)
@pytest.mark.skipif(
    not _admesh_domains_available(),
    reason="admesh-domains package not installed; `pip install admesh-domains`",
)
def test_uncovered_meshes_have_fixme() -> None:
    """Meta-test: every ADMESH-Domains entry without a captured MATLAB count
    must remain in the table with ``None`` so it is visible as a FIXME during
    benchmark review.
    """
    uncovered = [k for k, v in MATLAB_REFERENCE_LAYER_COUNTS.items() if v is None]
    if uncovered:
        pytest.skip(
            f"{len(uncovered)} ADMESH-Domains meshes lack a MATLAB reference layer "
            f"count. Capture via QuADMesh+ and add to MATLAB_REFERENCE_LAYER_COUNTS: "
            f"{', '.join(uncovered)}"
        )
