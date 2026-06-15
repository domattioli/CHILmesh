"""MATLAB-parity tests for skeletonization on external (Valence) meshes.

This is a sibling of ``tests/test_skeletonization_matlab_parity.py`` that pins
expected layer counts for meshes from the external Valence catalog,
rather than the bundled fixtures.

Why a separate file? Bundled-fixture parity is a fast, always-on guardrail
(every CI push runs it). External-mesh parity requires installing the
``valence-domains`` package and downloading mesh files, so it is opt-in. Keeping
the two concerns in separate files lets the fast tests stay cheap while still
documenting the broader correctness expectation for the maintainer.

The Python ``_skeletonize()`` method (see ``src/chilmesh/CHILmesh.py``) is a
faithful port of the original MATLAB ``meshLayers`` function from QuADMesh+
(``domattioli/QuADMesh-MATLAB/blob/master/00_CHILMesh_Class/@CHILmesh/CHILmesh.m``).
For each mesh below, the Python port must produce the same ``n_layers`` value
as the MATLAB reference.

Reference values were captured from the original QuADMesh+ ``meshLayers``
algorithm in ``00_CHILMesh_Class/@CHILmesh/CHILmesh.m``. The seven values tagged
"F12 captured 2026-05-23" were produced by running that MATLAB class under
GNU Octave 8.4 on the Valence catalog meshes (connectivity + points fed
to the 2-arg ``CHILmesh(ConnectivityList, Points)`` constructor, bypassing the
MATLAB ``readFort14`` reader). The harness was validated first against the three
already-known references — delaware-bay@default (17), lake-erie@5k (17),
misc-tests@wetting-drying (15) — which it reproduced exactly before the unknowns
were captured. The Python port's ``n_layers`` agreed with the Octave oracle on
all ten meshes. To run these tests locally:

.. code-block:: bash

   pip install valence-domains
   CHILMESH_RUN_EXTERNAL_PARITY=1 python -m pytest \\
       tests/test_skeletonization_matlab_parity_external.py -v

Contributing additional reference data
--------------------------------------
If you have a MATLAB QuADMesh+ install and a mesh that is not yet listed,
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
# Each entry is keyed by Valence identifier (catalog mesh name + variant tag).
# Values:
#   - int: exact n_layers known
#   - (int, int) tuple: range (lo, hi) when source mesh variant is unconfirmed
#   - None: known to be in the catalog but no MATLAB count captured yet
MATLAB_REFERENCE_LAYER_COUNTS: dict[str, Optional[int | Tuple[int, int]]] = {
    "italy@default-v1": 15,
    "lake-erie@5k-nodes-v1": 17,
    "lake-erie@refined-v1": 20,                  # F12 captured 2026-05-23 (Lake_Erie_mesh_refined.14)
    "delaware-bay@default-v1": 17,
    "delaware-bay@refined-100-20000-v1": 17,     # F12 captured 2026-05-23 (Deleware_Bay_hmin_100_hmax_20000.14)
    "wnat@hagen-v1": (35, 45),    # WNAT family observed at ~39; exact variant unconfirmed
    "wnat@onur-v1": (35, 45),
    "wnat@test-v1": (35, 45),
    "wnat@nc-inundation-v6c-v1": (35, 45),
    "misc-tests@wetting-drying-v1": 15,
    "chesapeake-bay@default-v1": 55,             # F12 captured 2026-05-23 (Chesapeake_Bay.14)
    "great-lakes@default-v1": 46,                # F12 captured 2026-05-23 (Great_Lakes.14)
    "lake-michigan@default-v1": 25,              # F12 captured 2026-05-23 (Lake_Michigan_mesh.14)
    "baranja-hill@default-v1": 12,               # F12 captured 2026-05-23 (Baranja_Hill.14)
    "baranja-hill@admesh-v2-v1": 10,             # F12 captured 2026-05-23 (Baranja_Hill_ADMESH_v2.14)
}


def _admesh_domains_available() -> bool:
    try:
        import admesh_domains  # noqa: F401
        return True
    except ImportError:
        return False


def _load_mesh(catalog_id: str):
    """Resolve an Valence catalog ID to a loaded CHILmesh.

    Catalog ID convention: ``"<mesh-name>@<variant>-v<version>"``, e.g.
    ``"italy@default-v1"``. The exact resolution depends on the
    ``valence-domains`` Python loader; we use a placeholder that will be wired
    up once the loader API is finalized.
    """
    import admesh_domains  # type: ignore

    mesh_name, variant_with_version = catalog_id.split("@", 1)
    variant, _, version = variant_with_version.rpartition("-v")

    record = admesh_domains.get(mesh_name, variant=variant, version=int(version))
    record.load()  # Valence lazy-load contract

    from chilmesh import CHILmesh
    return CHILmesh.from_admesh_domain(record, compute_layers=True)


_RUN_EXTERNAL = bool(os.environ.get("CHILMESH_RUN_EXTERNAL_PARITY"))
_SKIP_REASON_EXTERNAL = (
    "Set CHILMESH_RUN_EXTERNAL_PARITY=1 to run the external MATLAB-parity tests "
    "against the Valence catalog. Requires `pip install valence-domains`."
)


@pytest.mark.skipif(not _RUN_EXTERNAL, reason=_SKIP_REASON_EXTERNAL)
@pytest.mark.skipif(
    not _admesh_domains_available(),
    reason="valence-domains package not installed; `pip install valence-domains`",
)
@pytest.mark.parametrize(
    "catalog_id,expected",
    [(cid, exp) for cid, exp in MATLAB_REFERENCE_LAYER_COUNTS.items() if exp is not None],
)
def test_layer_count_matches_matlab_reference(catalog_id: str, expected) -> None:
    """For each Valence mesh with a known MATLAB ``n_layers``,
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


@pytest.mark.skipif(not _RUN_EXTERNAL, reason=_SKIP_REASON_EXTERNAL)
@pytest.mark.skipif(
    not _admesh_domains_available(),
    reason="valence-domains package not installed; `pip install valence-domains`",
)
def test_uncovered_meshes_have_fixme() -> None:
    """Meta-test: every Valence entry without a captured MATLAB count
    must remain in the table with ``None`` so it is visible as a FIXME during
    review.
    """
    uncovered = [k for k, v in MATLAB_REFERENCE_LAYER_COUNTS.items() if v is None]
    if uncovered:
        pytest.skip(
            f"{len(uncovered)} Valence meshes lack a MATLAB reference layer "
            f"count. Capture via QuADMesh+ and add to MATLAB_REFERENCE_LAYER_COUNTS: "
            f"{', '.join(uncovered)}"
        )
