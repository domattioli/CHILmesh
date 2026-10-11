"""Live parity for bundled meshes intentionally shared with Valence."""
import hashlib
import os
from importlib import resources

import pytest

from chilmesh import registry


@pytest.mark.skipif(not os.environ.get("VALENCE_TOKEN"), reason="VALENCE_TOKEN unset")
@pytest.mark.parametrize("name", ["Block_O.14", "structuredMesh1.14"])
def test_bundled_bytes_match_registry_hash(name):
    data = resources.files("chilmesh").joinpath("data", name).read_bytes()
    assert hashlib.sha256(data).hexdigest() == registry.fetch_hash(name)
