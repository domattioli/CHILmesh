"""Tests for partition_mesh (spec 005 US3, FAIL branch of T016)."""
import builtins

import numpy as np
import pytest

import chilmesh
from chilmesh import partition_mesh, partition_report

MESHES = ["donut", "structured", "block_o"]


@pytest.mark.parametrize("k", [4, 8])
@pytest.mark.parametrize("name", MESHES)
def test_spectral_threshold(name, k):  # SC-004
    mesh = getattr(chilmesh.examples, name)()
    labels = partition_mesh(mesh, k, method="spectral")
    assert labels.shape == (mesh.n_elems,)
    assert set(np.unique(labels)) == set(range(k))
    rep = partition_report(mesh, labels)
    assert rep.imbalance <= 1.25
    assert rep.pieces.max() == 1


def test_non_power_of_two():
    mesh = chilmesh.examples.structured()
    labels = partition_mesh(mesh, 5)
    rep = partition_report(mesh, labels)
    assert rep.imbalance <= 1.25 and rep.pieces.max() == 1


def test_deterministic():
    mesh = chilmesh.examples.donut()
    assert np.array_equal(partition_mesh(mesh, 4), partition_mesh(mesh, 4))


@pytest.mark.parametrize("k", [0, -1, 10**6])
def test_k_range(k):
    with pytest.raises(ValueError):
        partition_mesh(chilmesh.examples.donut(), k)


def test_unknown_method():
    with pytest.raises(ValueError):
        partition_mesh(chilmesh.examples.donut(), 4, method="bogus")


def test_metis_missing(monkeypatch):
    real = builtins.__import__

    def fake(name, *a, **kw):
        if name == "pymetis":
            raise ImportError("no pymetis")
        return real(name, *a, **kw)

    monkeypatch.setattr(builtins, "__import__", fake)
    with pytest.raises(ImportError, match="partition"):
        partition_mesh(chilmesh.examples.donut(), 4, method="metis")
