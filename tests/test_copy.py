"""Deep-copy invariants for :meth:`CHILmesh.copy`.

The test-suite fixture cache (``tests/conftest.py``) relies on ``.copy()`` to
hand out mutation-safe meshes. If ``copy()`` ever degrades to a shallow copy,
every mutating test would silently corrupt the cached mesh and bleed state
across parametrised runs. These tests guard that contract directly.

Tracks .planning/TEST-AUDIT.md F5.
"""
from __future__ import annotations

import numpy as np
import pytest


def test_copy_preserves_scalar_counts(mesh):
    clone = mesh.copy()
    assert clone.n_verts == mesh.n_verts
    assert clone.n_elems == mesh.n_elems
    assert clone.n_layers == mesh.n_layers
    assert clone.type == mesh.type


def test_copy_returns_distinct_instance(mesh):
    clone = mesh.copy()
    assert clone is not mesh


def test_copy_points_independent(fresh_mesh):
    clone = fresh_mesh.copy()
    assert clone.points is not fresh_mesh.points, "points array must not be aliased"
    np.testing.assert_array_equal(clone.points, fresh_mesh.points)

    snapshot = fresh_mesh.points.copy()
    clone.points[0, 0] += 12345.0
    np.testing.assert_array_equal(
        fresh_mesh.points, snapshot,
        err_msg="mutating clone.points leaked into the original",
    )


def test_copy_connectivity_independent(fresh_mesh):
    clone = fresh_mesh.copy()
    assert clone.connectivity_list is not fresh_mesh.connectivity_list
    np.testing.assert_array_equal(clone.connectivity_list, fresh_mesh.connectivity_list)

    snapshot = fresh_mesh.connectivity_list.copy()
    clone.connectivity_list[0, 0] = (clone.connectivity_list[0, 0] + 1) % fresh_mesh.n_verts
    np.testing.assert_array_equal(
        fresh_mesh.connectivity_list, snapshot,
        err_msg="mutating clone.connectivity_list leaked into the original",
    )


def test_copy_adjacencies_dict_not_aliased(mesh):
    clone = mesh.copy()
    assert clone.adjacencies is not mesh.adjacencies, "adjacencies dict must not be aliased"
    assert set(clone.adjacencies.keys()) == set(mesh.adjacencies.keys())


@pytest.mark.parametrize("key", ["Elem2Vert", "Edge2Vert", "Elem2Edge", "Edge2Elem"])
def test_copy_adjacency_arrays_not_aliased(mesh, key):
    clone = mesh.copy()
    arr_orig = mesh.adjacencies[key]
    arr_clone = clone.adjacencies[key]
    assert arr_clone is not arr_orig, f"adjacencies[{key!r}] is aliased"
    np.testing.assert_array_equal(arr_clone, arr_orig)


@pytest.mark.parametrize("key", ["Vert2Edge", "Vert2Elem"])
def test_copy_sparse_adjacencies_not_aliased(fresh_mesh, key):
    clone = fresh_mesh.copy()
    sparse_orig = fresh_mesh.adjacencies[key]
    sparse_clone = clone.adjacencies[key]
    assert sparse_clone is not sparse_orig, f"adjacencies[{key!r}] is aliased"
    assert sparse_clone == sparse_orig

    # Mutating an inner container on the clone must not propagate.
    probe_vert = next(iter(sparse_orig))
    orig_value = set(sparse_orig[probe_vert])
    sparse_clone[probe_vert].clear()
    assert set(sparse_orig[probe_vert]) == orig_value, (
        f"clearing clone.adjacencies[{key!r}][{probe_vert}] leaked into the original"
    )


def test_copy_layers_independent(mesh):
    clone = mesh.copy()
    assert clone.layers is not mesh.layers
    assert set(clone.layers.keys()) == set(mesh.layers.keys())


def test_copy_idempotent_on_repeat(mesh):
    once = mesh.copy()
    twice = once.copy()
    np.testing.assert_array_equal(twice.points, mesh.points)
    np.testing.assert_array_equal(twice.connectivity_list, mesh.connectivity_list)
    assert twice.n_layers == mesh.n_layers
