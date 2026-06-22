"""Regression guard: quality aspect-ratio path must not use the NumPy-2.0
deprecated 2-D ``np.cross`` (slated for removal). See development maintenance
fix, 2026-06-18. Mutation-verified: reverting the source to ``np.cross`` makes
the no-deprecation assertions fail under NumPy>=2.0.
"""
from __future__ import annotations

import warnings

import numpy as np

from chilmesh.quality import _triangle_quality


def _right_triangle():
    # legs 3 and 4, hypotenuse 5; area = 6.
    return np.array([0.0, 0.0]), np.array([3.0, 0.0]), np.array([0.0, 4.0])


def test_aspect_ratio_value_unchanged():
    v0, v1, v2 = _right_triangle()
    q = _triangle_quality(v0, v1, v2, metric="aspect_ratio")
    # aspect_ratio = 2 * r_in / r_circ; for a 3-4-5 triangle this is a fixed
    # value. Recompute independently from inradius/circumradius.
    la = np.linalg.norm(v1 - v0)
    lb = np.linalg.norm(v2 - v1)
    lc = np.linalg.norm(v0 - v2)
    s = (la + lb + lc) / 2
    area = 6.0
    r_in = area / s
    r_circ = (la * lb * lc) / (4 * area)
    expected = 2 * r_in / r_circ
    assert abs(q - expected) < 1e-12


def test_aspect_ratio_emits_no_numpy_cross_deprecation():
    v0, v1, v2 = _right_triangle()
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        _triangle_quality(v0, v1, v2, metric="aspect_ratio")


def test_chilmesh_element_quality_emits_no_numpy_cross_deprecation():
    # Guards the second fixed site: CHILmesh.element_quality (CHILmesh.py).
    from chilmesh.examples import annulus

    mesh = annulus()
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        mesh.element_quality(metric="aspect_ratio")
