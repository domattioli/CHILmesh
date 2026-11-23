"""Regression: find_element on a 1-element mesh (k=1 KD-tree query returns a bare scalar)."""
import numpy as np

from chilmesh import CHILmesh


def _single_triangle():
    pts = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    conn = np.array([[0, 1, 2]])
    return CHILmesh(connectivity=conn, points=pts, grid_name="single")


def test_point_inside_returns_element_id():
    mesh = _single_triangle()
    assert mesh.find_element(np.array([0.25, 0.25])) == 0


def test_point_outside_returns_minus_one():
    mesh = _single_triangle()
    assert mesh.find_element(np.array([5.0, 5.0])) == -1


def test_multi_element_behavior_unchanged():
    pts = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    conn = np.array([[0, 1, 2], [0, 2, 3]])
    mesh = CHILmesh(connectivity=conn, points=pts, grid_name="two")
    assert mesh.find_element(np.array([0.7, 0.2])) == 0
    assert mesh.find_element(np.array([0.2, 0.7])) == 1
    assert mesh.find_element(np.array([-5.0, 0.0])) == -1
