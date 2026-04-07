"""Utilities subpackage for chilmesh.

Re-exports the plotting mixin so that ``from chilmesh.utils import plot_utils``
continues to work and so the wheel always ships ``plot_utils.py`` regardless of
the setuptools auto-discovery rules in effect.
"""
from . import plot_utils  # noqa: F401

__all__ = ["plot_utils"]
