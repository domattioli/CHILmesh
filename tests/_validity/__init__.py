"""Mesh element validity test-suite helpers (spec 007).

This package is test-suite-internal. It is NOT part of the public chilmesh API.
Promotion to ``chilmesh.validate`` is tracked at issue #142.
"""
from tests._validity.types import (
    InformationalNote,
    MeshValidityReport,
    Violation,
)
from tests._validity.validator import validate_mesh_elements

__all__ = [
    "InformationalNote",
    "MeshValidityReport",
    "Violation",
    "validate_mesh_elements",
]
