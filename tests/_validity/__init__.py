"""Mesh element validity test-suite helpers.

Promoted to ``chilmesh.validation`` in v1.0.0. This module re-exports
from the canonical location for backward compatibility.
"""
from chilmesh.validation import (
    InformationalNote,
    MeshValidityReport,
    Violation,
    validate_mesh_elements,
)

__all__ = [
    "InformationalNote",
    "MeshValidityReport",
    "Violation",
    "validate_mesh_elements",
]
