"""Frozen contract for the mesh-element-validity entry point.

This file is the spec-side contract. The runtime implementation lives at
tests/_validity/validator.py and MUST match the signature below verbatim.
Any change here requires a spec.md amendment.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Violation:
    category: str
    element_ids: tuple[int, ...]
    edge_ids: tuple[int, ...] | None
    detail: str


@dataclass(frozen=True)
class InformationalNote:
    category: str
    element_ids: tuple[int, ...]
    detail: str


@dataclass(frozen=True)
class MeshValidityReport:
    ok: bool
    violations: tuple[Violation, ...]
    notes: tuple[InformationalNote, ...]
    n_elems_checked: int
    runtime_s: float


def validate_mesh_elements(
    mesh,  # CHILmesh
    *,
    tol: float | None = None,
) -> MeshValidityReport:
    """Verify mesh element validity per spec 007.

    Args:
        mesh: CHILmesh instance.
        tol: Override for the absolute tolerance. If None, uses
            1e-12 * bbox_diag (floored at 1e-15) per FR-013.

    Returns:
        MeshValidityReport. ok is True iff no violations were found.
    """
    raise NotImplementedError("Contract only; see tests/_validity/validator.py")
