"""Frozen result types for the mesh-element-validity test suite (spec 007)."""
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

    def categories(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for v in self.violations:
            counts[v.category] = counts.get(v.category, 0) + 1
        return counts
