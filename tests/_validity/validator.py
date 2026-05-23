"""Mesh element validity entry point (spec 007).

Test-suite-internal. Not public chilmesh API. See issue #142.
"""
from __future__ import annotations

import time
from typing import Any

import numpy as np

from chilmesh.validation.broadphase import UniformGridIndex, all_pairs_naive
from chilmesh.validation.predicates import (
    classify_element,
    effective_tol,
    element_bboxes,
    element_vertex_ids,
    is_self_intersecting_quad,
    point_strictly_in_polygon,
    segment_proper_cross,
    signed_area_polygon,
)
from chilmesh.validation.types import (
    InformationalNote,
    MeshValidityReport,
    Violation,
)

_NAIVE_PAIR_THRESHOLD = 5000


def validate_mesh_elements(
    mesh: Any,
    *,
    tol: float | None = None,
) -> MeshValidityReport:
    """Verify mesh element validity per spec 007.

    Args:
        mesh: A CHILmesh instance (or duck-typed equivalent: exposes
            ``connectivity_list``, ``points``, ``n_elems``,
            ``boundary_node_indices()``, ``layers``).
        tol: Override for the absolute geometric tolerance. If ``None``,
            uses ``1e-12 * bbox_diag`` floored at ``1e-15``.

    Returns:
        :class:`MeshValidityReport`. ``ok`` is True iff no violations.
    """
    t0 = time.perf_counter()
    violations: list[Violation] = []
    notes: list[InformationalNote] = []

    n_elems = int(getattr(mesh, "n_elems", 0))
    if n_elems == 0:
        notes.append(InformationalNote("EMPTY_MESH", (), "n_elems == 0"))
        return MeshValidityReport(
            ok=True,
            violations=(),
            notes=tuple(notes),
            n_elems_checked=0,
            runtime_s=time.perf_counter() - t0,
        )

    points = np.asarray(mesh.points, dtype=float)
    if points.shape[1] >= 3:
        z = points[:, 2]
        z_range = float(z.max() - z.min())
    else:
        z_range = 0.0
    points_xy = points[:, :2]

    tol_eff = effective_tol(points_xy, tol)

    if z_range > tol_eff:
        violations.append(
            Violation(
                category="NON_PLANAR_MESH",
                element_ids=(),
                edge_ids=None,
                detail=f"z-range {z_range:.6g} exceeds tol {tol_eff:.6g}",
            )
        )

    connectivity = np.asarray(mesh.connectivity_list)
    if connectivity.ndim != 2:
        violations.append(
            Violation(
                category="UNSUPPORTED_ELEMENT_ARITY",
                element_ids=(),
                edge_ids=None,
                detail=f"connectivity_list has shape {connectivity.shape}",
            )
        )
        return MeshValidityReport(
            ok=False,
            violations=tuple(violations),
            notes=tuple(notes),
            n_elems_checked=n_elems,
            runtime_s=time.perf_counter() - t0,
        )

    boundary_verts: set[int] = set()
    try:
        bv = mesh.boundary_node_indices()
        boundary_verts = {int(x) for x in np.asarray(bv).flatten()}
    except Exception as exc:  # pragma: no cover
        notes.append(
            InformationalNote(
                "BOUNDARY_LOOKUP_FAILED",
                (),
                f"boundary_node_indices() raised {type(exc).__name__}: {exc}",
            )
        )

    layers = getattr(mesh, "layers", None)
    layer0_elems: set[int] = set()
    if layers is None or not layers.get("OE"):
        if hasattr(mesh, "_skeletonize"):
            try:
                mesh._skeletonize()
                notes.append(
                    InformationalNote(
                        "LAYERS_AUTO_TRIGGERED",
                        (),
                        "validator triggered mesh._skeletonize() for FR-007",
                    )
                )
            except Exception as exc:  # pragma: no cover
                notes.append(
                    InformationalNote(
                        "LAYERS_AUTO_TRIGGER_FAILED",
                        (),
                        f"_skeletonize raised {type(exc).__name__}: {exc}",
                    )
                )
        layers = getattr(mesh, "layers", None)

    if layers and layers.get("OE"):
        oe0 = np.asarray(layers["OE"][0]).flatten()
        ie0 = np.asarray(layers["IE"][0]).flatten() if layers.get("IE") else np.array([], dtype=int)
        layer0_elems = {int(x) for x in oe0} | {int(x) for x in ie0}

    elem_types: list[str] = []
    elem_polys_xy: list[np.ndarray] = []
    for elem_id in range(n_elems):
        row = connectivity[elem_id]
        kind = classify_element(row)
        elem_types.append(kind)

        if kind == "UNSUPPORTED_ELEMENT_ARITY":
            violations.append(
                Violation(
                    category="UNSUPPORTED_ELEMENT_ARITY",
                    element_ids=(elem_id,),
                    edge_ids=None,
                    detail=f"row={row.tolist()}",
                )
            )
            elem_polys_xy.append(np.zeros((3, 2)))
            continue

        if kind == "DEGENERATE_QUAD_DUPLICATE_VERTEX":
            notes.append(
                InformationalNote(
                    category="DEGENERATE_QUAD_DUPLICATE_VERTEX",
                    element_ids=(elem_id,),
                    detail=f"row={row.tolist()}",
                )
            )

        verts = element_vertex_ids(row)
        poly = points_xy[list(verts)]
        elem_polys_xy.append(poly)

        if kind == "TRI":
            if not boundary_verts.isdisjoint(verts):
                pass
            else:
                violations.append(
                    Violation(
                        category="INTERIOR_TRIANGLE_FORBIDDEN",
                        element_ids=(elem_id,),
                        edge_ids=None,
                        detail=f"verts={verts} none in boundary_node_indices",
                    )
                )
            if layer0_elems and elem_id not in layer0_elems:
                violations.append(
                    Violation(
                        category="INTERIOR_LAYER_TRIANGLE_FORBIDDEN",
                        element_ids=(elem_id,),
                        edge_ids=None,
                        detail=f"verts={verts} not in layers[OE/IE][0]",
                    )
                )

        if kind in ("QUAD", "DEGENERATE_QUAD_DUPLICATE_VERTEX") and poly.shape[0] == 4:
            bowtie, edge_pair = is_self_intersecting_quad(poly, tol_eff)
            if bowtie:
                violations.append(
                    Violation(
                        category="SELF_INTERSECTING_QUAD",
                        element_ids=(elem_id,),
                        edge_ids=edge_pair,
                        detail=f"verts={verts} coords={poly.tolist()}",
                    )
                )

        area = signed_area_polygon(poly)
        bbox_d = float(
            np.linalg.norm(poly.max(axis=0) - poly.min(axis=0))
        )
        if abs(area) <= max(tol_eff, tol_eff * bbox_d):
            notes.append(
                InformationalNote(
                    category="DEGENERATE_ZERO_AREA",
                    element_ids=(elem_id,),
                    detail=f"signed_area={area:.6g}",
                )
            )

    edge_to_elems: dict[frozenset[int], list[int]] = {}
    for elem_id in range(n_elems):
        verts = element_vertex_ids(connectivity[elem_id])
        if not verts:
            continue
        n = len(verts)
        for i in range(n):
            key = frozenset({verts[i], verts[(i + 1) % n]})
            if len(key) < 2:
                continue
            edge_to_elems.setdefault(key, []).append(elem_id)

    edge_shared_pairs: set[tuple[int, int]] = set()
    for elems in edge_to_elems.values():
        if len(elems) < 2:
            continue
        for i in range(len(elems)):
            for j in range(i + 1, len(elems)):
                a, b = elems[i], elems[j]
                edge_shared_pairs.add((a, b) if a < b else (b, a))

    bboxes = element_bboxes(elem_polys_xy)
    if n_elems <= _NAIVE_PAIR_THRESHOLD:
        pair_iter = all_pairs_naive(n_elems, bboxes)
    else:
        pair_iter = UniformGridIndex(bboxes, n_elems).candidate_pairs()

    for i, j in pair_iter:
        if (i, j) in edge_shared_pairs:
            continue
        poly_i = elem_polys_xy[i]
        poly_j = elem_polys_xy[j]

        verts_i = element_vertex_ids(connectivity[i])
        verts_j = element_vertex_ids(connectivity[j])
        shared_verts = set(verts_i) & set(verts_j)

        crossed = False
        ni = poly_i.shape[0]
        nj = poly_j.shape[0]
        for ei in range(ni):
            a = poly_i[ei]
            b = poly_i[(ei + 1) % ni]
            va_id = verts_i[ei]
            vb_id = verts_i[(ei + 1) % ni]
            for ej in range(nj):
                c = poly_j[ej]
                d = poly_j[(ej + 1) % nj]
                vc_id = verts_j[ej]
                vd_id = verts_j[(ej + 1) % nj]
                if {va_id, vb_id} & {vc_id, vd_id}:
                    continue
                if segment_proper_cross(a, b, c, d, tol_eff):
                    violations.append(
                        Violation(
                            category="EDGE_CROSSING",
                            element_ids=(i, j),
                            edge_ids=(ei, ej),
                            detail=f"elem {i} edge ({va_id},{vb_id}) crosses elem {j} edge ({vc_id},{vd_id})",
                        )
                    )
                    crossed = True
                    break
            if crossed:
                break

        if not crossed:
            interior_hit = False
            for ki in range(ni):
                if verts_i[ki] in shared_verts:
                    continue
                if point_strictly_in_polygon(poly_i[ki], poly_j, tol_eff):
                    violations.append(
                        Violation(
                            category="INTERIOR_OVERLAP",
                            element_ids=(i, j),
                            edge_ids=None,
                            detail=f"vert {verts_i[ki]} of elem {i} lies inside elem {j}",
                        )
                    )
                    interior_hit = True
                    break
            if not interior_hit:
                for kj in range(nj):
                    if verts_j[kj] in shared_verts:
                        continue
                    if point_strictly_in_polygon(poly_j[kj], poly_i, tol_eff):
                        violations.append(
                            Violation(
                                category="INTERIOR_OVERLAP",
                                element_ids=(i, j),
                                edge_ids=None,
                                detail=f"vert {verts_j[kj]} of elem {j} lies inside elem {i}",
                            )
                        )
                        break

    runtime_s = time.perf_counter() - t0
    ok = not any(_is_failure(v.category) for v in violations)
    return MeshValidityReport(
        ok=ok,
        violations=tuple(violations),
        notes=tuple(notes),
        n_elems_checked=n_elems,
        runtime_s=runtime_s,
    )


def _is_failure(category: str) -> bool:
    return category in {
        "UNSUPPORTED_ELEMENT_ARITY",
        "INTERIOR_TRIANGLE_FORBIDDEN",
        "INTERIOR_LAYER_TRIANGLE_FORBIDDEN",
        "NON_PLANAR_MESH",
        "SELF_INTERSECTING_QUAD",
        "INTERIOR_OVERLAP",
        "EDGE_CROSSING",
    }


def format_failures(report: MeshValidityReport, cap_per_category: int = 10) -> str:
    """Build pytest-friendly assertion message aggregating all violation categories."""
    if report.ok:
        return "report.ok == True"
    by_cat: dict[str, list[Violation]] = {}
    for v in report.violations:
        by_cat.setdefault(v.category, []).append(v)
    lines = [
        f"validate_mesh_elements: {len(report.violations)} violations across "
        f"{len(by_cat)} categories on {report.n_elems_checked} elements "
        f"(runtime {report.runtime_s:.3f}s)"
    ]
    for cat, vs in by_cat.items():
        lines.append(f"  [{cat}] x {len(vs)}:")
        for v in vs[:cap_per_category]:
            extras = f" edges={v.edge_ids}" if v.edge_ids else ""
            lines.append(f"    elems={v.element_ids}{extras}: {v.detail}")
        if len(vs) > cap_per_category:
            lines.append(f"    ... and {len(vs) - cap_per_category} more")
    return "\n".join(lines)
