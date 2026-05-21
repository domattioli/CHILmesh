"""Uniform-grid broadphase for pair enumeration (spec FR-012)."""
from __future__ import annotations

from typing import Iterator

import numpy as np


class UniformGridIndex:
    def __init__(self, elem_bboxes: np.ndarray, n_elems: int):
        self.bboxes = elem_bboxes
        self.n_elems = n_elems

        bbox_min = elem_bboxes[:, :2].min(axis=0)
        bbox_max = elem_bboxes[:, 2:].max(axis=0)
        diag = float(np.linalg.norm(bbox_max - bbox_min))
        if diag == 0.0 or n_elems == 0:
            self.cell_size = 1.0
        else:
            mean_diag = float(
                np.mean(
                    np.linalg.norm(elem_bboxes[:, 2:] - elem_bboxes[:, :2], axis=1)
                )
            )
            self.cell_size = max(mean_diag * 2.0, 1e-12)
        self.bbox_min = bbox_min

        self.cells: dict[tuple[int, int], list[int]] = {}
        for elem_id in range(n_elems):
            i0 = int((elem_bboxes[elem_id, 0] - bbox_min[0]) // self.cell_size)
            j0 = int((elem_bboxes[elem_id, 1] - bbox_min[1]) // self.cell_size)
            i1 = int((elem_bboxes[elem_id, 2] - bbox_min[0]) // self.cell_size)
            j1 = int((elem_bboxes[elem_id, 3] - bbox_min[1]) // self.cell_size)
            for i in range(i0, i1 + 1):
                for j in range(j0, j1 + 1):
                    self.cells.setdefault((i, j), []).append(elem_id)

    def candidate_pairs(self) -> Iterator[tuple[int, int]]:
        seen: set[tuple[int, int]] = set()
        for elems in self.cells.values():
            k = len(elems)
            if k < 2:
                continue
            for a in range(k):
                ea = elems[a]
                for b in range(a + 1, k):
                    eb = elems[b]
                    pair = (ea, eb) if ea < eb else (eb, ea)
                    if pair in seen:
                        continue
                    seen.add(pair)
                    if _bbox_overlap(self.bboxes[pair[0]], self.bboxes[pair[1]]):
                        yield pair


def _bbox_overlap(a: np.ndarray, b: np.ndarray) -> bool:
    return not (a[2] < b[0] or b[2] < a[0] or a[3] < b[1] or b[3] < a[1])


def all_pairs_naive(n_elems: int, bboxes: np.ndarray) -> Iterator[tuple[int, int]]:
    for i in range(n_elems):
        for j in range(i + 1, n_elems):
            if _bbox_overlap(bboxes[i], bboxes[j]):
                yield (i, j)
