"""Mesh quality and coverage checks for planar triangle and quad meshes."""
from collections import Counter, defaultdict
from dataclasses import dataclass, fields
import math

import numpy as np
from matplotlib.path import Path
from scipy.spatial import cKDTree


def _cross(a, b):
    return float(a[0] * b[1] - a[1] * b[0])


def _area(p):
    return sum(_cross(p[i], p[(i + 1) % len(p)]) for i in range(len(p))) / 2


def _input(mesh_or_points, elements=None):
    if elements is None:
        points = np.asarray(mesh_or_points.points, dtype=float)[:, :2]
        rows = mesh_or_points.connectivity_list
        elements = [list(dict.fromkeys(map(int, row))) for row in rows]
    else:
        points = np.asarray(mesh_or_points, dtype=float)
        elements = [list(map(int, row)) for row in elements]
    if points.ndim != 2 or points.shape[1] != 2 or not np.isfinite(points).all():
        raise ValueError('points must be a finite (N, 2) array')
    if any(len(e) not in (3, 4) or len(set(e)) != len(e) or min(e) < 0 or max(e) >= len(points) for e in elements):
        raise ValueError('elements must contain three or four distinct valid vertex indices')
    return points, elements


def _contains(poly, p):
    return Path(np.vstack((poly, poly[0]))).contains_point(p)


def _faces(polygons, tol):
    """Split a planar segment arrangement into bounded faces and sample each face."""
    segments = [(a, b) for poly in polygons for a, b in zip(poly, np.roll(poly, -1, axis=0)) if np.linalg.norm(b-a) > tol]
    cuts = [[0., 1.] for _ in segments]
    for i, (a, b) in enumerate(segments):
        u = b-a
        for j in range(i):
            c, d = segments[j]; v = d-c
            den = _cross(u, v)
            if abs(den) > 1e-14 * max(np.linalg.norm(u)*np.linalg.norm(v), 1e-300):
                t, s = _cross(c-a, v)/den, _cross(c-a, u)/den
                if -1e-12 <= t <= 1+1e-12 and -1e-12 <= s <= 1+1e-12:
                    cuts[i].append(np.clip(t, 0, 1)); cuts[j].append(np.clip(s, 0, 1))
            elif abs(_cross(c-a, u)) <= tol*np.linalg.norm(u):
                for p in (c, d):
                    t = np.dot(p-a, u)/np.dot(u, u)
                    if 0 <= t <= 1: cuts[i].append(t)
                for p in (a, b):
                    t = np.dot(p-c, v)/np.dot(v, v)
                    if 0 <= t <= 1: cuts[j].append(t)
    vertices = []; adjacency = defaultdict(set)
    def vertex(p):
        for i, q in enumerate(vertices):
            if np.linalg.norm(p-q) <= max(tol, 1e-12): return i
        vertices.append(p); return len(vertices)-1
    for (a, b), ts in zip(segments, cuts):
        ts = sorted(set(ts))
        for x, y in zip(ts, ts[1:]):
            if (y-x)*np.linalg.norm(b-a) <= tol: continue
            i, j = vertex(a+x*(b-a)), vertex(a+y*(b-a))
            if i != j: adjacency[i].add(j); adjacency[j].add(i)
    ordered = {i: sorted(js, key=lambda j: math.atan2(*(vertices[j]-vertices[i])[::-1])) for i, js in adjacency.items()}
    visited = set(); result = []
    for i in ordered:
        for j in ordered[i]:
            if (i,j) in visited: continue
            cycle = []; a,b = i,j
            while (a,b) not in visited:
                visited.add((a,b)); cycle.append(a)
                around = ordered[b]; a,b = b, around[(around.index(a)-1) % len(around)]
            poly = np.array([vertices[k] for k in cycle])
            if len(poly) < 3 or _area(poly) <= tol*tol: continue
            edge = poly[1]-poly[0]; mid = (poly[0]+poly[1])/2
            offset = min(np.linalg.norm(edge)*1e-6, max(tol*2, np.linalg.norm(edge)*1e-9))
            sample = mid + offset*np.array([-edge[1], edge[0]])/np.linalg.norm(edge)
            result.append((poly, sample))
    return result


@dataclass
class ScoreSheet:
    schema_version: str
    points: np.ndarray
    elements: list
    element_measures: list
    counts: dict
    irregular_vertices: dict
    statistics: dict
    holes: int
    meshed_share: object
    defects: dict
    clean: bool
    checks: list
    hole_polygons: list

    def add_check(self, name, passed, detail=''):
        if any(c['name'] == name for c in self.checks):
            raise ValueError('check name already exists')
        self.checks.append(dict(name=str(name), passed=bool(passed), detail=str(detail)))

    def to_dict(self):
        def convert(x):
            if isinstance(x, np.ndarray): return convert(x.tolist())
            if isinstance(x, np.generic): return convert(x.item())
            if isinstance(x, float) and not math.isfinite(x): return None
            if isinstance(x, dict): return {k: convert(v) for k,v in x.items()}
            if isinstance(x, (list, tuple)): return [convert(v) for v in x]
            return x
        return {f.name: convert(getattr(self, f.name)) for f in fields(self)}


def score_mesh(mesh_or_points, elements=None, *, region=None, target_size=None, tol=None):
    points, elements = _input(mesh_or_points, elements)
    polygons = [points[e] for e in elements]
    edges = [np.linalg.norm(p[(i+1)%len(p)]-p[i]) for p in polygons for i in range(len(p))]
    if tol is None: tol = 1e-9 * min((x for x in edges if x > 0), default=1.)
    if tol < 0 or not np.isfinite(tol): raise ValueError('tol must be finite and nonnegative')
    if region is not None:
        region = np.asarray(region, dtype=float)
        if region.ndim != 2 or region.shape[1] != 2 or len(region) < 3 or _area(region) <= 0: raise ValueError('region must be a CCW polygon')
    measures = []; inverted = []; nonconvex = []; edge_count = Counter(); incident = defaultdict(list); neighbors = defaultdict(set)
    for eid, (e, p) in enumerate(zip(elements, polygons)):
        n = len(e); lengths = np.linalg.norm(np.roll(p,-1,axis=0)-p, axis=1); vectors = np.roll(p,-1,axis=0)-p
        signed = _area(p)
        if signed <= 0: inverted.append(eid)
        angles = []; jac = []
        for k in range(n):
            a, b = vectors[k], -vectors[k-1]; denom = lengths[k]*lengths[k-1]
            angles.append(float(np.degrees(np.arctan2(_cross(a,b), np.dot(a,b))) % 360) if denom else float('nan'))
            jac.append(_cross(a,b)/denom if denom else float('nan'))
            u,v = e[k],e[(k+1)%n]; edge_count[tuple(sorted((u,v)))] += 1; neighbors[u].add(v); neighbors[v].add(u); incident[u].append(eid)
        # Stimpson et al., The Verdict Geometric Quality Library, SAND2007-1751 (2007), triangle and quadrilateral metrics.
        scaled = float(np.clip(min(jac)*(2/math.sqrt(3) if n == 3 else 1), -1, 1))
        if n == 4 and signed > 0 and scaled <= 0:
            nonconvex.append(eid)
        stretch = taper = quality = float('nan')
        if n == 4:
            diagonal = max(np.linalg.norm(p[2]-p[0]), np.linalg.norm(p[3]-p[1]))
            stretch = math.sqrt(2)*min(lengths)/diagonal if diagonal else float('nan')
            x1 = p[1]-p[0]+p[2]-p[3]; x2 = p[2]-p[1]+p[3]-p[0]; x12 = p[0]-p[1]+p[2]-p[3]
            axis = min(np.linalg.norm(x1),np.linalg.norm(x2)); taper = np.linalg.norm(x12)/axis if axis else float('nan')
            # Folded or inverted quads have no positive edge-angle quality.
            quality = 0.0 if scaled <= 0 else (math.sqrt(stretch*min(angles)/max(angles)) if max(angles) else float('nan'))
        size = math.sqrt(abs(signed)); target = target_size(np.mean(p,axis=0)) if callable(target_size) else target_size
        if target is not None and (not np.isfinite(target) or target <= 0): raise ValueError('target size must be positive and finite')
        measures.append(dict(edge_angle_quality=quality, scaled_jacobian=scaled, stretch=stretch, taper=taper, min_angle=min(angles), max_angle=max(angles), angles=angles, size=size, size_ratio=size/target if target is not None else float('nan')))
    boundary = {v for edge,count in edge_count.items() if count == 1 for v in edge}
    irregular = [v for v in incident if v not in boundary and all(len(elements[i]) == 4 for i in incident[v]) and len(neighbors[v]) != 4]
    irregular_info = dict(count=len(irregular), degree_3=sorted(v for v in irregular if len(neighbors[v])==3), degree_5_plus=sorted(v for v in irregular if len(neighbors[v])>=5))
    duplicate = sorted([list(pair) for pair in cKDTree(points).query_pairs(tol)]) if len(points) else []
    # Query edge midpoints with a radius that encloses each edge. Exact distance
    # tests below decide membership; the tree only removes impossible candidates.
    edge_list = list(edge_count)
    edge_points = points[np.asarray(edge_list, dtype=int)] if edge_list else np.empty((0, 2, 2))
    edge_mid = edge_points.mean(axis=1)
    half_lengths = np.linalg.norm(edge_points[:, 1] - edge_points[:, 0], axis=1) / 2
    edge_tree = cKDTree(edge_mid) if len(edge_mid) else None
    hanging = []
    if edge_tree is not None:
        max_half_length = float(np.max(half_lengths))
        for v, p in enumerate(points):
            for idx in edge_tree.query_ball_point(p, max_half_length + tol):
                a, b = edge_list[idx]
                if v in (a,b) or np.linalg.norm(p-edge_mid[idx]) > half_lengths[idx]+tol: continue
                d = points[b]-points[a]; length2 = np.dot(d,d)
                if not length2: continue
                t = np.dot(p-points[a],d)/length2
                if tol/math.sqrt(length2) < t < 1-tol/math.sqrt(length2) and np.linalg.norm(p-(points[a]+t*d)) <= tol:
                    hanging.append(v); break
    overlap = []
    mins = np.array([np.min(p,axis=0) for p in polygons]).reshape(-1,2)
    maxs = np.array([np.max(p,axis=0) for p in polygons]).reshape(-1,2)
    centers = (mins+maxs)/2
    radii = np.linalg.norm((maxs-mins)/2,axis=1)
    element_tree = cKDTree(centers) if len(centers) else None
    if element_tree is not None:
        max_radius = float(np.max(radii))
        for i,p in enumerate(polygons):
            for j in sorted(element_tree.query_ball_point(centers[i], radii[i]+max_radius)):
                if j >= i: break
                if np.any(maxs[i] <= mins[j]+tol) or np.any(maxs[j] <= mins[i]+tol): continue
                q = polygons[j]
                faces = _faces([p,q], tol)
                if any(_contains(p,s) and _contains(q,s) for _,s in faces): overlap.append([j,i])
    # The arrangement gives a partition when all components are connected. Nested
    # disconnected loops are accounted for by subtracting their contained faces.
    # Conforming non-overlapping meshes need only their boundary segments for
    # coverage; internal edges do not change the covered set.
    if not overlap and not hanging and not duplicate:
        boundary_edges = [edge for edge, count in edge_count.items() if count == 1]
        adjacency = defaultdict(list)
        for a,b in boundary_edges:
            adjacency[a].append(b); adjacency[b].append(a)
        loops = []
        remaining = set(boundary_edges)
        while remaining:
            a,b = next(iter(remaining)); remaining.remove((a,b)); path = [a,b]
            while b != path[0] and len(adjacency[b]) == 2:
                c = next(v for v in adjacency[b] if v != a)
                edge = tuple(sorted((b,c)))
                if edge not in remaining: break
                remaining.remove(edge); path.append(c); a,b = b,c
            if path[-1] == path[0] and len(path) > 3: loops.append(points[path[:-1]])
        coverage_polygons = ([region] if region is not None else []) + loops
    else:
        coverage_polygons = ([region] if region is not None else []) + polygons
    all_faces = _faces(coverage_polygons, tol)
    # Disconnected nested boundaries produce one cycle per boundary. Deduct
    # immediately nested cycles to recover the area of each arrangement face.
    face_areas = [abs(_area(poly)) for poly, _ in all_faces]
    parents = [None] * len(all_faces)
    for i, (poly, sample) in enumerate(all_faces):
        containers = [j for j, (other, _) in enumerate(all_faces)
                      if face_areas[j] > face_areas[i] + tol * tol
                      and _contains(other, sample)]
        if containers:
            parents[i] = min(containers, key=lambda j: face_areas[j])
    net_areas = face_areas.copy()
    for i, parent in enumerate(parents):
        if parent is not None:
            net_areas[parent] -= face_areas[i]
    gaps = []
    uncovered = 0.
    for (poly, sample), area in zip(all_faces, net_areas):
        if area > tol * tol and (region is None or _contains(region, sample)) and not any(_contains(polygons[j], sample) for j in (element_tree.query_ball_point(sample, max_radius) if element_tree is not None else [])):
            gaps.append(poly)
            uncovered += area
    holes = len(gaps)
    share = None if region is None else max(0., min(1., 1-uncovered/_area(region)))
    quad = [m for e,m in zip(elements,measures) if len(e)==4]
    def stat(key, fn):
        values = [m[key] for m in quad if np.isfinite(m[key])]
        return float(fn(values)) if values else None
    stats = {f'{fn}_{key}': stat(key, func) for key in ('edge_angle_quality','scaled_jacobian') for fn,func in (('mean',np.mean),('min',np.min))}
    stats.update(mean_stretch=stat('stretch',np.mean), mean_taper=stat('taper',np.mean), mean_min_angle_deviation=float(np.mean([abs(m['min_angle']-90) for m in quad])) if quad else None, mean_max_angle_deviation=float(np.mean([abs(m['max_angle']-90) for m in quad])) if quad else None)
    ratios = [m['size_ratio'] for m in measures if np.isfinite(m['size_ratio'])]
    stats['size_ratio'] = dict(mean=float(np.mean(ratios)), p5=float(np.percentile(ratios,5)), p95=float(np.percentile(ratios,95)), share_within=float(np.mean([(0.7 <= x <= 1.4) for x in ratios]))) if ratios else None
    defects = {name: {'ids': ids, 'count': len(ids)} for name,ids in [('inverted',inverted),('nonconvex',nonconvex),('overlap',overlap),('hanging_node',hanging),('duplicate_vertex',duplicate),('gap',list(range(holes)))]}
    clean = not any(d['count'] for d in defects.values()); nq = len(quad); nt = len(elements)-nq
    return ScoreSheet('1.0', points, elements, measures, dict(elements=len(elements), triangles=nt, quads=nq, quad_share=nq/len(elements) if elements else 0.), irregular_info, stats, holes, share, defects, clean, [dict(name='no_defects',passed=clean,detail=f"{sum(d['count'] for d in defects.values())} defects")], gaps)
