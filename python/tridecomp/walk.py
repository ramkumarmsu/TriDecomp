"""PerWalk marking over a tessellation (Python)."""

from __future__ import annotations

from collections import deque

import numpy as np

from .geom import (
    AREA_TOL,
    INSIDE,
    OUTSIDE,
    POINT_TOL,
    UNMARKED,
    centroid,
    check_start,
    orient2d,
    prepare_triangle,
    ring_vertices,
    segment_on_triangle_side,
)


def tiles(events, t, start: str = "") -> dict:
    """Replay splits; return {id: triangle ndarray (3,2)} for current leaves."""
    start = check_start(start)
    t = prepare_triangle(t)
    out = {start: t.copy()}
    for e in sorted(events, key=lambda z: (len(z["id"]), z["id"])):
        parent = e["id"]
        if parent not in out:
            continue
        del out[parent]
        c0, c1 = e["child"]
        out[e["child0"]] = np.asarray(c0, dtype=np.float64)
        out[e["child1"]] = np.asarray(c1, dtype=np.float64)
    return {k: v for k, v in out.items() if _tri_area(v) > AREA_TOL}


def _tri_area(tri) -> float:
    tri = np.asarray(tri, dtype=np.float64)
    return 0.5 * abs(orient2d(tri[0], tri[1], tri[2]))


def _split_points_on_edge(a, b, events):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    ab = b - a
    len2 = float(ab @ ab)
    if len2 == 0.0:
        return []
    minxy = np.minimum(a, b) - 1e-9
    maxxy = np.maximum(a, b) + 1e-9
    found = []
    for e in events:
        for ls in e["line_splits"]:
            p = np.asarray(ls["point"], dtype=np.float64)
            if np.any(p < minxy) or np.any(p > maxxy):
                continue
            ap = p - a
            if abs(ab[0] * ap[1] - ab[1] * ap[0]) > POINT_TOL * (len2 ** 0.5):
                continue
            tt = float(ap @ ab) / len2
            if 0.0 < tt < 1.0:
                found.append((tt, p))
    found.sort(key=lambda z: z[0])
    uniq = []
    for tt, p in found:
        if not uniq or float((p - uniq[-1][1]) @ (p - uniq[-1][1])) > POINT_TOL * POINT_TOL:
            uniq.append((tt, p))
    return [p for _, p in uniq]


def refined_cycle(cycle, events) -> np.ndarray:
    ring = ring_vertices(cycle)
    n = ring.shape[0]
    pts = []
    for i in range(n):
        a, b = ring[i], ring[(i + 1) % n]
        pts.append(a)
        pts.extend(_split_points_on_edge(a, b, events))
    return np.asarray(pts, dtype=np.float64)


def _left_right_tiles(a, b, leaf_tiles):
    left_id = right_id = None
    for tid, tri in leaf_tiles.items():
        if not segment_on_triangle_side(a, b, tri):
            continue
        c = centroid(tri)
        o = orient2d(a, b, c)
        if o > POINT_TOL:
            left_id = tid
        elif o < -POINT_TOL:
            right_id = tid
    return left_id, right_id


def _try_mark(marks: dict, tid, want: int) -> bool:
    if tid is None:
        return False
    cur = marks.get(tid, UNMARKED)
    if cur == UNMARKED:
        marks[tid] = want
        return True
    if cur == want:
        return True
    return False


def _collinear_overlap(a, b, c, d, tol=POINT_TOL) -> bool:
    a, b, c, d = map(lambda z: np.asarray(z, float), (a, b, c, d))
    ab = b - a
    len2 = float(ab @ ab)
    if len2 < tol * tol:
        return False
    inv = len2 ** 0.5
    if abs(ab[0] * (c - a)[1] - ab[1] * (c - a)[0]) > tol * inv:
        return False
    if abs(ab[0] * (d - a)[1] - ab[1] * (d - a)[0]) > tol * inv:
        return False
    t0 = float((c - a) @ ab) / len2
    t1 = float((d - a) @ ab) / len2
    lo, hi = (min(t0, t1), max(t0, t1))
    return hi > 1e-12 and lo < 1.0 - 1e-12 and min(hi, 1.0) - max(lo, 0.0) > 1e-12


def _on_cycle_wall(a, b, cycle_edges) -> bool:
    for c, d in cycle_edges:
        if _collinear_overlap(a, b, c, d):
            return True
    return False


def _adjacency_without_cycle(leaf_tiles: dict, cycle_pts: np.ndarray) -> dict:
    n = cycle_pts.shape[0]
    walls = [(cycle_pts[i], cycle_pts[(i + 1) % n]) for i in range(n)]
    ids = list(leaf_tiles)
    sides = []
    for tid in ids:
        tri = leaf_tiles[tid]
        for s in range(3):
            sides.append((tid, tri[s], tri[(s + 1) % 3]))
    adj = {tid: set() for tid in ids}
    for i in range(len(sides)):
        ti, a, b = sides[i]
        if _on_cycle_wall(a, b, walls):
            continue
        for j in range(i + 1, len(sides)):
            tj, c, d = sides[j]
            if ti == tj:
                continue
            if _on_cycle_wall(c, d, walls):
                continue
            if _collinear_overlap(a, b, c, d):
                adj[ti].add(tj)
                adj[tj].add(ti)
    return adj


def _flood(marks: dict, adj: dict, seed_value: int) -> None:
    q = deque(tid for tid, m in marks.items() if m == seed_value)
    seen = set(q)
    while q:
        u = q.popleft()
        for v in adj.get(u, ()):
            if v in seen:
                continue
            cur = marks.get(v, UNMARKED)
            if cur == UNMARKED:
                marks[v] = seed_value
                seen.add(v)
                q.append(v)
            elif cur == seed_value:
                seen.add(v)
                q.append(v)


def walk(events, t, cycle, start: str = "") -> dict:
    """Mark leaves by walking ``cycle``. Returns marks and a simplicity flag."""
    leaf_tiles = tiles(events, t, start=start)
    marks = {tid: UNMARKED for tid in leaf_tiles}
    cycle_pts = refined_cycle(cycle, events)
    n = cycle_pts.shape[0]
    simple = True
    missing = 0

    for i in range(n):
        a, b = cycle_pts[i], cycle_pts[(i + 1) % n]
        if float((b - a) @ (b - a)) < POINT_TOL * POINT_TOL:
            continue
        left, right = _left_right_tiles(a, b, leaf_tiles)
        if left is None or right is None:
            missing += 1
            simple = False
            continue
        if not _try_mark(marks, left, INSIDE):
            simple = False
        if not _try_mark(marks, right, OUTSIDE):
            simple = False

    adj = _adjacency_without_cycle(leaf_tiles, cycle_pts)
    _flood(marks, adj, INSIDE)
    _flood(marks, adj, OUTSIDE)

    unmarked = [tid for tid, m in marks.items() if m == UNMARKED]
    if unmarked:
        simple = False

    n_in = sum(1 for m in marks.values() if m == INSIDE)
    n_out = sum(1 for m in marks.values() if m == OUTSIDE)
    return {
        "marks": marks,
        "simple": bool(simple),
        "tiles": leaf_tiles,
        "n_inside": n_in,
        "n_outside": n_out,
        "n_unmarked": len(unmarked),
        "missing_sides": missing,
    }
