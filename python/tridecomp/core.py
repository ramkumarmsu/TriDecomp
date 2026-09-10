"""EPIC stand-in for ``tessellate(segments, t)``.

Ports the C++ greedy angle-method path (lookahead 0). Geometry constructions
use the CGAL Python kernel; angle scoring stays in float64, as in C++.
"""

from __future__ import annotations

import bisect
from collections import deque

import numpy as np
from CGAL.CGAL_Kernel import (
    Line_2,
    ON_BOUNDARY,
    ON_BOUNDED_SIDE,
    ON_NEGATIVE_SIDE,
    ON_ORIENTED_BOUNDARY,
    ON_POSITIVE_SIDE,
    Point_2,
    Segment_2,
    Triangle_2,
    collinear,
    intersection,
    squared_distance,
)

POINT_TOL = 1e-5
AREA_TOL = 1e-6
MAX_PER_SIDE = 50
EMBED_BONUS = 1000.0
AEPS = 1e-9

REGION_BOUNDARY = 0
REGION_INSIDE = 1
REGION_OUTSIDE = 2

Q_LON_SCALE = 4294967296.0 / 360.0
Q_LAT_SCALE = 2147483647.0 / 90.0


# ---------------------------------------------------------------------------
# CGAL / numpy helpers
# ---------------------------------------------------------------------------

def _pt(xy) -> Point_2:
    return Point_2(float(xy[0]), float(xy[1]))


def _xy(p: Point_2) -> np.ndarray:
    return np.array([float(p.x()), float(p.y())], dtype=np.float64)


def _tri_xy(tri) -> np.ndarray:
    return np.stack([_xy(tri[0]), _xy(tri[1]), _xy(tri[2])])


def _edge_xy(e: Segment_2) -> np.ndarray:
    return np.stack([_xy(e.source()), _xy(e.target())])


def _point_eq(a: Point_2, b: Point_2, tol: float = POINT_TOL) -> bool:
    return float(squared_distance(a, b)) < tol * tol


def _mid(a: Point_2, b: Point_2) -> Point_2:
    return Point_2(0.5 * (float(a.x()) + float(b.x())),
                   0.5 * (float(a.y()) + float(b.y())))


def _area(tri) -> float:
    return abs(float(Triangle_2(tri[0], tri[1], tri[2]).area()))


def _shoelace(pts: np.ndarray) -> float:
    x, y = pts[:, 0], pts[:, 1]
    return 0.5 * float(np.abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def _prepare_polygon(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=np.float64)
    if p.ndim != 2 or p.shape[1] != 2:
        raise ValueError(f"p must have shape (N, 2), got {p.shape}")
    if p.shape[0] < 3:
        raise ValueError("p must have at least 3 vertices")
    if np.allclose(p[0], p[-1], atol=POINT_TOL, rtol=0.0):
        p = p[:-1]
        if p.shape[0] < 3:
            raise ValueError("p has fewer than 3 vertices after dropping the closing point")
    signed = float(np.dot(p[:, 0], np.roll(p[:, 1], -1)) - np.dot(p[:, 1], np.roll(p[:, 0], -1)))
    if signed < 0:
        p = p[::-1].copy()
    return p


def _prepare_triangle(t: np.ndarray) -> np.ndarray:
    t = np.asarray(t, dtype=np.float64)
    if t.shape != (3, 2):
        raise ValueError(f"t must have shape (3, 2), got {t.shape}")
    if _shoelace(t) < 0:
        t = t[[0, 2, 1]].copy()
    return t


def _quantize_lonlat(lon: float, lat: float) -> tuple[int, int]:
    xn = (lon + 180.0) * Q_LON_SCALE
    xn = min(max(xn, 0.0), 4294967295.0)
    yn = lat * Q_LAT_SCALE
    yn = min(max(yn, -2147483647.0), 2147483647.0)
    return int(round(xn)), int(round(yn))


def _quantize_pt(p: Point_2) -> tuple[int, int]:
    return _quantize_lonlat(float(p.x()), float(p.y()))


# ---------------------------------------------------------------------------
# Predicates used by the split loop
# ---------------------------------------------------------------------------

def _edge_on_triangle_side(edge: Segment_2, tri) -> bool:
    for i in range(3):
        side = Segment_2(tri[i], tri[(i + 1) % 3])
        if side.has_on(edge.source()) and side.has_on(edge.target()):
            return True
    return False


def _filter_degenerate(edges):
    return [e for e in edges if not _point_eq(e.source(), e.target())]


def _split_line_is_polygon_edge(split_line: Segment_2, polygon_edges) -> bool:
    sl = split_line.supporting_line()
    for pe in polygon_edges:
        if sl == pe.supporting_line():
            if pe.has_on(split_line.source()) and pe.has_on(split_line.target()):
                return True
            if split_line.has_on(pe.source()) and split_line.has_on(pe.target()):
                return True
    return False


def _edge_may_meet_cut(cut_line: Line_2, edge: Segment_2) -> bool:
    sa = cut_line.oriented_side(edge.source())
    sb = cut_line.oriented_side(edge.target())
    if sa == ON_ORIENTED_BOUNDARY or sb == ON_ORIENTED_BOUNDARY:
        return True
    return (sa == ON_POSITIVE_SIDE and sb == ON_NEGATIVE_SIDE) or (
        sa == ON_NEGATIVE_SIDE and sb == ON_POSITIVE_SIDE
    )


def _place_edge_in_child(seg: Segment_2, ct, child_edges) -> None:
    if _point_eq(seg.source(), seg.target()):
        return
    mid = _mid(seg.source(), seg.target())
    for i in range(2):
        side = ct[i].bounded_side(mid)
        if side == ON_BOUNDED_SIDE or side == ON_BOUNDARY:
            child_edges[i].append(seg)
            return


def _edge_lies_on_cut(edge: Segment_2, split_line: Segment_2) -> bool:
    if _point_eq(edge.source(), edge.target()):
        return False
    return split_line.has_on(edge.source()) and split_line.has_on(edge.target())


def _assign_edges(child_tris, edges, split_line: Segment_2, track: bool,
                  triangle_split_id: int, iteration: int, line_splits: list):
    child_edges = [[], []]
    ct = (
        Triangle_2(child_tris[0][0], child_tris[0][1], child_tris[0][2]),
        Triangle_2(child_tris[1][0], child_tris[1][1], child_tris[1][2]),
    )
    cut_line = split_line.supporting_line()
    for edge in edges:
        if _edge_lies_on_cut(edge, split_line):
            child_edges[0].append(edge)
            child_edges[1].append(edge)
            continue
        if not _edge_may_meet_cut(cut_line, edge):
            _place_edge_in_child(edge, ct, child_edges)
            continue
        obj = intersection(edge, split_line)
        if obj.empty():
            _place_edge_in_child(edge, ct, child_edges)
            continue
        if obj.is_Point_2():
            p = obj.get_Point_2()
            seg1 = Segment_2(edge.source(), p)
            seg2 = Segment_2(p, edge.target())
            s1_ok = not _point_eq(seg1.source(), seg1.target())
            s2_ok = not _point_eq(seg2.source(), seg2.target())
            if track and s1_ok and s2_ok:
                line_splits.append({
                    "split_id": len(line_splits) + 1,
                    "triangle_split_id": triangle_split_id,
                    "iteration": iteration,
                    "original_edge": edge,
                    "split_point": p,
                    "segment1": seg1,
                    "segment2": seg2,
                })
            if s1_ok:
                _place_edge_in_child(seg1, ct, child_edges)
            if s2_ok:
                _place_edge_in_child(seg2, ct, child_edges)
        elif obj.is_Segment_2():
            _place_edge_in_child(edge, ct, child_edges)
    return child_edges[0], child_edges[1]


# ---------------------------------------------------------------------------
# Angle method (greedy)
# ---------------------------------------------------------------------------

def _rel_angle(rx, ry, px, py) -> float:
    return abs(np.arctan2(rx * py - ry * px, rx * px + ry * py))


def _am_scored_cuts(triangle, edges):
    scored = []
    if not edges:
        return scored
    for c in range(3):
        v, u, w = triangle[c], triangle[(c + 1) % 3], triangle[(c + 2) % 3]
        vx, vy = float(v.x()), float(v.y())
        rx, ry = float(u.x()) - vx, float(u.y()) - vy
        aw = _rel_angle(rx, ry, float(w.x()) - vx, float(w.y()) - vy)
        if aw < 1e-12:
            continue
        seg_lo, seg_hi, cand = [], [], []
        for e in edges:
            sx, sy = float(e.source().x()) - vx, float(e.source().y()) - vy
            tx, ty = float(e.target().x()) - vx, float(e.target().y()) - vy
            a1 = _rel_angle(rx, ry, sx, sy)
            a2 = _rel_angle(rx, ry, tx, ty)
            src_at_v = sx * sx + sy * sy < 1e-24
            tgt_at_v = tx * tx + ty * ty < 1e-24
            if src_at_v:
                a1 = a2
            if tgt_at_v:
                a2 = a1
            seg_lo.append(min(a1, a2))
            seg_hi.append(max(a1, a2))
            # Never aim a ray at the corner itself (same angle as the far
            # endpoint after the collapse above); that makes Line_2 degenerate.
            if (not src_at_v) and 1e-9 < a1 < aw - 1e-9:
                cand.append((a1, e.source()))
            if (not tgt_at_v) and 1e-9 < a2 < aw - 1e-9:
                cand.append((a2, e.target()))
        if not cand:
            continue
        cand.sort(key=lambda z: z[0])
        uniq = []
        for ang, pt in cand:
            if not uniq or abs(ang - uniq[-1][0]) > 1e-9:
                uniq.append((ang, pt))
        if len(uniq) > 2 * MAX_PER_SIDE:
            mid = 0.5 * aw
            median_idx = min(range(len(uniq)), key=lambda i: abs(uniq[i][0] - mid))
            lo = max(0, median_idx - MAX_PER_SIDE)
            hi = min(len(uniq) - 1, median_idx + MAX_PER_SIDE)
            uniq = uniq[lo:hi + 1]
        sorted_hi = sorted(seg_hi)
        sorted_lo = sorted(seg_lo)
        n = len(edges)
        emb_angles = sorted(0.5 * (lo + hi) for lo, hi in zip(seg_lo, seg_hi) if hi - lo < AEPS)
        for ang, pt in uniq:
            raw_left = bisect.bisect_right(sorted_hi, ang)
            raw_right = n - bisect.bisect_left(sorted_lo, ang)
            embedded = (bisect.bisect_right(emb_angles, ang + AEPS)
                        - bisect.bisect_left(emb_angles, ang - AEPS))
            left = raw_left - embedded
            right = raw_right - embedded
            cuts = max(0, n - left - right - embedded)
            balance = 1.5 * abs(left - right) / (left + right) if (left + right) > 0 else 1.5
            cost = cuts + balance - EMBED_BONUS * embedded
            scored.append((cost, c, pt, cuts, embedded))
    scored.sort(key=lambda z: z[0])
    return scored


def _am_realize_cut(triangle, sc, all_polygon_edges):
    cost, corner, target, cuts, embedded = sc
    v = triangle[corner]
    u = triangle[(corner + 1) % 3]
    w = triangle[(corner + 2) % 3]
    if _point_eq(v, target):
        return None
    slope = Line_2(v, target)
    opposite = Segment_2(u, w)
    obj = intersection(opposite, slope)
    if obj.empty() or not obj.is_Point_2():
        return None
    newp = obj.get_Point_2()
    t1 = (v, u, newp)
    t2 = (v, newp, w)
    if _area(t1) <= AREA_TOL or _area(t2) <= AREA_TOL:
        return None
    split_line = Segment_2(v, newp)
    if embedded > 0 or _split_line_is_polygon_edge(split_line, all_polygon_edges):
        priority = 0
    else:
        priority = 1 if cuts == 0 else 2
    return (t1, t2), split_line, priority


def _find_best_split_angle(triangle, edges, all_polygon_edges):
    scored = _am_scored_cuts(triangle, edges)
    for sc in scored:
        realized = _am_realize_cut(triangle, sc, all_polygon_edges)
        if realized is not None:
            return realized
    return None


# ---------------------------------------------------------------------------
# PerWalk + MarkTriangles
# ---------------------------------------------------------------------------

def _pip_inside(ring_xy: np.ndarray, px: float, py: float) -> bool:
    x, y = ring_xy[:, 0], ring_xy[:, 1]
    n = len(x)
    inside = False
    j = n - 1
    for i in range(n):
        yi, yj = y[i], y[j]
        if (yi > py) != (yj > py) and (px < (x[j] - x[i]) * (py - yi) / (yj - yi) + x[i]):
            inside = not inside
        j = i
    return inside


def _mark_triangles(leaves, polygon_edges, ring_xy: np.ndarray):
    L = len(leaves)
    codes = [REGION_OUTSIDE] * L
    if L == 0:
        return {
            "codes": codes, "inside_count": 0, "outside_count": 0,
            "perwalk_seeds": 0, "flood_added": 0, "flood_unreached": 0,
            "area_ok": True, "inside_area": 0.0, "polygon_area": _shoelace(ring_xy),
            "leaf_side_keys": [], "poly_edge_by_key": {},
        }

    is_in = []
    for tri in leaves:
        c = (_xy(tri[0]) + _xy(tri[1]) + _xy(tri[2])) / 3.0
        is_in.append(_pip_inside(ring_xy, float(c[0]), float(c[1])))

    # Canonical vertex ids (POINT_TOL grid), then endpoint-key for full-side match.
    cell = POINT_TOL
    grid = {}
    dpt = []

    def canon(x, y):
        ix, iy = int(np.floor(x / cell)), int(np.floor(y / cell))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for vid in grid.get((ix + dx, iy + dy), ()):
                    ex, ey = dpt[vid][0] - x, dpt[vid][1] - y
                    if ex * ex + ey * ey < cell * cell:
                        return vid
        vid = len(dpt)
        dpt.append((x, y))
        grid.setdefault((ix, iy), []).append(vid)
        return vid

    def ekey(a, b):
        if a > b:
            a, b = b, a
        return (a << 32) | b

    leaf_side_keys = []
    for tri in leaves:
        v0 = canon(float(tri[0].x()), float(tri[0].y()))
        v1 = canon(float(tri[1].x()), float(tri[1].y()))
        v2 = canon(float(tri[2].x()), float(tri[2].y()))
        leaf_side_keys.append((ekey(v0, v1), ekey(v1, v2), ekey(v2, v0)))

    poly_edge_by_key = {}
    for pe_i, pe in enumerate(polygon_edges):
        a = canon(float(pe.source().x()), float(pe.source().y()))
        b = canon(float(pe.target().x()), float(pe.target().y()))
        poly_edge_by_key[ekey(a, b)] = pe_i

    # Segment-overlap adjacency (T-junctions).
    sref = []
    for t_i, tri in enumerate(leaves):
        for s in range(3):
            A, B = tri[s], tri[(s + 1) % 3]
            ax, ay = float(A.x()), float(A.y())
            dx, dy = float(B.x()) - ax, float(B.y()) - ay
            length = (dx * dx + dy * dy) ** 0.5
            if length < 1e-15:
                continue
            nx, ny = -dy / length, dx / length
            rho = ax * nx + ay * ny
            if nx < 0 or (nx == 0.0 and ny < 0):
                nx, ny, rho = -nx, -ny, -rho
            sref.append((rho, nx, ny, t_i, A, B))
    sref.sort(key=lambda z: (z[0], z[1], z[2]))

    adj = [[] for _ in range(L)]
    seed_touch = [False] * L
    RHO_TOL, N_TOL, OVL_TOL = 1e-9, 1e-6, 1e-9
    gs = 0
    NS = len(sref)
    while gs < NS:
        ge = gs + 1
        while (ge < NS and sref[ge][0] - sref[ge - 1][0] <= RHO_TOL
               and abs(sref[ge][1] - sref[gs][1]) <= N_TOL
               and abs(sref[ge][2] - sref[gs][2]) <= N_TOL):
            ge += 1
        K = ge - gs
        if K >= 2:
            _, _, _, _, ra, rb = sref[gs]
            ax, ay = float(ra.x()), float(ra.y())
            dx, dy = float(rb.x()) - ax, float(rb.y()) - ay
            dl = (dx * dx + dy * dy) ** 0.5
            if dl >= 1e-15:
                ux, uy = dx / dl, dy / dl
                iv, ivleaf = [], []
                for k in range(gs, ge):
                    _, _, _, leaf, sa, sb = sref[k]
                    if k != gs and (not collinear(ra, rb, sa) or not collinear(ra, rb, sb)):
                        continue
                    p0 = (float(sa.x()) - ax) * ux + (float(sa.y()) - ay) * uy
                    p1 = (float(sb.x()) - ax) * ux + (float(sb.y()) - ay) * uy
                    if p0 > p1:
                        p0, p1 = p1, p0
                    iv.append((p0, p1))
                    ivleaf.append(leaf)
                ord_i = sorted(range(len(iv)), key=lambda z: iv[z][0])
                for a in range(len(iv)):
                    ia = ord_i[a]
                    a1 = iv[ia][1]
                    for b in range(a + 1, len(iv)):
                        ib = ord_i[b]
                        if iv[ib][0] >= a1 - OVL_TOL:
                            break
                        ov = min(a1, iv[ib][1]) - max(iv[ia][0], iv[ib][0])
                        if ov <= OVL_TOL:
                            continue
                        u, v = ivleaf[ia], ivleaf[ib]
                        if u == v:
                            continue
                        iu, iv2 = is_in[u], is_in[v]
                        if iu and iv2:
                            adj[u].append(v)
                            adj[v].append(u)
                        elif iu != iv2:
                            seed_touch[u if iu else v] = True
        gs = ge

    stack = []
    reached = [False] * L
    perwalk_seeds = 0
    for t_i in range(L):
        if is_in[t_i] and seed_touch[t_i]:
            reached[t_i] = True
            stack.append(t_i)
            perwalk_seeds += 1
    flood_added = 0
    while stack:
        t_i = stack.pop()
        for nb in adj[t_i]:
            if not reached[nb]:
                reached[nb] = True
                flood_added += 1
                stack.append(nb)

    inside_count = outside_count = flood_unreached = 0
    inside_area = 0.0
    for t_i in range(L):
        if is_in[t_i]:
            codes[t_i] = REGION_INSIDE
            inside_count += 1
            if not reached[t_i]:
                flood_unreached += 1
            inside_area += _area(leaves[t_i])
        else:
            codes[t_i] = REGION_OUTSIDE
            outside_count += 1
    polygon_area = _shoelace(ring_xy)
    rel = abs(inside_area - polygon_area) / polygon_area if polygon_area > 0 else abs(inside_area)
    return {
        "codes": codes,
        "inside_count": inside_count,
        "outside_count": outside_count,
        "perwalk_seeds": perwalk_seeds,
        "flood_added": flood_added,
        "flood_unreached": flood_unreached,
        "area_ok": rel < 1e-3,
        "inside_area": inside_area,
        "polygon_area": polygon_area,
        "leaf_side_keys": leaf_side_keys,
        "poly_edge_by_key": poly_edge_by_key,
    }


# ---------------------------------------------------------------------------
# Walkable quantized boundary
# ---------------------------------------------------------------------------

def _param_on_edge(A: Point_2, B: Point_2, P: Point_2):
    ax, ay = float(A.x()), float(A.y())
    bx, by = float(B.x()), float(B.y())
    px, py = float(P.x()), float(P.y())
    dx, dy = bx - ax, by - ay
    dlen2 = dx * dx + dy * dy
    if dlen2 == 0.0:
        return None
    t = ((px - ax) * dx + (py - ay) * dy) / dlen2
    projx, projy = ax + t * dx, ay + t * dy
    if (px - projx) ** 2 + (py - projy) ** 2 > POINT_TOL * POINT_TOL:
        return None
    if t <= 0.0 or t >= 1.0:
        return None
    return t


def _build_walk(polygon_edges, line_splits):
    walk = []
    seq = 0
    for E in polygon_edges:
        A, B = E.source(), E.target()
        ax, ay = float(A.x()), float(A.y())
        bx, by = float(B.x()), float(B.y())
        minx, maxx = min(ax, bx) - 1e-9, max(ax, bx) + 1e-9
        miny, maxy = min(ay, by) - 1e-9, max(ay, by) + 1e-9
        mids = []
        seen = []
        for rec in line_splits:
            P = rec["split_point"]
            px, py = float(P.x()), float(P.y())
            if px < minx or px > maxx or py < miny or py > maxy:
                continue
            t = _param_on_edge(A, B, P)
            if t is None:
                continue
            if any(_point_eq(P, q) for q in seen):
                continue
            seen.append(P)
            mids.append((t, P))
        mids.sort(key=lambda z: z[0])
        prev = A
        pts = [p for _, p in mids] + [B]
        for nxt in pts:
            qs = _quantize_pt(prev)
            qt = _quantize_pt(nxt)
            walk.append({
                "seq": seq,
                "exact": np.stack([_xy(prev), _xy(nxt)]),
                "quantized": np.array([qs, qt], dtype=np.int64),
            })
            seq += 1
            prev = nxt
    connected = all(
        walk[i]["quantized"][1, 0] == walk[i + 1]["quantized"][0, 0]
        and walk[i]["quantized"][1, 1] == walk[i + 1]["quantized"][0, 1]
        for i in range(len(walk) - 1)
    ) if len(walk) > 1 else True
    closed = False
    if walk:
        closed = (walk[-1]["quantized"][1, 0] == walk[0]["quantized"][0, 0]
                  and walk[-1]["quantized"][1, 1] == walk[0]["quantized"][0, 1])
    walkable = bool(walk) and connected and closed
    return walk, walkable


# ---------------------------------------------------------------------------
# Public API (EPIC stand-in for C++ tessellate)
# ---------------------------------------------------------------------------

def _centroid_pt(tri) -> Point_2:
    return Point_2(
        (float(tri[0].x()) + float(tri[1].x()) + float(tri[2].x())) / 3.0,
        (float(tri[0].y()) + float(tri[1].y()) + float(tri[2].y())) / 3.0,
    )


def _order_children_left_right(split_line: Segment_2, t1, t2):
    side = split_line.supporting_line().oriented_side(_centroid_pt(t1))
    if side == ON_NEGATIVE_SIDE:
        return t2, t1
    return t1, t2


def _prepare_segments(segments: np.ndarray) -> np.ndarray:
    s = np.asarray(segments, dtype=np.float64)
    if s.size == 0:
        return np.zeros((0, 2, 2), dtype=np.float64)
    if s.ndim == 2 and s.shape[1] == 4:
        s = s.reshape(-1, 2, 2)
    if s.ndim != 3 or s.shape[1:] != (2, 2):
        raise ValueError(f"segments must have shape (M, 2, 2) or (M, 4), got {s.shape}")
    keep = np.linalg.norm(s[:, 0] - s[:, 1], axis=1) > POINT_TOL
    return s[keep]


def tessellate(segments: np.ndarray, t: np.ndarray, *, start: str = "",
               verbose: bool = False) -> list:
    """Greedy angle-method tessellation using CGAL Python / EPIC."""
    import sys
    from .geom import check_start

    start = check_start(start)
    segs_xy = _prepare_segments(segments)
    tri_xy = _prepare_triangle(t)
    constraint_edges = [Segment_2(_pt(s[0]), _pt(s[1])) for s in segs_xy]
    initial = (_pt(tri_xy[0]), _pt(tri_xy[1]), _pt(tri_xy[2]))

    queue = deque([(initial, list(constraint_edges), start)])
    max_iter = max(100_000, max(1, len(constraint_edges)) * 10)
    iteration = 0
    events = []

    while queue and iteration < max_iter:
        tri, edges, addr = queue.popleft()
        iteration += 1
        if verbose and (iteration == 1 or iteration % 100 == 0):
            print(f"[progress] iteration={iteration} queue={len(queue)} "
                  f"splits={len(events)}", file=sys.stderr)

        current = _filter_degenerate(
            [e for e in edges if not _edge_on_triangle_side(e, tri)]
        )
        if not current:
            continue

        best = _find_best_split_angle(tri, current, current)
        if best is None:
            continue
        (t1, t2), split_line, _priority = best
        t1, t2 = _order_children_left_right(split_line, t1, t2)
        line_splits = []
        ce1, ce2 = _assign_edges(
            (t1, t2), current, split_line, True, len(events) + 1, iteration, line_splits
        )
        events.append({
            "id": addr,
            "child0": addr + "0",
            "child1": addr + "1",
            "parent": _tri_xy(tri),
            "cut": _edge_xy(split_line),
            "child": (_tri_xy(t1), _tri_xy(t2)),
            "line_splits": [
                {
                    "edge": _edge_xy(rec["original_edge"]),
                    "point": _xy(rec["split_point"]),
                    "seg": (_edge_xy(rec["segment1"]), _edge_xy(rec["segment2"])),
                }
                for rec in line_splits
            ],
        })
        if _area(t1) > AREA_TOL:
            queue.append((t1, _filter_degenerate(ce1), addr + "0"))
        if _area(t2) > AREA_TOL:
            queue.append((t2, _filter_degenerate(ce2), addr + "1"))

    return events
