"""Shared geometry helpers for tessellate sugar and walk."""

from __future__ import annotations

import numpy as np

POINT_TOL = 1e-5
AREA_TOL = 1e-6


def check_start(start) -> str:
    s = "" if start is None else str(start)
    if s and any(c not in "01" for c in s):
        raise ValueError("start must be a bitstring of 0/1 (empty for a fresh triangle)")
    return s


UNMARKED = 0
INSIDE = 1
OUTSIDE = 2


def prepare_triangle(t) -> np.ndarray:
    t = np.asarray(t, dtype=np.float64)
    if t.shape != (3, 2):
        raise ValueError(f"t must have shape (3, 2), got {t.shape}")
    signed = float(
        (t[1, 0] - t[0, 0]) * (t[2, 1] - t[0, 1])
        - (t[1, 1] - t[0, 1]) * (t[2, 0] - t[0, 0])
    )
    if signed < 0:
        t = t[[0, 2, 1]].copy()
    return t


def ring_vertices(p) -> np.ndarray:
    """Drop a repeated closing vertex; reverse clockwise rings to CCW."""
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


def edges_of(p) -> np.ndarray:
    """Closed ring ``p`` (N, 2) → constraint segments (N, 2, 2), CCW.

    ``out[i] = [[p[i,0], p[i,1]], [p[i+1,0], p[i+1,1]]]``  (last wraps to first).
    """
    ring = ring_vertices(p)
    return np.stack([ring, np.roll(ring, -1, axis=0)], axis=1)


def orient2d(a, b, c) -> float:
    return float((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))


def centroid(tri) -> np.ndarray:
    return np.mean(np.asarray(tri, dtype=np.float64), axis=0)


def point_on_segment(p, a, b, tol: float = POINT_TOL) -> bool:
    p, a, b = np.asarray(p, float), np.asarray(a, float), np.asarray(b, float)
    ab = b - a
    ap = p - a
    len2 = float(ab @ ab)
    if len2 < tol * tol:
        return float(ap @ ap) < tol * tol
    cross = abs(ab[0] * ap[1] - ab[1] * ap[0])
    if cross > tol * (len2 ** 0.5):
        return False
    t = float(ap @ ab) / len2
    return -1e-9 <= t <= 1.0 + 1e-9


def segment_on_triangle_side(a, b, tri, tol: float = POINT_TOL) -> bool:
    tri = np.asarray(tri, dtype=np.float64)
    for i in range(3):
        u, v = tri[i], tri[(i + 1) % 3]
        if point_on_segment(a, u, v, tol) and point_on_segment(b, u, v, tol):
            return True
    return False
