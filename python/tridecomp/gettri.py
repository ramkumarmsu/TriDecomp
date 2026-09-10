"""Enclosing triangle from a 2-D point set."""

from __future__ import annotations

import numpy as np


def gettri(p: np.ndarray, rat: float = 1.5) -> np.ndarray:
    """Return a (3, 2) enclosing triangle for polygon vertices ``p``.

    ``rat`` pulls the apex west of ``min x`` by ``(max y - min y) / rat``.
    Default is 1.5, as used in the county tests.
    """
    p = np.asarray(p, dtype=np.float64)
    slope, intercept = np.polyfit(p[:, 0], p[:, 1], deg=1)
    minx, miny, maxx, maxy = np.min(p[:, 0]), np.min(p[:, 1]), np.max(p[:, 0]), np.max(p[:, 1])
    ydiff = maxy - miny
    px = minx - ydiff / rat
    py = intercept + slope * px
    P = np.array([px, py])
    vec = p - P
    angles = np.unwrap(np.arctan2(vec[:, 1], vec[:, 0]))
    mina, maxa = np.min(angles) - 0.05, np.max(angles) + 0.05
    mida = (mina + maxa) / 2.0
    dirl = np.array([np.cos(mina), np.sin(mina)])
    dirh = np.array([np.cos(maxa), np.sin(maxa)])
    dirm = np.array([np.cos(mida), np.sin(mida)])
    distances = vec @ dirm
    mxdist = np.max(distances) * 1.05
    tl = mxdist / np.dot(dirl, dirm)
    th = mxdist / np.dot(dirh, dirm)
    Q, R = P + tl * dirl, P + th * dirh
    return np.array([P, Q, R])
