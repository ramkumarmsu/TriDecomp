"""Triangular decomposition.

C++ is used only for ``tessellate(segments, t)`` (exact EPEC kernel).
Walk, marks, and polygon sugar live in Python.
"""

from .core import tessellate as tessellate_python
from .geom import INSIDE, OUTSIDE, UNMARKED, edges_of
from .gettri import gettri
from .read_tiger import read_tiger_polygon
from .walk import tiles, walk

__all__ = [
    "tessellate",
    "tessellate_python",
    "tridecomp",
    "edges_of",
    "walk",
    "tiles",
    "gettri",
    "read_tiger_polygon",
    "UNMARKED",
    "INSIDE",
    "OUTSIDE",
]


def tessellate(segments, t, *, start: str = "", engine: str = "cpp",
               verbose: bool = False):
    """Split triangle ``t`` against constraint ``segments``.

    segments : array-like, shape (M, 2, 2) or (M, 4)
        Untyped constraint edges. Each 2×2 block is two endpoints::

            segments[i] = [[x0, y0],   # start
                           [x1, y1]]   # end

        Row-major this is the same as ``(M, 4)``: ``x0, y0, x1, y1``.
        Tessellation treats the edge as undirected; walk uses cycle order
        separately. A polygon ring is one source (see ``edges_of``).
    t : array-like, shape (3, 2)
        Triangle to split: the enclosing triangle, or a leaf being refined.
    start : str
        Tree address of ``t``. Empty for a fresh polygon / enclosing triangle.
        Pass the leaf's code (e.g. ``"01011"``) to tessellate that tile further
        after adding more segments inside it; child ids are ``start+'0'`` / ``start+'1'``.
    engine : {'cpp', 'python'}
        'cpp' is CGAL EPEC (requires ``make python``). 'python' is the EPIC port.
    """
    key = engine.lower()
    if key in ("python", "py", "epic"):
        return tessellate_python(segments, t, start=start, verbose=verbose)
    if key in ("cpp", "c++", "epec", "cgal"):
        try:
            from . import _cpp
        except ImportError as e:
            raise ImportError(
                "C++ engine not built. From the repo root run:  make python\n"
                f"(original error: {e})"
            ) from e
        return _cpp.tessellate(segments, t, start=start, verbose=verbose)
    raise ValueError(f"engine must be 'cpp' or 'python', got {engine!r}")


def tridecomp(p, t, *, start: str = "", engine: str = "cpp", verbose: bool = False):
    """Sugar: ``tessellate(edges_of(p), t, start=start)``."""
    return tessellate(edges_of(p), t, start=start, engine=engine, verbose=verbose)
