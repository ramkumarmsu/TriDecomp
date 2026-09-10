"""Python triangular decomposition: polygon + enclosing triangle in, lists out.

Default engine is the C++ EPEC library (same algorithm as ./TriDecomp).
Pass engine='python' for the CGAL-SWIG / EPIC port.
"""

from .core import tridecomp as tridecomp_python
from .gettri import gettri
from .read_tiger import read_tiger_polygon

__all__ = ["tridecomp", "tridecomp_python", "gettri", "read_tiger_polygon"]


def tridecomp(p, t, *, engine: str = "cpp", verbose: bool = False):
    """Decompose polygon ``p`` inside enclosing triangle ``t``.

    engine:
      'cpp'    — CGAL Exact_predicates_exact_constructions (requires ``make python``)
      'python' — CGAL Python / EPIC port
    """
    key = engine.lower()
    if key in ("python", "py", "epic"):
        return tridecomp_python(p, t, verbose=verbose)
    if key in ("cpp", "c++", "epec", "cgal"):
        try:
            from . import _cpp
        except ImportError as e:
            raise ImportError(
                "C++ engine not built. From the repo root run:  make python\n"
                f"(original error: {e})"
            ) from e
        return _cpp.tridecomp(p, t, verbose=verbose)
    raise ValueError(f"engine must be 'cpp' or 'python', got {engine!r}")
