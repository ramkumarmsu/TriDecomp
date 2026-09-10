"""Read a TIGER/Line polygon ring from an ESRI shapefile or zip."""

from __future__ import annotations

import os
from typing import Optional

import numpy as np


def read_tiger_polygon(
    path: str,
    *,
    state: Optional[str] = None,
    county: Optional[str] = None,
) -> tuple[np.ndarray, dict]:
    """Return the first matching exterior ring as an (N, 2) array.

    ``state`` / ``county`` are 2- and 3-digit FIPS strings (e.g. ``"31"``, ``"039"``).
    If the record is a MultiPolygon, the part with largest area is used.
    The ring is not forced CCW here; ``tridecomp`` does that.
    """
    try:
        import pyproj
        os.environ.setdefault("PROJ_DATA", pyproj.datadir.get_data_dir())
        os.environ.setdefault("PROJ_LIB", os.environ["PROJ_DATA"])
    except Exception:
        pass
    import pyogrio

    where = None
    clauses = []
    if state is not None:
        clauses.append(f"STATEFP='{str(state).zfill(2)}'")
    if county is not None:
        clauses.append(f"COUNTYFP='{str(county).zfill(3)}'")
    if clauses:
        where = " AND ".join(clauses)

    df = pyogrio.read_dataframe(path, where=where)
    if df.empty:
        raise ValueError(f"no polygon matched {path} where={where}")

    row = df.iloc[0]
    g = row.geometry
    if g.geom_type == "MultiPolygon":
        g = max(g.geoms, key=lambda p: p.area)
    elif g.geom_type != "Polygon":
        raise ValueError(f"expected Polygon, got {g.geom_type}")

    coords = np.asarray(g.exterior.coords, dtype=np.float64)
    meta = {
        "name": str(row["NAME"]) if "NAME" in df.columns else "",
        "statefp": str(row["STATEFP"]) if "STATEFP" in df.columns else "",
        "countyfp": str(row["COUNTYFP"]) if "COUNTYFP" in df.columns else "",
        "geoid": str(row["GEOID"]) if "GEOID" in df.columns else "",
        "n_src": int(coords.shape[0]),
    }
    return coords, meta
