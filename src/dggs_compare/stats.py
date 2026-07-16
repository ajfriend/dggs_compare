"""Per-cell statistics: aspect ratio (csar) + spherical area (sparea).

Computed once, at generation time, and stored as columns in the tables —
downstream consumers read them; they never re-solve. Solver settings come
from `config` and are recorded in each table's metadata (see `cache.py`).
"""

import numpy as np
import csar
import sparea

from . import config


def cell_stats(latlng):
    """(ar, area_sr) for one open (lat, lng)-degree ring.

    ar: enclosing-cone aspect ratio, NaN if the solve did not certify at
    `config.GAP_TOL` (the f64 gap floor at the finest sub-metre resolutions).
    area_sr: spherical polygon area in steradians (multiply by
    `config.SR2KM2` for km^2).
    """
    v = csar.to_vec3(latlng, geo='latlng_deg')
    r = csar.solve(v, geo='vec3', gap_tol=config.GAP_TOL,
                   method=config.CSAR_METHOD)
    ar = r.aspect_ratio if isinstance(r, csar.Converged) else float('nan')
    return ar, float(sparea.area(v, geo='vec3'))


def sample_uniform_lnglat(n, rng):
    """Uniform-on-sphere samples as (lng_deg, lat_deg), shape (n, 2)."""
    lng = 360.0 * rng.random(n) - 180.0
    lat = np.degrees(np.arcsin(2.0 * rng.random(n) - 1.0))  # equal-area in lat
    return np.column_stack([lng, lat])
