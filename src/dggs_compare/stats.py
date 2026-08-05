"""Per-cell statistics: aspect ratio (csar) + spherical area (sparea).

Computed once, at generation time, and stored as columns in the tables —
downstream consumers read them; they never re-solve. Solver settings come
from `config` and are recorded in each table's metadata (see `cache.py`).
"""

import numpy as np
import csar
import sparea

from . import config


def orient_ccw(ring):
    """`ring` ([(lat, lng), ...]) in CCW-seen-from-outside order.

    Binding-neutral home (both engines need it; importing one engine from
    the other would load its native library). sparea's unsigned area turns
    a CW ring into its ~4pi complement, so every ring an engine hands out
    is normalized. Decided per ring: DGGAL winds two pentagons per level
    opposite to the rest of their grid, so a per-grid flip would corrupt
    exactly those cells.
    """
    if sparea.area(ring, signed=True) < 0:
        ring = ring[::-1]
    return ring


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
    area = float(sparea.area(v, geo='vec3'))
    if area > 2 * np.pi:
        # No cell approaches a hemisphere: this is a wound-backwards ring
        # (sparea returns the ~4pi complement) or garbage geometry. Never
        # let it reach a table silently.
        raise ValueError(f'cell area {area:.4f} sr > 2*pi — reversed ring?')
    return ar, area


def sample_uniform_lnglat(n, rng):
    """Uniform-on-sphere samples as (lng_deg, lat_deg), shape (n, 2)."""
    lng = 360.0 * rng.random(n) - 180.0
    lat = np.degrees(np.arcsin(2.0 * rng.random(n) - 1.0))  # equal-area in lat
    return np.column_stack([lng, lat])
