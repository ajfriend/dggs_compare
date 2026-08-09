"""Per-cell statistics: aspect ratio (csar) + spherical area (sparea).

Computed by the metrics stage and stored as columns in the tables —
downstream consumers read them; they never re-solve. Solver settings come
from `config` and are recorded in each table's metadata.
"""

import numpy as np
import csar
import sparea

from . import config


def refine_geodesic(ring, k):
    """Insert `k` slerp points along each great-circle edge of an open
    (lat, lng)-degree vertex list.

    For systems whose edges ARE geodesics between their vertices (h3, s2),
    this is what higher sampling density means — the enclosing cone can't
    change (convexity), so the convergence check confirms the
    implementation matches the theorem.
    """
    la, lo = np.radians(np.asarray(ring, dtype=float)).T
    v = np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo),
                         np.sin(la)])
    out = []
    n = len(v)
    for i in range(n):
        a, b = v[i], v[(i + 1) % n]
        out.append(a)
        ang = np.arccos(np.clip(a @ b, -1.0, 1.0))
        if ang > 0:
            s = np.sin(ang)
            for j in range(1, k + 1):
                t = j / (k + 1)
                out.append((np.sin((1 - t) * ang) * a + np.sin(t * ang) * b) / s)
    out = np.asarray(out)
    lat = np.degrees(np.arcsin(np.clip(out[:, 2], -1.0, 1.0)))
    lng = np.degrees(np.arctan2(out[:, 1], out[:, 0]))
    return list(zip(lat.tolist(), lng.tolist()))


def is_ccw(ring):
    """True if `ring` ([(lat, lng), ...]) is CCW seen from outside.

    THE winding predicate — the sign convention lives here and nowhere
    else (binding-neutral home: both engines need it, and importing one
    engine from the other would load its native library).
    """
    return sparea.area(ring, signed=True) >= 0


def orient_ccw(ring):
    """`ring` ([(lat, lng), ...]) in CCW-seen-from-outside order.

    sparea's unsigned area turns a CW ring into its ~4pi complement, so
    every ring an engine hands out is normalized. Decided per ring: DGGAL
    winds two pentagons per level opposite to the rest of their grid, so
    a per-grid flip would corrupt exactly those cells.
    """
    return ring if is_ccw(ring) else ring[::-1]


def ar(latlng):
    """Enclosing-cone aspect ratio of one open (lat, lng)-degree vertex
    list at the CONFIG solver settings, or None if the solve did not
    certify. The one AR entry point, so every published number (tables,
    convergence residuals) shares the same solver provenance."""
    r = csar.solve(csar.to_vec3(latlng, geo='latlng_deg'), geo='vec3',
                   gap_tol=config.GAP_TOL, method=config.CSAR_METHOD)
    return r.aspect_ratio if isinstance(r, csar.Converged) else None


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


def sample_uniform_latlng(n, rng):
    """Uniform-on-sphere samples as (lat_deg, lng_deg), shape (n, 2) —
    the registry contract's point order. (lng is drawn first to preserve
    the historical RNG stream: same seed, same points, same tables.)"""
    lng = 360.0 * rng.random(n) - 180.0
    lat = np.degrees(np.arcsin(2.0 * rng.random(n) - 1.0))  # equal-area in lat
    return np.column_stack([lat, lng])
