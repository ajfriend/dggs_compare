"""Per-cell statistics (csar AR, sparea area) + shared sphere helpers.

The statistics are computed by the metrics stage and stored as columns in
the tables — downstream consumers read them; they never re-solve. Solver
settings come from `config` and are recorded in each table's metadata.
The helpers (winding, edge refinement, uniform sampling, the authalic
map) serve the generation scripts.
"""

import numpy as np
import csar
import sparea

from . import config


# WGS84 shape constants (authalic_lat only — no datum concept anywhere
# else in the library; whether and how to use the transform is each
# implementation script's decision and record).
_F = 1 / 298.257223563
_E2 = _F * (2 - _F)
_E = np.sqrt(_E2)


def _q(sinlat):
    return (1 - _E2) * (sinlat / (1 - _E2 * sinlat * sinlat)
                        - np.log((1 - _E * sinlat) / (1 + _E * sinlat))
                        / (2 * _E))


_Q_POLE = _q(1.0)


def authalic_lat(lat_deg):
    """Authalic latitude (degrees) for WGS84 geodetic latitude (degrees);
    scalar or array. Longitude is unchanged by the transform.

    The exact closed form (Snyder 1987, eq. 3-11/3-12), not a series: the
    latitude on the equal-AREA sphere, so a system that is exactly
    equal-area on the WGS84 ellipsoid is exactly equal-area on the sphere
    after mapping each vertex through this. Max |shift| is ~0.13 deg at
    mid-latitudes. Evaluated on |lat| and mirrored (the function is odd),
    which keeps both poles exact."""
    lat = np.asarray(lat_deg, dtype=float)
    s = np.sin(np.radians(np.abs(lat)))
    xi = np.degrees(np.arcsin(np.clip(_q(s) / _Q_POLE, -1.0, 1.0)))
    return np.copysign(xi, lat)


def authalic_rings(rings):
    """[[(lat, lng), ...], ...] with every vertex latitude mapped through
    `authalic_lat` — the one-call form for a `boundaries()` return value.
    An implementation whose native coordinates are WGS84 geodetic applies
    this as the last step of `boundaries`; doing so IS the script's record
    of how its system got to the sphere."""
    out = []
    for ring in rings:
        r = np.asarray(ring, dtype=float)
        out.append(list(zip(authalic_lat(r[:, 0]).tolist(),
                            r[:, 1].tolist())))
    return out


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
