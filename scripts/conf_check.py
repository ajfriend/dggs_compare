# /// script
# requires-python = ">=3.11,<3.14"
# dependencies = [
#     "dggal>=0.0.6",
#     "a5_fast",
#     "hex9>=2.3.1",
#     "dggs-compare",
# ]
#
# [tool.uv.sources]
# dggs-compare = { path = "..", editable = true }
# ///
"""Authalic vs conformal spherization (#82): the analytic Tissot bound,
and a paired-sample check that real grids sit inside it.

MATH. Both maps are latitude-only ellipsoid->sphere reparameterizations:
phi -> xi(phi) (authalic) and phi -> chi(phi) (conformal). The choice
between the two spherizations of ANY geodetic-addressed grid is governed
by the sphere->sphere composite T = chi . xi^{-1}, whose Tissot scales
at geodetic latitude phi are

    h = chi'(phi) / xi'(phi)      (meridional)
    k = cos(chi) / cos(xi)        (parallel)

For an infinitesimal cell, a linear map with singular-value ratio
sigma = max(h/k, k/h) changes any convex body's aspect ratio by a factor
in [1/sigma, sigma] — so max_phi sigma is the sharp uniform bound on the
AR change, for every grid and every cell. The conformal-sphere area
field is J = h*k (area-weighted mean exactly 1, since T is a sphere
bijection): its spread is exactly what an authalic-equal-area grid pays
for choosing the conformal map.

SAMPLING. For each system, point-drawn cells at its working resolution;
each cell's RAW GEODETIC boundary is mapped through BOTH transforms and
measured identically (stats.cell_stats). Paired per-cell ratios overlay
the curves: AR ratios must fall inside the sigma envelope, and the
paired area ratios must trace J — for every system, since pairing
cancels each grid's intrinsic area spread.

One-off analysis for issue #82 — not part of the pipeline. Bindings are
deliberately co-resolved here (a one-env exception to the per-script-env
rule) because the script needs each grid's pre-transform geodetic skin.

Run with:  uv run scripts/conf_check.py
(macOS: x86_64 python for dggal — the justfile's UV_PYTHON pin.)
"""

import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent / 'systems'))

from dggs_compare import config, stats
from dggs_compare.stats import _E   # WGS84 eccentricity (shared constant)

# ----- knobs -------------------------------------------------------------
N_CELLS = 100_000                   # paired cells per system
DISPLAY_N = 8000                    # scatter subsample (stats use all cells)
SEED = 42
SYSTEMS = ('a5', 'isea7h', 'ivea7h', 'hex9')   # geodetic-addressed + control
OUT = Path(__file__).resolve().parent.parent / 'out' / 'conf_check.png'
# -------------------------------------------------------------------------


# ----- the two latitude maps (radians in, radians out) -------------------
def chi(phi):
    """Conformal latitude for WGS84 geodetic latitude (Snyder 1987,
    eq. 3-1). Odd; evaluated on |phi| and mirrored, like authalic_lat."""
    p = np.abs(phi)
    s = np.sin(p)
    x = 2 * np.arctan(np.tan(np.pi / 4 + p / 2)
                      * ((1 - _E * s) / (1 + _E * s)) ** (_E / 2)) - np.pi / 2
    return np.copysign(x, phi)


def xi(phi):
    """Authalic latitude, radians (the library's exact closed form)."""
    return np.radians(stats.authalic_lat(np.degrees(phi)))


def tissot(phi, dphi=1e-7):
    """(sigma, J) of T = chi . xi^{-1} at geodetic latitude phi (radians):
    anisotropy ratio (>= 1) and area scale. Central differences at 1e-7
    rad on the exact closed forms (~1e-9 relative accuracy)."""
    h = ((chi(phi + dphi) - chi(phi - dphi))
         / (xi(phi + dphi) - xi(phi - dphi)))
    k = np.cos(chi(phi)) / np.cos(xi(phi))
    return np.maximum(h / k, k / h), h * k


def conformal_rings(rings):
    """authalic_rings' twin: every vertex latitude through chi."""
    out = []
    for ring in rings:
        r = np.asarray(ring, dtype=float)
        lat = np.degrees(chi(np.radians(r[:, 0])))
        out.append(list(zip(lat.tolist(), r[:, 1].tolist())))
    return out


# ----- per-system raw-geodetic cell sources ------------------------------
def _a5_source():
    import a5_fast as a5
    res = config.TARGET_RES['a5']

    def corners(latlng, turn_deg=5.0):
        """Corner reduction copied from scripts/systems/a5-a5.py (a5_fast
        returns adaptively densified boundaries)."""
        la, lo = np.radians(np.asarray(latlng, dtype=float)).T
        v = np.column_stack([np.cos(la) * np.cos(lo),
                             np.cos(la) * np.sin(lo), np.sin(la)])
        e = np.roll(v, -1, axis=0) - v
        e /= np.linalg.norm(e, axis=1, keepdims=True)
        cosang = np.clip(np.einsum('ij,ij->i', np.roll(e, 1, axis=0), e),
                         -1, 1)
        idx = np.nonzero(np.degrees(np.arccos(cosang)) > turn_deg)[0]
        return [latlng[i] for i in idx] if len(idx) >= 3 else latlng

    def cells(points):
        return [a5.lonlat_to_cell(lng, lat, res) for lat, lng in points]

    def ring(c):
        verts = a5.cell_to_boundary(c)   # closed (lng, lat), densified
        return corners([(lat, lng) for lng, lat in verts[:-1]])

    return res, cells, ring


def _dggal_source(cls):
    from _dggal_engine import Adapter
    ad = Adapter(cls)
    res = config.TARGET_RES[cls.lower()]

    def cells(points):
        return ad.cells_at(res, points)

    return res, cells, ad.cell_boundary


def _hex9_source():
    import hex9
    res = config.TARGET_RES['hex9']

    def cells(points):
        pts = np.asarray(points, dtype=float).reshape(-1, 2)
        full = hex9.encode(np.ascontiguousarray(pts[:, 1]),
                           np.ascontiguousarray(pts[:, 0]))
        return [b.tobytes() for b in hex9.bin(full, res)]

    def ring(c):
        # Open CCW (lat, lng) geodetic ring (see scripts/systems/hex9-hex9.py)
        u = np.frombuffer(c, dtype=np.uint8)
        return hex9.cell(u, res, 0)[-2::-1, ::-1]

    return res, cells, ring


SOURCES = {'a5': _a5_source,
           'isea7h': lambda: _dggal_source('ISEA7H'),
           'ivea7h': lambda: _dggal_source('IVEA7H'),
           'hex9': _hex9_source}


def sample_system(name):
    """N_CELLS distinct point-drawn cells: per-cell (centroid |lat| deg,
    ar_auth, ar_conf, area_auth_norm, area_conf_norm)."""
    res, cells_at, ring_of = SOURCES[name]()
    mean_area = 4.0 * np.pi / config.CELLS_PER_RES[name](res)
    rng = np.random.default_rng(SEED)
    seen, picked = set(), []
    while len(picked) < N_CELLS:
        pts = stats.sample_uniform_latlng(2 * (N_CELLS - len(picked)), rng)
        for c in cells_at(pts):
            if c not in seen:
                seen.add(c)
                picked.append(c)
                if len(picked) == N_CELLS:
                    break
    rows = []
    for c in picked:
        ring = [(float(la), float(ln)) for la, ln in ring_of(c)]
        ar_a, area_a = stats.cell_stats(stats.authalic_rings([ring])[0])
        ar_c, area_c = stats.cell_stats(conformal_rings([ring])[0])
        lat = abs(np.mean([la for la, _ in ring]))
        rows.append((lat, ar_a, ar_c, area_a / mean_area, area_c / mean_area))
    return np.array(rows)


def main():
    # ----- the curves ----------------------------------------------------
    phi = np.radians(np.linspace(0.0, 89.99, 20_000))
    sigma, J = tissot(phi)
    i = int(np.argmax(sigma))
    print(f'analytic: max sigma = {sigma[i]:.6f} at |lat| '
          f'{np.degrees(phi[i]):.2f} deg;  '
          f'J spread max/min = {J.max() / J.min():.6f} '
          f'(J in [{J.min():.6f}, {J.max():.6f}])')

    # ----- the samples ---------------------------------------------------
    data = {}
    for s in SYSTEMS:
        data[s] = d = sample_system(s)
        lat, ar_a, ar_c, na, nc = d.T
        r = np.maximum(ar_c / ar_a, ar_a / ar_c)
        env = np.interp(lat, np.degrees(phi), sigma)
        print(f'{s:7} n {len(d):>7}  AR ratio: median {np.median(r):.6f} '
              f'max {r.max():.6f} (bound {env.max():.6f}, '
              f'violations {(r > env * (1 + 1e-6)).sum()})  '
              f'area_conf max/min {nc.max() / nc.min():.6f} '
              f'(auth {na.max() / na.min():.6f})')

    # ----- the figure ----------------------------------------------------
    # Scatter shows a DISPLAY_N subsample per system (stats above use every
    # cell); the paired area ratio cancels each grid's intrinsic spread, so
    # panel 2's dots should collapse onto J for every system.
    rng = np.random.default_rng(SEED)
    shown = {s: data[s][rng.choice(len(data[s]),
                                   min(DISPLAY_N, len(data[s])),
                                   replace=False)]
             for s in SYSTEMS}
    label = {s: f'{s.upper()} {config.RES_PREFIX[s]}{config.TARGET_RES[s]}'
             for s in SYSTEMS}
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))
    latdeg = np.degrees(phi)
    ax1.plot(latdeg, sigma, color='0.2', lw=1.5, zorder=3,
             label='σ(φ) — analytic envelope')
    for s in SYSTEMS:
        d = shown[s]
        r = np.maximum(d[:, 2] / d[:, 1], d[:, 1] / d[:, 2])
        ax1.plot(d[:, 0], r, '.', ms=2, alpha=0.25,
                 color=matplotlib.colors.to_hex(config.SYS_COLOR[s]),
                 label=label[s])
    ax1.annotate(f'max σ = {sigma.max():.6f} (equator)',
                 (0.03, 0.03), xycoords='axes fraction', fontsize=9)
    ax1.set_xlabel('|latitude| (deg)')
    ax1.set_ylabel(r'per-cell  $\max(\mathrm{AR}_\chi/\mathrm{AR}_\xi,'
                   r'\ \mathrm{AR}_\xi/\mathrm{AR}_\chi)$')
    ax1.set_title('shape: the per-cell AR ratio is bounded by σ(φ)',
                  fontsize=10)
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    ax2.plot(latdeg, J, color='0.2', lw=1.5, zorder=3,
             label='J(φ) — analytic area scale')
    for s in SYSTEMS:
        d = shown[s]
        ax2.plot(d[:, 0], d[:, 4] / d[:, 3], '.', ms=2, alpha=0.25,
                 color=matplotlib.colors.to_hex(config.SYS_COLOR[s]),
                 label=label[s])
    ax2.annotate(f'max/min J = {J.max() / J.min():.6f}',
                 (0.03, 0.93), xycoords='axes fraction', fontsize=9)
    ax2.set_xlabel('|latitude| (deg)')
    ax2.set_ylabel(r'per-cell  $\mathrm{area}_\chi \,/\, \mathrm{area}_\xi$')
    ax2.set_title('area: the per-cell ratio equals J(φ) exactly',
                  fontsize=10)
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    fig.suptitle('Authalic (ξ) vs conformal (χ) spherization: the Tissot '
                 'curves of T = χ∘ξ⁻¹, with paired samples at the working '
                 'resolutions', fontsize=12)
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=200)
    print(f'wrote {OUT}')


if __name__ == '__main__':
    main()
