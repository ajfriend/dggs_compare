"""DGGS aspect-ratio + authalicity survey — reads the `ar` and `area`
columns of the tables (no solving) and writes the comparison plots to out/:

  histograms.png     cross-system AR distributions at the working resolutions
  extremes.png       best/worst cell per system, drawn with its ellipse
  by_res_<sys>.png   per-system AR distribution stacked by resolution
  authalicity.png    cross-system cell-area ratio distributions

Aspect ratio (AR) = major/minor semi-axis ratio (a/b, a>=1) of each cell's
enclosing-cone ellipse — the discrete, per-cell analogue of Tissot's
indicatrix. AR == 1 is isotropic; AR > 1 is anisotropic ("squished"). For an
equal-area DGGS the areal factor a*b is ~fixed by construction, so AR is
where the unavoidable distortion shows up.

Authalicity: each cell's area as a ratio to the ideal equal-area share,
area * N(res) / 4pi, with N(res) from the closed-form counts. Ratio 1 is
perfect equal-area; the spread (CV%, max/min) is the areal counterpart of
the AR story. Two honesty notes, both handled below: (1) the tables' `area`
is measured on the coordinate sphere (geodetic lat/lng read as spherical) —
a shared datum, but the systems that address the WGS84 authalic sphere or
ellipsoid internally pick up ~0.1-0.4% of apparent CV that is datum
mismatch, not grid distortion, so a native-datum CV (geodetic->authalic
latitudes, exact closed form) is reported alongside; (2) point-sampled
tables include cells with probability ~proportional to area, so their raw
mean ratio sits at ~1+CV^2 — sampled-regime stats are inverse-area weighted
back to per-cell uniform.

Run with:  just survey
No CLI args (project convention).
"""

import time
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

import csar

from dggs_compare import cache, config

# ----- knobs -------------------------------------------------------------
RES = config.TARGET_RES
SYSTEMS = [s for s in RES if s in cache.available_systems()]
SYS_LABEL = {s: f'{s.upper()} {config.RES_PREFIX[s]}{RES[s]}' for s in SYSTEMS}
SYS_COLOR = config.SYS_COLOR
OUT_DIR = Path(__file__).resolve().parent.parent / 'out'
N_BINS = 60
DPI = 200
# Systems whose native addressing surface IS the coordinate sphere. Everyone
# else converts geodetic latitude internally and does its equal-area work on
# the WGS84 authalic sphere (dggal's Snyder family, DGGRID's ISEA4T, hex9)
# or the ellipsoid (a5) — for those the native-datum column re-reads the
# boundary through the exact authalic latitude, which also covers the
# ellipsoid case (the authalic sphere preserves ellipsoidal areas, and the
# ratio normalization cancels the radius).
SPHERE_NATIVE = frozenset({'h3', 's2'})
# -------------------------------------------------------------------------


def sweep_system(name):
    """Per-resolution AR arrays + DNC counts from the tables, plus the
    best/worst cell at the target resolution (re-solved — two cells — so the
    extremes plot can draw the certified ellipse)."""
    target = RES[name]
    by_res = {}
    for res in cache.available_resolutions(name):
        a = cache.load_columns(name, res, ['ar'])['ar']
        by_res[res] = {'ars': a[~np.isnan(a)], 'dnc': int(np.isnan(a).sum())}

    cols = cache.load_columns(name, target, ['cid', 'verts', 'ar'])

    def record(idx):
        verts = csar.to_vec3(cols['verts'][idx], geo='latlng_deg')
        r = csar.solve(verts, geo='vec3', gap_tol=config.GAP_TOL,
                       method=config.CSAR_METHOD)
        return {'ar': float(cols['ar'][idx]), 'id': cols['cid'][idx],
                'verts': verts, 'result': r}

    return {'by_res': by_res,
            'best': record(int(np.nanargmin(cols['ar']))),
            'worst': record(int(np.nanargmax(cols['ar'])))}


# ----- authalicity -------------------------------------------------------
_E2 = 0.00669437999014133                # WGS84 first eccentricity squared


def _authalic_lat(lat_deg):
    """Exact authalic latitude (degrees), closed-form q-function."""
    e = np.sqrt(_E2)

    def q(s):
        return (1 - _E2) * (s / (1 - _E2 * s * s)
                            - (0.5 / e) * np.log((1 - e * s) / (1 + e * s)))

    s = np.sin(np.radians(lat_deg))
    return np.degrees(np.arcsin(np.clip(q(s) / q(1.0), -1.0, 1.0)))


def _chord_areas(verts_list, authalic):
    """Cell areas (steradians) from corner rings, vectorized by ring length.

    Planar 3D-polygon area 0.5*|sum(w_i x w_i+1)| with vertices centred on
    the ring mean — at working-resolution cell sizes this matches the
    spherical area to ~1e-15 relative while staying at full precision (the
    centring keeps the edge vectors O(cell size), not O(1)). With
    `authalic`, latitudes are mapped geodetic->authalic first — the sphere
    the non-SPHERE_NATIVE systems actually address.
    """
    out = np.empty(len(verts_list))
    by_len = {}
    for i, v in enumerate(verts_list):
        by_len.setdefault(len(v), []).append(i)
    for _, idx in by_len.items():
        v = np.asarray([verts_list[i] for i in idx], float)   # (m, k, 2)
        lat = _authalic_lat(v[..., 0]) if authalic else v[..., 0]
        la, lo = np.radians(lat), np.radians(v[..., 1])
        xyz = np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo),
                        np.sin(la)], axis=-1)
        w = xyz - xyz.mean(axis=1, keepdims=True)
        cr = np.cross(w, np.roll(w, -1, axis=1)).sum(axis=1)
        out[np.asarray(idx)] = 0.5 * np.linalg.norm(cr, axis=1)
    return out


def sweep_area(name):
    """The target-resolution cell-area ratios (area / ideal equal-area share)
    plus per-cell weights and the native-datum CV.

    Weights: the point-draw selection regime includes a cell with probability
    ~proportional to its area, so those tables are inverse-area weighted back
    to per-cell uniform (raw mean sits at ~1+CV^2 otherwise). Enumerated and
    subsampled tables (total <= SUBSAMPLE_MAX_RATIO * N_CELLS) are already
    uniform — the regime is derivable from the closed-form count.
    """
    res = RES[name]
    total = config.CELLS_PER_RES[name](res)
    cols = cache.load_columns(name, res, ['area', 'verts'])
    ratio = cols['area'] * total / (4.0 * np.pi)
    sampled = total > config.SUBSAMPLE_MAX_RATIO * config.N_CELLS
    w = 1.0 / ratio if sampled else np.ones_like(ratio)
    native_cv = None
    if name not in SPHERE_NATIVE:
        nat = _chord_areas(cols['verts'], authalic=True)
        nr = nat * total / (4.0 * np.pi)
        mu = np.average(nr, weights=w)
        native_cv = 100.0 * np.sqrt(np.average((nr - mu) ** 2, weights=w)) / mu
    return {'ratio': ratio, 'w': w, 'sampled': sampled, 'native_cv': native_cv}


def plot_authalicity(area_results):
    """Cross-system cell-area-ratio distributions at the working resolutions,
    stacked on shared log2 bins (the spreads differ by ~100x across systems,
    so linear shared bins would flatten the tight ones)."""
    lo = min(float(d['ratio'].min()) for d in area_results.values())
    hi = max(float(d['ratio'].max()) for d in area_results.values())
    bins = np.geomspace(lo * 0.999, hi * 1.001, N_BINS + 1)

    print(f'\n{"sys":5} {"n":>8} {"mean":>8} {"cv%":>8} {"natCV%":>8} '
          f'{"min":>8} {"median":>8} {"max":>8} {"max/min":>8}')
    for s in SYSTEMS:
        d = area_results[s]
        r, w = d['ratio'], d['w']
        mu = np.average(r, weights=w)
        cv = 100.0 * np.sqrt(np.average((r - mu) ** 2, weights=w)) / mu
        nat = f'{d["native_cv"]:>7.4f}†' if d['native_cv'] is not None else f'{cv:>7.4f} '
        print(f'{s:5} {r.size:>8} {mu:>8.5f} {cv:>8.4f} {nat} '
              f'{r.min():>8.4f} {np.median(r):>8.4f} {r.max():>8.4f} '
              f'{r.max() / r.min():>8.4f}')
    print('† native-datum CV: boundary re-read through the exact authalic '
          'latitude —\n  the surface these grids address internally; the '
          'plain cv% keeps the shared\n  coordinate-sphere datum for '
          'comparability.')

    fig, axes = plt.subplots(len(SYSTEMS), 1, figsize=(8, 9), sharex=True)
    for ax, s in zip(axes, SYSTEMS):
        d = area_results[s]
        r, w = d['ratio'], d['w']
        mu = np.average(r, weights=w)
        cv = 100.0 * np.sqrt(np.average((r - mu) ** 2, weights=w)) / mu
        ax.hist(r, bins=bins, weights=w, color=SYS_COLOR[s],
                edgecolor='white', linewidth=0.3)
        ax.axvline(1.0, color='0.25', lw=0.8, ls='--', alpha=0.7)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylabel('count (log)')
        nat = (f', native CV {d["native_cv"]:.3f}%'
               if d['native_cv'] is not None else '')
        ax.set_title(f'{SYS_LABEL[s]}  (CV {cv:.3f}%{nat}, '
                     f'max/min {r.max() / r.min():.4f})', fontsize=10)
        ax.grid(True, alpha=0.3, which='both')
    ticks = [t for t in (0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4)
             if lo * 0.99 <= t <= hi * 1.01]
    axes[-1].set_xticks(ticks, [f'{t:g}' for t in ticks])
    axes[-1].xaxis.set_minor_formatter(NullFormatter())
    axes[-1].set_xlabel('cell area / ideal equal-area share  '
                        '(shared log bins; 1.0 = perfect)')
    fig.suptitle('DGGS cell-area ratio distributions (~H3 r9 cell size)',
                 fontsize=12)
    fig.tight_layout()
    out = OUT_DIR / 'authalicity.png'
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print(f'wrote {out}')


# ----- plotting ----------------------------------------------------------
def plot_histograms(results):
    """Cross-system AR distributions at each system's working resolution."""
    ars = {s: results[s]['by_res'][RES[s]]['ars'] for s in SYSTEMS}
    dnc = {s: results[s]['by_res'][RES[s]]['dnc'] for s in SYSTEMS}

    def stat(a):
        return dict(n=a.size, min=float(a.min()), median=float(np.median(a)),
                    p99=float(np.percentile(a, 99)), max=float(a.max()))

    st = {s: stat(ars[s]) for s in SYSTEMS}

    print(f'{"sys":5} {"n_conv":>8} {"n_dnc":>7} {"min":>10} {"median":>10} {"p99":>10} {"max":>10}')
    for s in SYSTEMS:
        d = st[s]
        print(f'{s:5} {d["n"]:>8} {dnc[s]:>7} {d["min"]:>10.6f} '
              f'{d["median"]:>10.6f} {d["p99"]:>10.6f} {d["max"]:>10.6f}')

    bins = np.linspace(1.0, max(a.max() for a in ars.values()), N_BINS + 1)
    fig, axes = plt.subplots(len(SYSTEMS), 1, figsize=(8, 9), sharex=True)
    for ax, s in zip(axes, SYSTEMS):
        d = st[s]
        ax.hist(ars[s], bins=bins, color=SYS_COLOR[s], edgecolor='white', linewidth=0.3)
        ax.set_yscale('log')
        ax.set_ylabel('count (log)')
        ax.set_title(f'{SYS_LABEL[s]}  (median {d["median"]:.4f}, max {d["max"]:.4f}, DNC {dnc[s]})',
                     fontsize=10)
        ax.grid(True, alpha=0.3)
    axes[-1].set_xlabel(f'aspect ratio (shared bins, gap_tol = {config.GAP_TOL:g})')
    fig.suptitle('DGGS aspect-ratio distributions (~H3 r9 cell size)', fontsize=12)
    fig.tight_layout()
    out = OUT_DIR / 'histograms.png'
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print(f'wrote {out}')


def draw_cell(ax, rec, color):
    """Draw a cell's boundary + enclosing ellipse, major axis horizontal."""
    xy, semi = csar.project_to_cone(rec['result'], rec['verts'], up=None)
    ring = np.vstack([xy, xy[:1]])
    t = np.linspace(0.0, 2.0 * np.pi, 400)
    ax.plot(ring[:, 0], ring[:, 1], '-o', color=color, lw=1.3, ms=4, label='cell')
    ax.plot(semi[0] * np.cos(t), semi[1] * np.sin(t), '-', color='0.25', lw=1.5,
            label='enclosing ellipse')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.text(0.03, 0.95, f'AR {rec["ar"]:.4f}\nid {rec["id"]}', transform=ax.transAxes,
            va='top', ha='left', fontsize=8,
            bbox=dict(boxstyle='round', fc='white', ec='0.7', alpha=0.85))


def plot_extremes(results):
    fig, axes = plt.subplots(len(SYSTEMS), 2, figsize=(11, 4.3 * len(SYSTEMS)))
    axes[0, 0].set_title('best AR (most circular)', fontsize=12, pad=10)
    axes[0, 1].set_title('worst AR', fontsize=12, pad=10)
    for row, s in enumerate(SYSTEMS):
        for col, kind in ((0, 'best'), (1, 'worst')):
            ax = axes[row, col]
            draw_cell(ax, results[s][kind], SYS_COLOR[s])
            ax.set_xlabel('major axis (m)')
            if col == 0:
                ax.set_ylabel(f'{SYS_LABEL[s]}\nminor axis (m)')
    axes[0, 0].legend(loc='lower right', fontsize=8)
    fig.suptitle('DGGS cells (~H3 r9 cell size): best vs worst aspect ratio\n'
                 '(enclosing-cone cross-section ‖Ax‖ <= b·x; major axis horizontal)',
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = OUT_DIR / 'extremes.png'
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print(f'wrote {out}')
    for s in SYSTEMS:
        print(f'  {s}: best AR {results[s]["best"]["ar"]:.4f}  ·  '
              f'worst AR {results[s]["worst"]["ar"]:.4f}')


def plot_by_resolution(name, by_res):
    """One tall file per system: AR distribution at every resolution stacked
    vertically on a shared aspect-ratio axis. Resolutions with DNC failures
    are labelled in red."""
    res_list = sorted(by_res)
    allars = np.concatenate([d['ars'] for d in by_res.values() if d['ars'].size])
    amax = float(np.percentile(allars, 99.9))
    bins = np.linspace(1.0, amax, N_BINS + 1)
    # Put the n/DNC note on the emptier horizontal half.
    mass = np.histogram(allars, bins=bins)[0]
    third = max(len(mass) // 3, 1)
    note_x, note_ha = ((0.985, 'right') if mass[:third].sum() >= mass[-third:].sum()
                       else (0.015, 'left'))

    n = len(res_list)
    fig_h = 1.4 * n + 1.4
    fig, axes = plt.subplots(n, 1, figsize=(10, fig_h), sharex=True, squeeze=False)
    for ax, res in zip(axes[:, 0], res_list):
        d = by_res[res]
        a, dnc = d['ars'], d['dnc']
        red = bool(dnc)
        counts = (ax.hist(a, bins=bins, color=SYS_COLOR[name],
                          edgecolor='white', linewidth=0.3)[0]
                  if a.size else np.zeros(1))
        ax.set_yscale('log')
        # Pin major y ticks to exact powers of ten; sub-decade ticks are clutter.
        maxc = max(int(counts.max()), 1)
        ax.set_yticks([10.0 ** k for k in range(int(np.floor(np.log10(maxc))) + 1)])
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.set_ylabel(f'r{res}', rotation=0, ha='right', va='center', labelpad=12,
                      fontsize=13, fontweight='bold', color='red' if red else '0.2')
        ax.tick_params(labelsize=10)
        ax.grid(True, alpha=0.25)
        note = f'n = {a.size:,}' + (f'      DNC {dnc:,}' if red else '')
        ax.text(note_x, 0.9, note, transform=ax.transAxes, ha=note_ha, va='top',
                fontsize=11, color='red' if red else '0.4')
    axes[-1, 0].set_xlabel(f'aspect ratio (shared bins, gap_tol = {config.GAP_TOL:g})',
                           fontsize=12)
    fig.suptitle(f'{name.upper()} aspect-ratio distribution by resolution '
                 f'(coarsest at top; shared bins 1.00–{amax:.2f}, log y)',
                 fontsize=15, y=1 - 0.4 / fig_h, va='top')
    fig.tight_layout(rect=(0, 0, 1, 1 - 1.0 / fig_h))
    out = OUT_DIR / f'by_res_{name}.png'
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f'wrote {out}')


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results = {}
    for s in SYSTEMS:
        t0 = time.perf_counter()
        results[s] = sweep_system(s)
        nres = len(results[s]['by_res'])
        ncells = sum(d['ars'].size + d['dnc'] for d in results[s]['by_res'].values())
        print(f'[{s}] {ncells:,} cells over {nres} resolutions '
              f'in {time.perf_counter() - t0:.1f}s')

    plot_histograms(results)
    plot_extremes(results)
    for s in SYSTEMS:
        plot_by_resolution(s, results[s]['by_res'])

    area_results = {}
    for s in SYSTEMS:
        t0 = time.perf_counter()
        area_results[s] = sweep_area(s)
        print(f'[{s}] area sweep in {time.perf_counter() - t0:.1f}s')
    plot_authalicity(area_results)


if __name__ == '__main__':
    main()
