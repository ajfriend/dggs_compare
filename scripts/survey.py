"""DGGS aspect-ratio + cell-area survey — reads the `ar`, `area`, and
`irregular` columns of the tables and writes the comparison plots to out/
(the only csar calls re-solve each system's two extreme cells, to draw
their certified ellipses):

  histograms.png        cross-system AR distributions at the working resolutions
  area_histograms.png   cross-system cell-area distributions (area / (4π/N),
                        regular cells; pentagon ratio as a labeled tick)
  area_ratio_by_res.png worst-case area ratio (max/min, regular cells) vs
                        resolution — the area twin of the AR curves
  tradeoff.png          worst-case AR vs worst-case area ratio per system
  extremes.png          best/worst cell per system, drawn with its ellipse
  by_res_<sys>.png      per-system AR distribution stacked by resolution

Aspect ratio (AR) = major/minor semi-axis ratio (a/b, a>=1) of each cell's
enclosing-cone ellipse — the discrete, per-cell analogue of Tissot's
indicatrix. AR == 1 is isotropic; AR > 1 is anisotropic ("squished"). For an
equal-area DGGS the areal factor a*b is ~fixed by construction, so AR is
where the unavoidable distortion shows up.

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
# -------------------------------------------------------------------------


def sweep_system(name):
    """Per-resolution AR arrays, DNC counts, and area-ratio scalars from
    the tables, plus the best/worst cell at the target resolution
    (re-solved — two cells — so the extremes plot can draw the certified
    ellipse). Full per-cell area arrays are kept ONLY at the target
    resolution (the histogram input): retained everywhere they'd cost
    ~2.5 GB at the full budget for values only ever reduced to scalars."""
    target = RES[name]
    by_res = {}
    for res in cache.available_resolutions(name):
        cols = cache.load_columns(name, res, ['ar', 'area', 'irregular'])
        a = cols['ar']
        rel = cols['area'] / config.mean_cell_area(name, res)
        irr = cols['irregular']
        d = {'ars': a[~np.isnan(a)], 'dnc': int(np.isnan(a).sum()),
             'ratio': area_ratio(rel, irr),
             # The pentagons' relative areas, kept only where both cell
             # classes are present (empty for all-pentagon coarse levels
             # and for samples that caught none).
             'pent': rel[irr] if irr.any() and not irr.all() else rel[:0]}
        if res == target:
            d['area'], d['irr'] = rel, irr
        by_res[res] = d

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


def regular_areas(area, irr):
    """The REGULAR cells' normalized areas — the pentagon deficit is a
    design constant (exactly 5/6 for the ISEA family), not distortion,
    so the area statistics exclude it. A resolution with no regular
    cells (ISEA-family r0: all 12 cells ARE the pentagons) uses every
    cell. THE selection rule: every area stat and histogram goes through
    here."""
    a = area[~irr]
    return a if a.size else area


def area_ratio(area, irr):
    """max/min normalized area over regular cells (see regular_areas)."""
    a = regular_areas(area, irr)
    return float(a.max() / a.min())


def _save(fig, name, dpi=DPI):
    out = OUT_DIR / name
    fig.savefig(out, dpi=dpi)
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
    _save(fig, 'histograms.png')


def plot_area_histograms(results):
    """Cross-system normalized-area distributions at the working
    resolutions, regular cells only. The pentagon tick is sourced from
    the deepest resolution that carries pentagons — 12 cells among
    millions are ~never sampled at the working resolution."""
    data = {}
    for s in SYSTEMS:
        by_res = results[s]['by_res']
        d = by_res[RES[s]]
        reg = regular_areas(d['area'], d['irr'])
        pent_res = next((r for r in sorted(by_res, reverse=True)
                         if by_res[r]['pent'].size), None)
        data[s] = {'reg': reg, 'med': float(np.median(reg)),
                   'ratio': d['ratio'],
                   'all_ratio': float(d['area'].max() / d['area'].min()),
                   'pent_res': pent_res,
                   'pent': (by_res[pent_res]['pent']
                            if pent_res is not None else reg[:0])}

    print(f'{"sys":5} {"n":>8} {"min":>10} {"median":>10} {"max":>10} '
          f'{"max/min":>9} {"(all)":>9}')
    for s in SYSTEMS:
        v = data[s]
        print(f'{s:5} {v["reg"].size:>8} {v["reg"].min():>10.6f} '
              f'{v["med"]:>10.6f} {v["reg"].max():>10.6f} '
              f'{v["ratio"]:>9.6f} {v["all_ratio"]:>9.6f}')

    lo = min(a.min() for v in data.values()
             for a in (v['reg'], v['pent']) if a.size)
    hi = max(v['reg'].max() for v in data.values())
    # Pad so 1.0 stays in frame even when every system is a spike at 1.
    lo, hi = min(lo, 0.98), max(hi, 1.02)
    bins = np.linspace(lo, hi, N_BINS + 1)
    fig, axes = plt.subplots(len(SYSTEMS), 1, figsize=(8, 9), sharex=True)
    for ax, s in zip(axes, SYSTEMS):
        v = data[s]
        ax.hist(v['reg'], bins=bins, color=SYS_COLOR[s], edgecolor='white',
                linewidth=0.3)
        ax.set_yscale('log')
        ax.set_ylabel('count (log)')
        if v['pent'].size:
            p = float(v['pent'].mean())
            ax.axvline(p, color='0.25', ls='--', lw=1.2)
            ax.text(p, 0.92, f' pentagons {p:.4f} (r{v["pent_res"]})',
                    transform=ax.get_xaxis_transform(),
                    fontsize=8, color='0.25', va='top')
        ax.set_title(f'{SYS_LABEL[s]}  (median {v["med"]:.4f}, '
                     f'max/min {v["ratio"]:.4f})', fontsize=10)
        ax.grid(True, alpha=0.3)
    axes[-1].set_xlabel('cell area / (4π/N) — exact-mean normalized; '
                        'regular cells (shared bins)')
    fig.suptitle('DGGS cell-area distributions (~H3 r9 cell size)', fontsize=12)
    fig.tight_layout()
    _save(fig, 'area_histograms.png')


def plot_area_ratio_by_res(results):
    """One figure: worst-case area ratio vs resolution per system — the
    area twin of the worst-case-AR-by-resolution story. Equal-area
    systems hold 1.0 to float precision; sampled resolutions use sample
    extremes (area is a smooth field, so 1M uniform samples reach close
    to the true range)."""
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for s in SYSTEMS:
        by_res = results[s]['by_res']
        res_list = sorted(by_res)
        ax.plot(res_list, [by_res[r]['ratio'] for r in res_list],
                '-o', ms=3.5, lw=1.4, color=SYS_COLOR[s], label=SYS_LABEL[s])
    ax.axhline(1.0, color='0.8', lw=1, zorder=0)
    ax.set_xlabel('resolution')
    ax.set_ylabel('max/min cell area (regular cells, log)')
    # Log y: the equal-area family hugs 1.0 while the finest levels of
    # sub-attosteradian cells can spike on the f64 area floor — log keeps
    # the 2x band readable either way.
    ax.set_yscale('log')
    ax.yaxis.set_major_formatter('{x:g}')
    ax.set_yticks([1.0, 1.1, 1.25, 1.5, 2.0, 3.0, 4.0])
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_title('Worst-case area ratio by resolution\n'
                 '(pentagons excluded — their 5/6 deficit is a design '
                 'constant, not distortion)', fontsize=11)
    fig.tight_layout()
    _save(fig, 'area_ratio_by_res.png')


def plot_tradeoff(results):
    """Pareto scatter at the working resolutions: worst-case AR (all
    cells, matching the published AR stats) against worst-case area
    ratio (regular cells). The ideal grid sits at (1, 1)."""
    fig, ax = plt.subplots(figsize=(7.5, 6))
    for s in SYSTEMS:
        d = results[s]['by_res'][RES[s]]
        x, y = d['ars'].max(), d['ratio']
        ax.plot(x, y, 'o', ms=9, color=SYS_COLOR[s])
        ax.annotate(SYS_LABEL[s], (x, y), textcoords='offset points',
                    xytext=(7, 4), fontsize=9)
    ax.axhline(1.0, color='0.85', lw=1, zorder=0)
    ax.axvline(1.0, color='0.85', lw=1, zorder=0)
    ax.set_xlabel('worst-case aspect ratio (all cells)')
    ax.set_ylabel('worst-case area ratio, max/min (regular cells)')
    ax.grid(True, alpha=0.3)
    ax.set_title('Shape vs area: worst case at the working resolutions\n'
                 '(ideal = (1, 1); AR keeps pentagons, area excludes them)',
                 fontsize=11)
    fig.tight_layout()
    _save(fig, 'tradeoff.png')


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
    _save(fig, 'extremes.png')
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
    _save(fig, f'by_res_{name}.png', dpi=130)


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
    plot_area_histograms(results)
    plot_area_ratio_by_res(results)
    plot_tradeoff(results)
    plot_extremes(results)
    for s in SYSTEMS:
        plot_by_resolution(s, results[s]['by_res'])


if __name__ == '__main__':
    main()
