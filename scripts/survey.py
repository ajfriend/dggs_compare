"""DGGS aspect-ratio + cell-area survey — reads the `ar` and `area`
columns of the tables and writes the comparison plots to out/ (the only
csar calls re-solve each system's two extreme cells, to draw their
certified ellipses):

  histograms.png        cross-system AR distributions at the working resolutions
  area_histograms.png   cross-system cell-area distributions (area / (4π/N),
                        all cells)
  tradeoff.png          shape-vs-area scatter, all-cells and central-99.9% panels
  extremes.png          best/worst cell per system, drawn with its ellipse
  by_res_<sys>.png      per-system AR distribution stacked by resolution
  summary_table.html    the working-resolution stat lines as one table
                        (an HTML fragment; injected into the page by globe.js)

Aspect ratio (AR) = major/minor semi-axis ratio (a/b, a>=1) of each cell's
enclosing-cone ellipse — the discrete, per-cell analogue of Tissot's
indicatrix. AR == 1 is isotropic; AR > 1 is anisotropic (elongated).

Two summary statistics per metric (#83): the ALL-CELLS extreme (max, or
max/min for area — the worst case an application can meet) and the
CENTRAL-99.9% version (p_9995, or p_9995/p_0005 — one uniform trim,
identical for every system, no declared cell classes). The 0.05%
per-side cut sits above the vanishing defect classes (point defects
~1/N, seam lines ~1/sqrt(N)) and below the bulk, and a 0.05% quantile
of a 1M-row sample is statistically stable where a sample extreme
is not.

Every distribution and quantile is PER-CELL UNIFORM: point-draw-sampled
tables include cells ~proportionally to their area, so they are
inverse-area weighted back (see sweep_system). All-cells extremes are
weighting-free.

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
P_LO, P_HI = 0.0005, 0.9995   # the central-99.9% trim (#83)
# -------------------------------------------------------------------------


def sweep_system(name):
    """Per-resolution AR arrays and DNC counts from the tables, plus
    the best/worst cell at the target resolution
    (re-solved — two cells — so the extremes plot can draw the certified
    ellipse). Full per-cell area arrays are kept ONLY at the target
    resolution (the area-histogram and tradeoff input): retained
    everywhere they'd cost ~2.5 GB at the full budget when their only
    remaining job is deriving the sampling weights. Weights are the
    exception — the by-resolution histograms consume them at every
    sampled resolution — kept lean as float32 and None where uniform
    (~0.6 GB at the full budget)."""
    target = RES[name]
    impl = config.PRIMARY_IMPL[name]
    by_res = {}
    for res in cache.available_resolutions(name):
        cols = cache.load_columns(name, res, ['ar', 'area'])
        a = cols['ar']
        rel = cols['area'] / config.mean_cell_area(name, res)
        # Point-draw-sampled tables include a cell with probability
        # ~proportional to its area; inverse-area weights (None = the
        # table is already per-cell uniform; float32 — histogram bars and
        # quantile crossings don't need better) restore the per-cell
        # measure every distribution here reports. The raw size-biased
        # mean of rel sits at ~1+CV^2. Classified at the table's own
        # STAMPED budget, so override-budget builds weight correctly too.
        budget = int(cache.table_metadata(name, impl, res)['n_cells'])
        sampled = config.sampling_regime(name, res, budget) == 'sample'
        w = (1.0 / rel).astype(np.float32) if sampled else None
        bad = np.isnan(a)
        d = {'ars': a[~bad], 'w_ars': None if w is None else w[~bad],
             'dnc': int(bad.sum())}
        if res == target:
            d['area'], d['w'] = rel, w
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


def _wq(vals, q, w):
    """Weighted quantile(s), per-cell uniform (inverted_cdf is the only
    method numpy allows with weights; w=None on full enumerations)."""
    return np.quantile(vals, q, weights=w, method='inverted_cdf')


def _save(fig, name, dpi=DPI):
    out = OUT_DIR / name
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    print(f'wrote {out}')


# ----- plotting ----------------------------------------------------------
def plot_histograms(results, rows):
    """Cross-system AR distributions at each system's working resolution;
    the stat lines come from the shared summary rows."""
    ars = {s: results[s]['by_res'][RES[s]]['ars'] for s in SYSTEMS}
    ws = {s: results[s]['by_res'][RES[s]]['w_ars'] for s in SYSTEMS}

    print(f'{"sys":5} {"n_conv":>8} {"n_dnc":>7} {"min":>10} {"median":>10} {"p99.95":>10} {"max":>10}')
    for r in rows:
        a = ars[r['sys']]
        print(f'{r["sys"]:5} {a.size:>8} {r["dnc"]:>7} {a.min():>10.6f} '
              f'{r["ar_median"]:>10.6f} {r["ar_p9995"]:>10.6f} {r["ar_max"]:>10.6f}')

    bins = np.linspace(1.0, max(a.max() for a in ars.values()), N_BINS + 1)
    fig, axes = plt.subplots(len(SYSTEMS), 1, figsize=(8, 9), sharex=True)
    for ax, r in zip(axes, rows):
        s = r['sys']
        ax.hist(ars[s], bins=bins, weights=ws[s], color=SYS_COLOR[s],
                edgecolor='white', linewidth=0.3)
        ax.set_yscale('log')
        ax.set_ylabel('count (log)')
        ax.set_title(f'{SYS_LABEL[s]}  (median {r["ar_median"]:.4f}, '
                     f'p99.95 {r["ar_p9995"]:.4f}, max {r["ar_max"]:.4f}, '
                     f'DNC {r["dnc"]})', fontsize=10)
        ax.grid(True, alpha=0.3)
    axes[-1].set_xlabel(f'aspect ratio (shared bins, gap_tol = {config.GAP_TOL:g})')
    fig.suptitle('DGGS aspect-ratio distributions (~H3 r9 cell size)', fontsize=12)
    fig.tight_layout()
    _save(fig, 'histograms.png')


def plot_area_histograms(results, rows):
    """Cross-system normalized-area distributions at the working
    resolutions — every cell included. The stat line pairs the all-cells
    max/min with the central-99.9% ratio (#83), from the shared summary
    rows."""
    areas = {s: results[s]['by_res'][RES[s]]['area'] for s in SYSTEMS}
    ws = {s: results[s]['by_res'][RES[s]]['w'] for s in SYSTEMS}

    print(f'{"sys":5} {"n":>8} {"min":>10} {"median":>10} {"max":>10} '
          f'{"max/min":>9} {"central":>9}')
    for r in rows:
        area = areas[r['sys']]
        print(f'{r["sys"]:5} {area.size:>8} {area.min():>10.6f} '
              f'{r["area_median"]:>10.6f} {area.max():>10.6f} '
              f'{r["area_all"]:>9.6f} {r["area_trim"]:>9.6f}')

    lo = min(a.min() for a in areas.values())
    hi = max(a.max() for a in areas.values())
    # Pad so 1.0 stays in frame even when every system is a spike at 1.
    lo, hi = min(lo, 0.98), max(hi, 1.02)
    bins = np.linspace(lo, hi, N_BINS + 1)
    fig, axes = plt.subplots(len(SYSTEMS), 1, figsize=(8, 9), sharex=True)
    for ax, r in zip(axes, rows):
        s = r['sys']
        ax.hist(areas[s], bins=bins, weights=ws[s], color=SYS_COLOR[s],
                edgecolor='white', linewidth=0.3)
        ax.set_yscale('log')
        ax.set_ylabel('count (log)')
        ax.set_title(f'{SYS_LABEL[s]}  (median {r["area_median"]:.4f}, '
                     f'max/min {r["area_all"]:.4f}, '
                     f'central-99.9% {r["area_trim"]:.4f})', fontsize=10)
        ax.grid(True, alpha=0.3)
    axes[-1].set_xlabel('cell area / (4π/N) — exact-mean normalized; '
                        'all cells (shared bins)')
    fig.suptitle('DGGS cell-area distributions (~H3 r9 cell size)', fontsize=12)
    fig.tight_layout()
    _save(fig, 'area_histograms.png')


def _scatter_panel(ax, pts):
    """One tradeoff panel: colored points, stacked labels for
    near-coincident points (within 3% of the axis range of a cluster's
    first point; the stack anchors there so near-misses can't smear it),
    and the ideal star at (1, 1)."""
    for s, x, y in pts:
        ax.plot(x, y, 'o', ms=9, color=SYS_COLOR[s], label=SYS_LABEL[s])
    xs, ys = [p[1] for p in pts], [p[2] for p in pts]
    ex = 0.03 * ((max(xs) - min(xs)) or 1)
    ey = 0.03 * ((max(ys) - min(ys)) or 1)
    clusters = []
    for pt in pts:
        for c in clusters:
            if abs(pt[1] - c[0][1]) < ex and abs(pt[2] - c[0][2]) < ey:
                c.append(pt)
                break
        else:
            clusters.append([pt])
    for c in clusters:
        x0, y0 = c[0][1], c[0][2]
        for i, (s, _, _) in enumerate(c):
            ax.annotate(s.upper(), (x0, y0), textcoords='offset points',
                        xytext=(7, 4 + 12 * i), fontsize=9)
    ax.plot(1, 1, marker='*', ms=15, color='0.25', zorder=3)
    ax.annotate('ideal', (1, 1), textcoords='offset points',
                xytext=(8, -4), fontsize=9, color='0.25')
    ax.axhline(1.0, color='0.85', lw=1, zorder=0)
    ax.axvline(1.0, color='0.85', lw=1, zorder=0)
    ax.grid(True, alpha=0.3)


def plot_tradeoff(rows):
    """Pareto scatter at the working resolutions, shown both ways side
    by side (#83): all-cells extremes and the central-99.9% quantile
    statistics — the same construction in each panel, neither favored.
    Consumes the summary_stats rows. The ideal grid sits at (1, 1),
    starred."""
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.8),
                             sharex=True, sharey=True)
    panels = (('ar_max', 'area_all',
               'all cells:  max(AR)  vs  max/min(area)'),
              ('ar_p9995', 'area_trim',
               'central 99.9%:  p99.95(AR)  vs  p99.95/p0.05(area)'))
    for ax, (kx, ky, sub) in zip(axes, panels):
        _scatter_panel(ax, [(r['sys'], r[kx], r[ky]) for r in rows])
        ax.set_title(sub, fontsize=10)
        ax.set_xlabel('aspect ratio')
    axes[0].set_ylabel('cell-area ratio')
    axes[0].legend(fontsize=8, ncols=2, framealpha=0.9)   # system -> working res
    fig.suptitle('Shape vs area at the working resolutions (★ = the ideal grid)',
                 fontsize=12)
    fig.tight_layout()
    _save(fig, 'tradeoff.png')


def summary_stats(results):
    """The working-resolution stats, one row per system: for each metric
    the all-cells extreme and its central-99.9% counterpart (#83), plus
    the medians and the DNC count. The single computation of these
    numbers — the histogram stat lines, the tradeoff panels, and the
    summary table all consume the rows, so they agree by construction.
    (AR's floor is 1, so its central-99.9% counterpart is just the
    upper edge p99.95.)"""
    rows = []
    for s in SYSTEMS:
        d = results[s]['by_res'][RES[s]]
        a, w = d['ars'], d['w_ars']
        area, wa = d['area'], d['w']
        med, hi = _wq(a, [0.5, P_HI], w)
        alo, amed, ahi = _wq(area, [P_LO, 0.5, P_HI], wa)
        rows.append({
            'sys': s, 'n': area.size, 'dnc': d['dnc'],
            'ar_median': float(med),
            'ar_p9995': float(hi),
            'ar_max': float(a.max()),
            'area_median': float(amed),
            'area_trim': float(ahi / alo),
            'area_all': float(area.max() / area.min()),
        })
    return rows


def write_summary_table(rows):
    """The summary stats as an HTML fragment, styled by the site's own
    CSS (table.stats in web/style.css)."""
    def sci(n):   # 1000000 -> '1.0e6'
        m, e = f'{n:.1e}'.split('e')
        return f'{m}e{int(e)}'

    body = []
    for r in rows:
        s = r['sys']
        cells = [sci(r['n']),
                 *(f'{r[k]:.4f}' for k in ('ar_median', 'ar_p9995', 'ar_max')),
                 str(r['dnc']),
                 *(f'{r[k]:.4f}' for k in ('area_trim', 'area_all'))]
        body.append(
            '<tr><td class="sys"><span class="swatch" '
            f'style="background:{matplotlib.colors.to_hex(SYS_COLOR[s])}">'
            f'</span>{SYS_LABEL[s]}</td>'
            + ''.join(f'<td class="num">{c}</td>' for c in cells) + '</tr>')
    html = (
        '<table class="stats">\n<thead>\n'
        '<tr><th rowspan="2" class="sys">system</th>'
        '<th rowspan="2" class="num">cells</th>'
        '<th colspan="4" class="grp">aspect ratio</th>'
        '<th colspan="2" class="grp">area ratio</th></tr>\n'
        '<tr><th class="num">median</th><th class="num">p99.95</th>'
        '<th class="num">max</th><th class="num">DNC</th>'
        '<th class="num">central&#8209;99.9%</th><th class="num">max/min</th></tr>\n'
        '</thead>\n<tbody>\n' + '\n'.join(body) + '\n</tbody>\n</table>\n')
    out = OUT_DIR / 'summary_table.html'
    out.write_text(html)
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
        counts = (ax.hist(a, bins=bins, weights=d['w_ars'], color=SYS_COLOR[name],
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

    rows = summary_stats(results)
    plot_histograms(results, rows)
    plot_area_histograms(results, rows)
    plot_tradeoff(rows)
    write_summary_table(rows)
    plot_extremes(results)
    for s in SYSTEMS:
        plot_by_resolution(s, results[s]['by_res'])


if __name__ == '__main__':
    main()
