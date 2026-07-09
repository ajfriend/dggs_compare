"""DGGS aspect-ratio survey — reads the `ar` column of the tables (no
solving) and writes the comparison plots to out/:

  histograms.png     cross-system AR distributions at the working resolutions
  extremes.png       best/worst cell per system, drawn with its ellipse
  by_res_<sys>.png   per-system AR distribution stacked by resolution

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

import skar

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
        verts = skar.to_vec3(cols['verts'][idx], geo='latlng_deg')
        r = skar.solve(verts, geo='vec3', gap_tol=config.GAP_TOL,
                       method=config.SKAR_METHOD)
        return {'ar': float(cols['ar'][idx]), 'id': cols['cid'][idx],
                'verts': verts, 'result': r}

    return {'by_res': by_res,
            'best': record(int(np.nanargmin(cols['ar']))),
            'worst': record(int(np.nanargmax(cols['ar'])))}


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
    xy, semi = skar.project_to_cone(rec['result'], rec['verts'], up=None)
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


if __name__ == '__main__':
    main()
