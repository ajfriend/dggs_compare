"""Validate the corners-only enclosing-cone metric for the DGGAL grids.

The tables feed csar only a cell's *corner* vertices. For geodesic edges
that's exact: the min-enclosing ellipse of the corners already contains them
(convexity). The equal-area DGGAL grids have slightly *non-geodesic* edges
that could bow outward at coarse levels, so this checks it empirically for
every DGGAL-backed system in the registry: across levels — including the
coarsest and the 12 pentagons (level 0) — it compares the aspect ratio from
corners against the ratio from edge-refined vertices. If the max delta is
within solver tolerance, corners-only is confirmed.

Prints a report; writes nothing. Run whenever a new DGGAL grid is added —
this is the check that admits it to the pipeline.

Run with:  just validate-corners
No CLI args (project convention) — edit the knobs below.
"""

import numpy as np

import csar

from dggs_compare import registry

# ----- knobs -------------------------------------------------------------
SEED = 0xC0FFEE
REFINE = 20                 # edge-refinement points per edge for the reference
LEVELS = [0, 1, 2, 3, 5, 8, 11]   # incl coarsest + the 12 pentagons (level 0)
K = 300                     # cells tested per level (enumerate if fewer exist)
# -------------------------------------------------------------------------


def ar(verts):
    r = csar.solve(verts, geo='vec3')
    return r.aspect_ratio if isinstance(r, csar.Converged) else None


def cells_for_level(ad, level, rng):
    if ad.count(level) <= K:
        return list(ad.enumerate(level))
    seen, out = set(), []
    for zone in ad.sample(level, K * 4, rng):
        if zone not in seen:
            seen.add(zone)
            out.append(zone)
            if len(out) >= K:
                break
    return out


def check(name, ad):
    """Report corners-vs-refined max |dAR| per level; return the overall max."""
    rng = np.random.default_rng(SEED)
    print(f'\n{name} corners-vs-refined (edgeRefinement={REFINE})')
    print(f'{"lvl":>3} {"cells":>6} {"pents":>6} {"max|dAR|":>10} '
          f'{"max_rel":>10} {"corners_AR_range":>22}')
    overall = 0.0
    for level in LEVELS:
        if level > ad.max_level():
            continue
        zones = cells_for_level(ad, level, rng)
        npent = 0
        max_abs = max_rel = 0.0
        ars = []
        for z in zones:
            corners = ad.verts(z)
            if corners.shape[0] == 5:
                npent += 1
            a_c = ar(corners)
            a_r = ar(csar.to_vec3(ad.refined_boundary(z, REFINE),
                                  geo='latlng_deg'))
            if a_c is None or a_r is None:
                continue
            ars.append(a_c)
            d = abs(a_c - a_r)
            max_abs = max(max_abs, d)
            max_rel = max(max_rel, d / a_r)
        overall = max(overall, max_abs)
        rng_txt = f'[{min(ars):.4f}, {max(ars):.4f}]' if ars else 'n/a'
        print(f'{level:>3} {len(zones):>6} {npent:>6} {max_abs:>10.2e} '
              f'{max_rel:>10.2e} {rng_txt:>22}')
    return overall


def main():
    worst = 0.0
    for name in registry.names():
        sysmod = registry.get(name)
        if not hasattr(sysmod, 'adapter'):     # DGGAL-backed systems only
            continue
        worst = max(worst, check(name, sysmod.adapter()))
    print(f'\noverall max |dAR| across all grids = {worst:.3e}')
    print('corners-only CONFIRMED (within solver tolerance)' if worst < 1e-3
          else 'corners-only delta NOT negligible — investigate')


if __name__ == '__main__':
    main()
