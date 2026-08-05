"""Validate the stats-input ring against a refined reference, per DGGAL grid.

The tables feed csar a cell's *corner* vertices by default. For geodesic
edges that's exact: the min-enclosing ellipse of the corners already contains
them (convexity). The equal-area DGGAL grids have slightly *non-geodesic*
edges that could bow outward at coarse levels — and isea3h's odd-level cells
genuinely kink at icosahedron edges, which is why that system feeds the
solvers edge-refined `stats_rings` instead. So this checks, for every
DGGAL-backed system in the registry, whatever ring the pipeline ACTUALLY
feeds csar (corners, or the system's stats_rings override) against a
doubly-refined reference boundary, across levels — including the coarsest
and the 12 pentagons (level 0). If the max delta is within solver tolerance,
the stats input is confirmed faithful.

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
REF_REFINE = 40             # reference-boundary refinement; must stay finer
                            # than any system's stats ring (isea3h uses 20)
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


def check(name, ad, stats_rings=None):
    """Report stats-ring-vs-reference max |dAR| per level; return the max."""
    rng = np.random.default_rng(SEED)
    src = 'corners' if stats_rings is None else 'stats_rings'
    print(f'\n{name} {src}-vs-refined (edgeRefinement={REF_REFINE})')
    print(f'{"lvl":>3} {"cells":>6} {"pents":>6} {"max|dAR|":>10} '
          f'{"max_rel":>10} {"stats_AR_range":>22}')
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
            ring = stats_rings([z])[0] if stats_rings else None
            stats_in = corners if ring is None else csar.to_vec3(
                ring, geo='latlng_deg')
            a_c = ar(stats_in)
            a_r = ar(csar.to_vec3(ad.refined_boundary(z, REF_REFINE),
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
        worst = max(worst, check(name, sysmod.adapter(),
                                 getattr(sysmod, 'stats_rings', None)))
    print(f'\noverall max |dAR| across all grids = {worst:.3e}')
    print('stats input CONFIRMED (within solver tolerance)' if worst < 1e-3
          else 'stats-input delta NOT negligible — investigate')


if __name__ == '__main__':
    main()
