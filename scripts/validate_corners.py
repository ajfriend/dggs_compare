"""Validate every system's stats-input rings against refined references.

The tables feed csar a cell's *corner* vertices by default. For geodesic
edges that's exact (the min-enclosing ellipse of the corners already
contains them, by convexity), but several systems have non-geodesic or
even kinked edges, so this checks — for EVERY system in the registry —
whatever ring the pipeline actually feeds csar (corners, or the system's
stats_rings override) against that system's refined_boundaries reference,
across levels including the coarsest and the pentagons. isea3h is
expected red (~4e-3): its odd-level corners-only approximation is a
documented decision (issue #25). If the max delta elsewhere is within
solver tolerance, the stats inputs are confirmed faithful.

Prints a report; writes nothing. Run whenever a new grid is added — this
is the check that admits it to the pipeline.

Run with:  just validate-corners
No CLI args (project convention) — edit the knobs below.
"""

import numpy as np

import csar

from dggs_compare import registry, stats
from dggs_compare.cache import open_ring

# ----- knobs -------------------------------------------------------------
SEED = 0xC0FFEE
REF_REFINE = 40             # reference-boundary refinement; must stay finer
                            # than any system's stats ring
LEVELS = [0, 1, 2, 3, 5, 8, 11]   # incl coarsest + the pentagons (level 0)
K = 300                     # cells tested per level (enumerate if fewer exist)
# -------------------------------------------------------------------------


def ar(ring):
    r = csar.solve(csar.to_vec3(ring, geo='latlng_deg'), geo='vec3')
    return r.aspect_ratio if isinstance(r, csar.Converged) else None


def cells_for_level(mod, level, rng):
    total = mod.num_cells(level)
    if total <= K:
        return list(mod.enumerate_cells(level))
    if total <= 4 * K:
        # Coupon-collector territory: enumerate + subsample instead of
        # point-sampling (mirrors cache._select_zones's middle regime).
        zones = list(mod.enumerate_cells(level))
        idx = rng.choice(len(zones), K, replace=False)
        return [zones[i] for i in idx]
    seen, out = set(), []
    drawn = 0
    while len(out) < K:
        if drawn >= 60 * K:   # cap, like cache.MAX_DRAW_FACTOR
            raise RuntimeError(
                f'{drawn:,} draws yielded only {len(out)}/{K} distinct cells')
        pts = stats.sample_uniform_latlng(K * 4, rng).tolist()
        drawn += len(pts)
        for z in mod.cells_at(level, pts):
            if z is not None and z not in seen:
                seen.add(z)
                out.append(z)
                if len(out) >= K:
                    break
    return out


def check(name, mod):
    """Report stats-ring-vs-reference max |dAR| per level; return the max."""
    rng = np.random.default_rng(SEED)
    stats_rings = getattr(mod, 'stats_rings', None)
    src = 'corners' if stats_rings is None else 'stats_rings'
    print(f'\n{name} {src}-vs-refined (refine={REF_REFINE})')
    print(f'{"lvl":>3} {"cells":>6} {"pents":>6} {"max|dAR|":>10} '
          f'{"max_rel":>10} {"stats_AR_range":>22}')
    overall = 0.0
    for level in LEVELS:
        if level not in mod.resolutions():
            continue
        zones = cells_for_level(mod, level, rng)
        corners = mod.boundaries(level, zones)
        refs = mod.refined_boundaries(level, zones, REF_REFINE)
        override = (stats_rings(level, zones) if stats_rings
                    else [None] * len(zones))
        npent = 0
        max_abs = max_rel = 0.0
        ars = []
        for ring, ref, sring in zip(corners, refs, override):
            if len(ring) == 5:
                npent += 1
            # open_ring, exactly as cache.build_table feeds the solvers —
            # the validator must test the pipeline's true stats input.
            a_c = ar(open_ring(ring if sring is None else sring))
            a_r = ar(open_ring(ref))
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
    worst = {}
    for name in registry.names():
        worst[name] = check(name, registry.get(name))
    print()
    bad = {n: w for n, w in worst.items() if w >= 1e-3 and n != 'isea3h'}
    print(f'overall max |dAR| per system: '
          + '  '.join(f'{n}={w:.1e}' for n, w in sorted(worst.items())))
    if bad:
        print(f'stats-input delta NOT negligible for {sorted(bad)} — investigate')
    else:
        print('stats inputs CONFIRMED (within solver tolerance; isea3h red '
              'by design, see issue #25)')


if __name__ == '__main__':
    main()
