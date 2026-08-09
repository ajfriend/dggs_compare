"""Verify every system's metrics are CONVERGED in edge-sampling density.

A cell's boundary is a list of vertices sampled along its edges; the
sampling density is the only knob (see the registry contract). The tables
feed the solvers density-0 vertex lists. For great-circle-edge systems
that's exact (the min-enclosing ellipse of the vertices already contains
the edges, by convexity); for systems with other native edge models this
script measures the residual: for EVERY system in the registry, compare
the metric at density 0 against a high-density reference, across levels
including the coarsest and the pentagons. isea3h is expected red (~4e-3):
its odd-level density-0 approximation is a documented decision (issue
#25). If the max delta elsewhere is within solver tolerance, the stats
inputs are confirmed faithful.

Prints a report; writes nothing. Run whenever a new grid is added — this
is the check that admits it to the pipeline.

Run with:  just convergence
No CLI args (project convention) — edit the knobs below.
"""

import numpy as np

import csar

from dggs_compare import registry, stats
from dggs_compare.cache import open_ring

# ----- knobs (shared with the two-stage pipeline's runner) ---------------
from dggs_compare import config  # noqa: E402

SEED = config.CONV_SEED
REF_SAMPLES = config.CONV_SAMPLES   # must stay denser than the tables use
LEVELS = config.CONV_LEVELS         # incl coarsest + the pentagons
K = config.CONV_K
# -------------------------------------------------------------------------


def ar(verts):
    r = csar.solve(csar.to_vec3(verts, geo='latlng_deg'), geo='vec3')
    return r.aspect_ratio if isinstance(r, csar.Converged) else None


def cells_for_level(mod, level, rng):
    total = mod.num_cells(level)
    if total <= K:
        return list(mod.enumerate_cells(level))
    if total <= 4 * K:
        # Coupon-collector territory: enumerate + subsample instead of
        # point-sampling (mirrors cache._select_cells's middle regime).
        cells = list(mod.enumerate_cells(level))
        idx = rng.choice(len(cells), K, replace=False)
        return [cells[i] for i in idx]
    seen, out = set(), []
    drawn = 0
    while len(out) < K:
        if drawn >= 60 * K:   # cap, like cache.MAX_DRAW_FACTOR
            raise RuntimeError(
                f'{drawn:,} draws yielded only {len(out)}/{K} distinct cells')
        pts = stats.sample_uniform_latlng(K * 4, rng).tolist()
        drawn += len(pts)
        for c in mod.cells_at(level, pts):
            if c is not None and c not in seen:
                seen.add(c)
                out.append(c)
                if len(out) >= K:
                    break
    return out


def check(name, mod):
    """Report density-0-vs-reference max |dAR| per level; return the max."""
    rng = np.random.default_rng(SEED)
    print(f'\n{name} density-0 vs density-{REF_SAMPLES}')
    print(f'{"lvl":>3} {"cells":>6} {"pents":>6} {"max|dAR|":>10} '
          f'{"max_rel":>10} {"density0_AR_range":>22}')
    overall = 0.0
    for level in LEVELS:
        if level not in mod.resolutions():
            continue
        cells = cells_for_level(mod, level, rng)
        base = mod.boundaries(level, cells)
        refs = mod.boundaries(level, cells, REF_SAMPLES)
        npent = 0
        max_abs = max_rel = 0.0
        ars = []
        for verts, ref in zip(base, refs):
            if len(verts) == 5:
                npent += 1
            # open_ring, exactly as cache.build_table feeds the solvers —
            # this must test the pipeline's true stats input.
            a_0 = ar(open_ring(verts))
            a_r = ar(open_ring(ref))
            if a_0 is None or a_r is None:
                continue
            ars.append(a_0)
            d = abs(a_0 - a_r)
            max_abs = max(max_abs, d)
            max_rel = max(max_rel, d / a_r)
        overall = max(overall, max_abs)
        rng_txt = f'[{min(ars):.4f}, {max(ars):.4f}]' if ars else 'n/a'
        print(f'{level:>3} {len(cells):>6} {npent:>6} {max_abs:>10.2e} '
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
        print(f'density-0 delta NOT negligible for {sorted(bad)} — investigate')
    else:
        print('stats inputs CONFIRMED converged (within solver tolerance; '
              'isea3h red by design, see issue #25)')


if __name__ == '__main__':
    main()
