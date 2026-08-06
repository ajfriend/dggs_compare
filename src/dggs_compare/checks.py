"""Data-quality checks: the DNC invariants and the corners-only validation.

Both are library functions returning structured results; the thin
scripts/dnc_check.py and scripts/validate_corners.py print/plot and set exit
codes.

DNC sweep modes (`resolve=`):
  False (default) — read the cached `ar` column; DNC = NaN. Fast: checks the
      published artifact itself.
  True — re-solve every cell's `verts` with the *installed* csar at the
      config solver settings, ignoring the cached stats. This is the csar
      pre-release regression gate: point the pyproject csar pin at a release
      candidate, `uv sync`, and run the gate — no table rebuild needed.
"""

import numpy as np

from . import cache, config, registry

# DNC-fraction noise floor — sampled resolutions are N_CELLS cells
# (sampling noise), so a stray cell or two at the f64 floor isn’t a
# real band.
NOISE_TOL = 1e-2
MAX_EXAMPLES = 5      # offending cell ids reported per failing resolution


def missing_systems():
    """Registry systems absent from the data or the per-system config.

    Returns (no_tables, no_config): names with no tables in data/cells/,
    and 'name (DICT)' entries for each config.PER_SYSTEM dict a name is
    missing from. Both empty = the artifact covers the whole registry —
    the release-gate completeness check.
    """
    names = registry.names()
    have = set(cache.available_systems())
    no_tables = [s for s in names if s not in have]
    no_config = [f'{s} ({k})' for s in names
                 for k, d in config.PER_SYSTEM.items() if s not in d]
    return no_tables, no_config


def stale_tables():
    """Tables on disk OUTSIDE their system's declared resolutions().

    Disk == contract is an invariant every table consumer (survey,
    calibrate, webdata, the sweeps below) relies on without checking — a
    stale table (e.g. left behind when a system's MAX_RES was lowered, or
    a partial write from a crashed run) silently pollutes all of them, so
    the gate fails loudly on any rather than one consumer filtering.
    Costs a lazy module import per system, so gate-only.
    """
    return [f'{name}_r{res}.parquet'
            for name in cache.available_systems()
            for declared in [set(registry.get(name).resolutions())]
            for res in cache.available_resolutions(name)
            if res not in declared]


def target_res_problems():
    """TARGET_RES entries that are drifted, out of contract, or untabled.

    TARGET_RES has a defined correct answer (the count match, issue #32) and
    three consumers that trust it blindly (check_system's clean-where-used
    bound, the site manifest, survey) — so the gate asserts it rather than
    relying on someone reading calibrate's output. Also requires the entry
    to be a declared, on-disk resolution: a target finer than the finest
    table would make the DNC invariant pass vacuously.
    """
    anchor = config.CELLS_PER_RES['h3'](config.TARGET_RES['h3'])
    problems = []
    for s, baked in config.TARGET_RES.items():
        if s != 'h3' and baked != (pick := config.count_match_res(s, anchor)):
            problems.append(f'{s}: TARGET_RES r{baked} != count-match r{pick}')
        if baked not in registry.get(s).resolutions():
            problems.append(f'{s}: TARGET_RES r{baked} outside declared '
                            f'resolutions')
        elif baked not in cache.available_resolutions(s):
            problems.append(f'{s}: no table at TARGET_RES r{baked}')
    return problems


def sweep_system(name, *, resolve=False):
    """[(res, tested, dnc, [example cids])] over the system's tables."""
    rows = []
    for res in cache.available_resolutions(name):
        if resolve:
            import csar
            tested = dnc = 0
            examples = []
            for cid, latlng in cache.load_cells(name, res):
                tested += 1
                r = csar.solve(csar.to_vec3(latlng, geo='latlng_deg'),
                               geo='vec3', gap_tol=config.GAP_TOL,
                               method=config.CSAR_METHOD)
                if not isinstance(r, csar.Converged):
                    dnc += 1
                    if len(examples) < MAX_EXAMPLES:
                        examples.append(cid)
        else:
            cols = cache.load_columns(name, res, ['cid', 'ar'])
            bad = np.isnan(cols['ar'])
            tested = len(cols['ar'])
            dnc = int(bad.sum())
            examples = [c for c, b in zip(cols['cid'], bad) if b][:MAX_EXAMPLES]
        rows.append((res, tested, dnc, examples))
    return rows


def check_system(name, rows):
    """Return (failures, onset_res, finest_frac) for one system's sweep rows.

    Invariants:
      1. clean where it's used — 0 DNC at the working (target) resolution and
         all coarser;
      2. monotone — DNC only grows toward the finest resolutions: no
         meaningful DNC band with a clean finer resolution (no "islands"),
         and the fraction never meaningfully drops as resolution rises.
    """
    target = config.TARGET_RES[name]
    # An empty table (a truncated/failed write) is its own failure; count it
    # as all-DNC so the monotonicity logic needs no special cases.
    frac = {res: dnc / tested if tested else 1.0
            for res, tested, dnc, _ in rows}
    reslist = [res for res, *_ in rows]
    failures = [f'r{res}: empty table' for res, tested, _, _ in rows
                if not tested]

    for i, (res, tested, dnc, ex) in enumerate(rows):
        if res <= target and dnc:
            failures.append(f'r{res}: {dnc}/{tested} DNC at a working '
                            f'resolution (<= target r{target}); e.g. {ex}')
        if frac[res] >= NOISE_TOL and any(frac[r] == 0 for r in reslist[i + 1:]):
            clean = [r for r in reslist[i + 1:] if frac[r] == 0]
            failures.append(f'r{res}: {100*frac[res]:.1f}% DNC but finer res '
                            f'{clean} clean (non-monotone island)')
        if i + 1 < len(rows):
            nxt = reslist[i + 1]
            if frac[nxt] + NOISE_TOL < frac[res]:
                failures.append(f'r{nxt}: {100*frac[nxt]:.1f}% DNC < r{res} '
                                f'{100*frac[res]:.1f}% (non-monotone drop)')

    onset = next((res for res, _, dnc, _ in rows if dnc), None)
    return failures, onset, frac[reslist[-1]]
