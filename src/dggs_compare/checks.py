"""Data-quality checks: the DNC invariants and artifact/config coherence.

All checks are pure table/metadata/filename reads — no DGGS binding is
ever imported. The implementation registry is the scripts/systems/ file
listing (one script per (grid, impl), named '{grid}-{impl}.py'); each
table's declared resolution coverage comes from the 'resolutions' key the
runner stamps into its metadata.

Library functions returning structured results; the thin
scripts/dnc_check.py prints/plots and sets exit codes.

DNC sweep modes (`resolve=`):
  False (default) — read the cached `ar` column; DNC = NaN. Fast: checks the
      published artifact itself.
  True — re-solve every cell's `verts` with the *installed* csar at the
      config solver settings, ignoring the cached stats. This is the csar
      pre-release regression gate: point a [tool.uv.sources] entry at a
      csar_py branch/rev, `uv sync`, and run the gate — no table rebuild
      needed.
"""

import re
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

from . import cache, config

SCRIPTS_DIR = Path(__file__).resolve().parents[2] / 'scripts' / 'systems'

# DNC-fraction noise floor — sampled resolutions are N_CELLS cells
# (sampling noise), so a stray cell or two at the f64 floor isn’t a
# real band.
NOISE_TOL = 1e-2
MAX_EXAMPLES = 5      # offending cell ids reported per failing resolution


def implementations():
    """Sorted (grid, impl) pairs from the scripts/systems/ listing — the
    registry. Scripts are named '{grid}-{impl}.py'."""
    pat = re.compile(r'^([^-]+)-([^-]+)\.py$')
    return sorted((m.group(1), m.group(2))
                  for p in SCRIPTS_DIR.glob('*.py')
                  if (m := pat.match(p.name)))


def _declared_resolutions(grid, impl):
    """The 'resolutions' metadata from one of the pair's tables (they all
    carry the same value), or None if no table exists."""
    res_list = cache.available_tables().get((grid, impl))
    if not res_list:
        return None
    path = cache.DATA_DIR / f'{grid}-{impl}_r{res_list[0]}.parquet'
    meta = pq.ParquetFile(path).metadata.metadata or {}
    raw = meta.get(b'resolutions')
    return None if raw is None else {int(r) for r in raw.decode().split(',')}


def missing_systems():
    """Registry implementations absent from the data or the per-grid config.

    Returns (no_tables, no_config): (grid, impl) pairs with no tables in
    data/cells/, and 'grid (DICT)' entries for each config.PER_SYSTEM dict
    a grid is missing from. Both empty = the artifact covers the whole
    registry — the release-gate completeness check.
    """
    have = cache.available_tables()
    impls = implementations()
    no_tables = [gi for gi in impls if gi not in have]
    grids = sorted({g for g, _ in impls})
    no_config = [f'{g} ({k})' for g in grids
                 for k, d in config.PER_SYSTEM.items() if g not in d]
    return no_tables, no_config


def stale_tables():
    """Tables on disk OUTSIDE their own declared resolution coverage.

    Disk == declaration is an invariant every consumer relies on without
    checking — a stale table (left behind when an implementation's max
    resolution was lowered, or a partial write from a crashed run) would
    silently pollute them all, so the gate fails loudly on any. Files not
    matching the {grid}-{impl} naming scheme are ignored (not part of the
    artifact)."""
    out = []
    for (grid, impl), res_list in sorted(cache.available_tables().items()):
        declared = _declared_resolutions(grid, impl)
        out += [f'{grid}-{impl}_r{res}.parquet' for res in res_list
                if declared is not None and res not in declared]
    return out


def target_res_problems():
    """TARGET_RES entries that are drifted or untabled.

    TARGET_RES has a defined correct answer (the count match, issue #32)
    and consumers that trust it blindly (check_system's clean-where-used
    bound, the site manifest, survey) — so the gate asserts it. Also
    requires a primary-implementation table at the target: a target finer
    than the finest table would make the DNC invariant pass vacuously."""
    anchor = config.CELLS_PER_RES['h3'](config.TARGET_RES['h3'])
    problems = []
    for g, baked in config.TARGET_RES.items():
        if g != 'h3' and baked != (pick := config.count_match_res(g, anchor)):
            problems.append(f'{g}: TARGET_RES r{baked} != count-match r{pick}')
        if baked not in cache.available_resolutions(g):
            problems.append(f'{g}: no table at TARGET_RES r{baked}')
    return problems


def sweep_system(name, *, resolve=False):
    """[(res, tested, dnc, [example cids])] over a grid's primary tables."""
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
    """Return (failures, onset_res, finest_frac) for one grid's sweep rows.

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
