"""Data-quality checks: the DNC invariants and the corners-only validation.

Both are library functions returning structured results; the thin
scripts/dnc_check.py and scripts/validate_corners.py print/plot and set exit
codes.

DNC sweep modes (`resolve=`):
  False (default) — read the cached `ar` column; DNC = NaN. Fast: checks the
      published artifact itself.
  True — re-solve every cell's `verts` with the *installed* skar at the
      config solver settings, ignoring the cached stats. This is the skar
      pre-release regression gate: point the pyproject skar pin at a release
      candidate, `uv sync`, and run the gate — no table rebuild needed.
"""

import numpy as np

from . import cache, config

# DNC-fraction noise floor — the finest resolutions are N_SMALL=25k cells
# (~0.3% sampling noise), so a stray cell or two at the f64 floor isn't a
# real band.
NOISE_TOL = 1e-2
MAX_EXAMPLES = 5      # offending cell ids reported per failing resolution


def sweep_system(name, *, resolve=False):
    """[(res, tested, dnc, [example cids])] over the system's tables."""
    rows = []
    for res in cache.available_resolutions(name):
        if resolve:
            import skar
            tested = dnc = 0
            examples = []
            for cid, latlng in cache.load_cells(name, res):
                tested += 1
                r = skar.solve(skar.to_vec3(latlng, geo='latlng_deg'),
                               geo='vec3', gap_tol=config.GAP_TOL,
                               method=config.SKAR_METHOD)
                if not isinstance(r, skar.Converged):
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
    frac = {res: dnc / tested for res, tested, dnc, _ in rows}
    reslist = [res for res, *_ in rows]
    failures = []

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
    finest_res, finest_tested, finest_dnc, _ = rows[-1]
    return failures, onset, finest_dnc / finest_tested
