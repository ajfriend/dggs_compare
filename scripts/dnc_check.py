"""DNC invariant check across every DGGS and resolution — pass/fail.

Reads the cached `ar` column (or re-solves, see RESOLVE) and verifies:
  1. clean where it's used — 0 DNC at each system's working resolution and
     all coarser;
  2. monotone — DNC only grows toward the finest, sub-metre resolutions
     (no islands, no drops).

Exits non-zero (printing the offenders) if either invariant breaks. Writes
one diagnostic plot (out/dnc_check.png: DNC % vs resolution per system).

RESOLVE=True re-solves every cell with the *installed* csar instead of
reading the cached stats — the csar pre-release regression gate (point the
pyproject csar pin at the candidate, `uv sync`, run this).

Run with:  just dnc-check
No CLI args (project convention) — edit RESOLVE below.
"""

import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

from dggs_compare import cache, checks, config

# ----- knobs -------------------------------------------------------------
RESOLVE = False           # True: re-solve with the installed csar (release gate)
OUT = Path(__file__).resolve().parent.parent / 'out' / 'dnc_check.png'
# -------------------------------------------------------------------------


def plot(per_system):
    """DNC % vs resolution, a line per system: flat-zero until the finest
    resolutions then a monotone rise — both invariants at a glance."""
    fig, ax = plt.subplots(figsize=(9, 5))
    for key, rows in per_system.items():
        grid, _ = cache.parse_key(key)
        res = [r for r, *_ in rows]
        pct = [100 * dnc / tested if tested else 100.0
               for _, tested, dnc, _ in rows]
        ax.plot(res, pct, '-o', ms=4, color=config.SYS_COLOR.get(grid), label=key)
    ax.set_xlabel('resolution / level')
    ax.set_ylabel('did-not-converge (%)')
    ax.set_title(f'DGGS DNC fraction vs resolution (gap_tol = {config.GAP_TOL:g})')
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=150)
    plt.close(fig)
    print(f'wrote {OUT}')


def main():
    # The release gate must see EVERY grid: a registry system without tables
    # or per-system config is a broken artifact, not a smaller one. (On a
    # local checkout this also means: fetch a release that matches the code.)
    no_tables, no_config = checks.missing_systems()
    if no_tables:
        print(f'FAIL — no tables in data/cells/: {no_tables} '
              '(fetch an up-to-date release?)')
    if no_config:
        print(f'FAIL — missing config.PER_SYSTEM entries: {no_config}')
    if no_tables or no_config:
        sys.exit(1)
    coverage = checks.coverage_problems()
    if coverage:
        print('FAIL — disk/declaration coverage problems (a releasable '
              'artifact has none; a partial local build reports its holes '
              'here):')
        for p in coverage:
            print(f'  {p}')
        sys.exit(1)
    bad_targets = checks.target_res_problems()
    if bad_targets:
        print(f'FAIL — config.TARGET_RES problems: {bad_targets} '
              '(run `just calibrate`)')
        sys.exit(1)

    print(f'{"implementation":16} {"onset":>6} {"finest %DNC":>12} {"result":>8}')
    all_failures = []
    per_system = {}
    for grid, impl in checks.implementations():
        key = cache.key(grid, impl)
        t0 = time.perf_counter()
        rows = checks.sweep_system(grid, impl=impl, resolve=RESOLVE)
        per_system[key] = rows
        failures, onset, finest_frac = checks.check_system(grid, rows)
        all_failures += [(key, f) for f in failures]
        print(f'{key:16} {("r%d" % onset if onset is not None else "none"):>6} '
              f'{100*finest_frac:>11.1f}% {"FAIL" if failures else "OK":>8} '
              f'  ({time.perf_counter() - t0:.1f}s)')

    plot(per_system)
    if all_failures:
        print('\nFAIL — DNC invariant(s) broken:')
        for name, f in all_failures:
            print(f'  {name:8} {f}')
        sys.exit(1)
    print('\nPASS: working resolutions clean; DNC only at the finest, '
          'sub-metre resolutions, and monotone.')


if __name__ == '__main__':
    main()
