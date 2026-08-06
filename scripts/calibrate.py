"""Match each system's working resolution to an H3-res-9 cell, by cell count.

Cells partition the sphere, so a resolution's MEAN cell area is exactly
4*pi/N(res) — matching counts matches average cell size, no tables needed
(issue #32): config.CELLS_PER_RES holds closed-form counts for every
system. The pick is argmin_r |log(N_sys(r)/N_h3(9))| — the symmetric
size-ratio; a plain count difference would bias toward coarser picks.

Adding a new DGGS: add its count formula to config.CELLS_PER_RES, run
this, bake the pick into config.TARGET_RES — before any tables exist.

Run with:  just calibrate
No CLI args (project convention) — edit the constants below in place.
"""

import math

from dggs_compare import config

# ----- knobs -------------------------------------------------------------
TARGET = ('h3', config.TARGET_RES['h3'])   # reference system + resolution
PRINT_WINDOW = 2            # rows shown either side of each pick
# -------------------------------------------------------------------------


def main():
    tsys, tres = TARGET
    anchor = config.CELLS_PER_RES[tsys](tres)
    print(f'target: {tsys} r{tres} = {anchor:,} cells '
          f'(mean cell area {4 * math.pi / anchor:.4e} sr)\n')

    for s, n_of in config.CELLS_PER_RES.items():
        if s == tsys:                        # the reference, not a candidate
            continue

        def mismatch(r):
            return abs(math.log(n_of(r) / anchor))

        best = 0
        while mismatch(best + 1) < mismatch(best):   # unimodal: walk down
            best += 1

        print(f'--- {s} (target {tsys} r{tres}) ---')
        print(f'{"res":>4} {"cells":>18} {"area_ratio":>10}')
        for r in range(max(0, best - PRINT_WINDOW), best + PRINT_WINDOW + 1):
            mark = '  <== pick' if r == best else ''
            print(f'{r:>4} {n_of(r):>18,} {anchor / n_of(r):>10.3f}{mark}')
        baked = config.TARGET_RES.get(s)
        note = '' if baked == best else f'  (config.TARGET_RES has r{baked}!)'
        print(f'-> {s} r{best}  ({anchor / n_of(best):.3f}x {tsys} r{tres} '
              f'cell area){note}\n')


if __name__ == '__main__':
    main()
