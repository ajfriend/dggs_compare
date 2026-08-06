"""Match each system's working resolution to an H3-res-9 cell, by cell count.

Matching counts matches mean cell area with no tables read (see the
config.CELLS_PER_RES comment for the derivation and config.count_match_res
for the criterion — issue #32). Prints each pick next to the baked
config.TARGET_RES and flags disagreements; the dnc-check gate asserts the
same agreement, so this script is the human-readable view, not the
enforcement.

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
        best = config.count_match_res(s, anchor)
        print(f'--- {s} ---')
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
