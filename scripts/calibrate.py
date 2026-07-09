"""Match S2 / A5 / DGGAL resolutions to an H3-res-9 cell by area.

Picks, for each DGGS, the resolution whose median cell area is closest (in
log-ratio) to the reference H3 res-9 cell — the resolution baked into
config.TARGET_RES so all systems compare cells of roughly the same size.

Reads the `area` column straight off the tables (computed at generation
time by sparea) — nothing is measured here.

Adding a new DGGS: drop its module in systems/, `just gen`, run this (it
scans every cached system across all resolutions automatically), then bake
the pick into config.TARGET_RES.

Run with:  just calibrate
No CLI args (project convention) — edit the constants below in place.
"""

import numpy as np

from dggs_compare import cache, config

# ----- knobs -------------------------------------------------------------
TARGET = ('h3', config.TARGET_RES['h3'])   # reference system + resolution
PRINT_WINDOW = 3            # rows shown either side of each pick
# -------------------------------------------------------------------------


def median_area(dggs, res):
    """Median cell area (steradians) from the table's area column."""
    return float(np.median(cache.load_columns(dggs, res, ['area'])['area']))


def main():
    tsys, tres = TARGET
    target = median_area(tsys, tres)
    print(f'target: {tsys} r{tres} median area = {target:.4e} sr  '
          f'(N_CELLS={config.N_CELLS:,}, '
          f'seed={config.SEED:#x})\n')

    for s in cache.available_systems():
        if s == tsys:                        # the reference, not a candidate
            continue
        rows = [(res, median_area(s, res))
                for res in cache.available_resolutions(s)]
        best = min(rows, key=lambda r: abs(np.log(r[1] / target)))
        bi = next(i for i, (res, _) in enumerate(rows) if res == best[0])
        print(f'--- {s} (target {tsys} r{tres}) ---')
        print(f'{"res":>4} {"area_sr":>11} {"ratio":>8}')
        for res, area in rows[max(0, bi - PRINT_WINDOW):bi + PRINT_WINDOW + 1]:
            mark = '  <== pick' if res == best[0] else ''
            print(f'{res:>4} {area:>11.4e} {area / target:>8.3f}{mark}')
        print(f'-> {s} r{best[0]}  ({best[1] / target:.3f}x {tsys} r{tres})\n')


if __name__ == '__main__':
    main()
