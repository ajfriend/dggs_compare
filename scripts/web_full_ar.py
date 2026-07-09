"""Solve aspect ratios for the full-enumeration geometry (pass two of the
full-globe pipeline; run web_full_geom.py first).

Levels are discovered from the *_pos.f32 files. A level whose _ar.f32
already exists is skipped — delete it to re-solve (e.g. after a solver
bump). Writes {SYSTEM}_r{res}_ar.f32 (Float32 per cell, NaN = DNC, order
matching _idx.u32).

Run with:  just web-full
No CLI args (project convention) — edit SYSTEM below.
"""

import re
import time
from pathlib import Path

import numpy as np

import skar

from dggs_compare import config

# ----- knobs -------------------------------------------------------------
SYSTEM = 'ivea7h'
OUT = Path(__file__).resolve().parent.parent / 'web' / 'out' / 'full'
# -------------------------------------------------------------------------


def main():
    pos_paths = sorted(OUT.glob(f'{SYSTEM}_r*_pos.f32'),
                       key=lambda p: int(re.search(r'_r(\d+)_', p.name)[1]))
    if not pos_paths:
        print(f'no {SYSTEM} geometry in {OUT} — run `just web-full-geom` first.')
        return
    for pos_path in pos_paths:
        res = int(re.search(r'_r(\d+)_', pos_path.name)[1])
        ar_path = OUT / f'{SYSTEM}_r{res}_ar.f32'
        if ar_path.exists():
            print(f'r{res}: {ar_path.name} exists — skipping (delete to re-solve).')
            continue
        pos = np.fromfile(pos_path, dtype='<f4').reshape(-1, 2)   # [lng, lat]
        idx = np.fromfile(OUT / f'{SYSTEM}_r{res}_idx.u32', dtype='<u4')
        n = len(idx) - 1
        print(f'{SYSTEM} r{res}: solving {n:,} cells...', flush=True)
        t0 = time.perf_counter()
        # One vectorized deg->vec3 pass over all vertices; per-cell slices are
        # contiguous views, so the hot loop allocates nothing.
        xyz = skar.to_vec3(pos[:, ::-1].astype(float), geo='latlng_deg')
        ar = np.empty(n, dtype='<f4')
        dnc = 0
        for i in range(n):
            r = skar.solve(xyz[idx[i]:idx[i + 1]], geo='vec3',
                           gap_tol=config.GAP_TOL, method=config.SKAR_METHOD)
            if isinstance(r, skar.Converged):
                ar[i] = r.aspect_ratio
            else:
                ar[i] = np.nan
                dnc += 1
            if (i + 1) % 200_000 == 0:
                print(f'  {i + 1:,}/{n:,}', flush=True)
        ar.tofile(ar_path)
        finite = ar[np.isfinite(ar)]
        print(f'  wrote {ar_path.name}: DNC {dnc:,}, '
              f'AR min {finite.min():.4f} / median {np.median(finite):.4f} / '
              f'max {finite.max():.4f}, {time.perf_counter() - t0:.1f}s', flush=True)


if __name__ == '__main__':
    main()
