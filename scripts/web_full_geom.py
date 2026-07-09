"""Enumerate EVERY cell of SYSTEM at each RES_LIST level and stream its
geometry to binary — the complete-coverage companion to the sampled tables
(up to 168,072 cells at r5 and 1,176,492 at r6 for ivea7h).

Per level it writes two little-endian binaries to web/out/full/ (gitignored):
  {SYSTEM}_r{res}_pos.f32   Float32 [lng, lat] vertex pairs, flattened, open
                            rings — ajglobe converts each vertex to a unit
                            vector, so no antimeridian/winding preprocessing.
  {SYSTEM}_r{res}_idx.u32   Uint32 start indices (len = n_cells + 1).
Counts are implicit: n_verts = pos.size/2, n_cells = idx.size - 1. Levels
whose binaries already exist are skipped (geometry is deterministic).

Aspect ratios are solved afterwards by web_full_ar.py.

Run with:  just web-full-geom
No CLI args (project convention) — edit SYSTEM/RES_LIST below.
"""

import time
from pathlib import Path

import numpy as np

from dggs_compare import registry

# ----- knobs -------------------------------------------------------------
SYSTEM = 'ivea7h'
RES_LIST = [1, 2, 3, 5, 6]
OUT = Path(__file__).resolve().parent.parent / 'web' / 'out' / 'full'
# -------------------------------------------------------------------------


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sysmod = registry.get(SYSTEM)
    for res in RES_LIST:
        pos_path = OUT / f'{SYSTEM}_r{res}_pos.f32'
        idx_path = OUT / f'{SYSTEM}_r{res}_idx.u32'
        if pos_path.exists() and idx_path.exists():
            print(f'r{res}: {pos_path.name} exists — skipping (geometry is '
                  f'deterministic; delete to regenerate).', flush=True)
            continue
        n_cells = sysmod.num_cells(res)
        print(f'{SYSTEM} r{res}: {n_cells:,} cells — enumerating...', flush=True)
        t0 = time.perf_counter()
        starts = [0]
        nverts = 0
        done = 0
        with open(pos_path, 'wb') as fpos:
            for z in sysmod.enumerate_cells(res):
                ring = sysmod.cell_boundary(z)
                # [lat, lng] rows -> [lng, lat] (the browser axis order)
                arr = np.asarray(ring, dtype='<f4')[:, ::-1]
                fpos.write(arr.tobytes())
                nverts += len(ring)
                starts.append(nverts)
                done += 1
                if done % 100_000 == 0:
                    print(f'  {done:,}/{n_cells:,}', flush=True)
        np.asarray(starts, dtype='<u4').tofile(idx_path)
        mb = (pos_path.stat().st_size + idx_path.stat().st_size) / 1e6
        print(f'  wrote {pos_path.name} + {idx_path.name}: {nverts:,} verts, '
              f'{mb:.1f} MB, {time.perf_counter() - t0:.1f}s', flush=True)


if __name__ == '__main__':
    main()
