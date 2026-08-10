"""Cross-implementation agreement report for one grid (issue #57).

For every pair of implementations with raw geometry in data/raw/, compare
their FULL-ENUMERATION resolutions (rows < the file's stamped n_cells
budget — the guaranteed-complete ones): cell counts, nearest-centroid
cell matching, vertex agreement (Hausdorff over matched cells), per-cell
AR and area, and the best rigid rotation between the constructions
(Kabsch over matched vertices). A pure convention difference in grid
placement shows up as a clean rotation whose removal collapses the
residual to numerics; anything that survives the fit is a real geometric
disagreement.

Sampled resolutions are deliberately excluded: implementations may
interpret the same sample point on different surfaces (geodetic vs
spherical), so their sampled cell sets differ by up to the authalic
shift without any grid disagreement.

Manual tool — generate local raw for each implementation first (the tiny
budget suffices; full enumerations are identical at any budget):

    DGGS_COMPARE_N_CELLS=1000 just gen isea3h-dggal
    DGGS_COMPARE_N_CELLS=1000 just gen isea3h-dggrid
    just cross-impl isea3h
"""

import os
from itertools import combinations

import numpy as np
import pyarrow.parquet as pq

from dggs_compare import cache, runner, stats
from dggs_compare.cache import open_ring

GRID = os.environ.get('DGGS_COMPARE_GRID', '')
MAX_N = 20_000    # skip full enumerations larger than this (memory/csar)
KM = 111.32       # degrees -> km, for readability of the residuals


def unit(latlng):
    la, lo = np.radians(np.asarray(latlng, float)).T
    return np.column_stack([np.cos(la) * np.cos(lo),
                            np.cos(la) * np.sin(lo), np.sin(la)])


def full_resolutions(key):
    """{res: rows} for `key`'s guaranteed-complete raw files: fewer rows
    than the stamped budget means enumerate-all (a file exactly AT the
    budget is ambiguous and skipped, with MAX_N the practical bound
    anyway)."""
    out = {}
    for p in runner.RAW_DIR.glob('*.parquet'):
        parsed = cache.parse_table_name(p.name)
        if not parsed or cache.key(parsed[0], parsed[1]) != key:
            continue
        f = pq.ParquetFile(p)
        n_cells = int(f.metadata.metadata[b'n_cells'].decode())
        if f.metadata.num_rows < min(n_cells, MAX_N + 1):
            out[parsed[2]] = f.metadata.num_rows
    return out


def load(key, res):
    t = pq.read_table(runner.raw_path(key, res), columns=['verts'])
    return [np.asarray(open_ring(v), float) for v in t['verts'].to_pylist()]


def centroids(rings):
    c = np.array([unit(r).mean(0) for r in rings])
    return c / np.linalg.norm(c, axis=1, keepdims=True)


def nearest(A, B, block=1000):
    """Index of the nearest row of B for each row of A (unit vectors)."""
    out = np.empty(len(A), dtype=int)
    for lo in range(0, len(A), block):
        out[lo:lo + block] = np.argmax(A[lo:lo + block] @ B.T, axis=1)
    return out


def kabsch(A, B):
    """Rotation minimizing ||R A - B|| over matched unit vectors."""
    U, _, Vt = np.linalg.svd(A.T @ B)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return Vt.T @ np.diag([1, 1, d]) @ U.T


def compare(key_a, key_b, res):
    ra, rb = load(key_a, res), load(key_b, res)
    nn = nearest(centroids(ra), centroids(rb))
    haus = 0.0          # worst vertex mismatch over matched cells (deg)
    d_ar = d_area = 0.0
    pairs_a, pairs_b = [], []
    for i, j in enumerate(nn):
        a, b = unit(ra[i]), unit(rb[j])
        d = np.arccos(np.clip(a @ b.T, -1, 1))
        haus = max(haus, max(d.min(1).max(), d.min(0).max()))
        ar_a, area_a = stats.cell_stats(ra[i])
        ar_b, area_b = stats.cell_stats(rb[j])
        if not (np.isnan(ar_a) or np.isnan(ar_b)):
            d_ar = max(d_ar, abs(ar_a - ar_b))
        d_area = max(d_area, abs(area_a - area_b) / area_b)
        # Kabsch wants 1:1 vertex pairs: same count, unambiguous match.
        k = np.argmax(a @ b.T, axis=1)
        if len(a) == len(b) and len(set(k.tolist())) == len(k):
            pairs_a.append(a)
            pairs_b.append(b[k])
    A, B = np.vstack(pairs_a), np.vstack(pairs_b)
    R = kabsch(A, B)
    angle = np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1)))
    resid = np.degrees(np.arccos(np.clip(((A @ R.T) * B).sum(1), -1, 1)))
    print(f'  r{res:<2} {len(ra):>6} vs {len(rb):>6} cells | vertex max '
          f'{np.degrees(haus) * KM:9.4f} km | max |dAR| {d_ar:.1e} | max '
          f'rel dArea {d_area:.1e} | fit: rot {angle:.6f} deg, residual '
          f'{resid.max() * KM:.4f} km')


def main():
    assert GRID, 'set the grid: just cross-impl <grid>'
    keys = sorted(k for k in {cache.key(p[0], p[1])
                              for p in map(cache.parse_table_name,
                                           (f.name for f in
                                            runner.RAW_DIR.glob('*.parquet')))
                              if p} if cache.parse_key(k)[0] == GRID)
    assert len(keys) >= 2, f'{GRID}: need raw from >=2 implementations, ' \
                           f'found {keys or "none"} (run `just gen` first)'
    for key_a, key_b in combinations(keys, 2):
        shared = sorted(set(full_resolutions(key_a))
                        & set(full_resolutions(key_b)))
        print(f'{key_a} vs {key_b} — full-enumeration resolutions '
              f'{shared or "NONE (nothing guaranteed-complete on both sides)"}')
        for res in shared:
            compare(key_a, key_b, res)


if __name__ == '__main__':
    main()
