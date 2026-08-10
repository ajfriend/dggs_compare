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

import csar
import numpy as np
import pyarrow.parquet as pq

from dggs_compare import cache, config, metrics, runner, stats
from dggs_compare.cache import open_ring

GRID = os.environ.get('DGGS_COMPARE_GRID', '')
MAX_N = 20_000    # skip full enumerations larger than this (memory/csar)
KM = np.radians(1.0) * config.R_KM     # degrees -> km along a great circle

# Known, explained divergences, annotated on the affected report rows —
# this tool's analog of config.CONV_EXPECTED_RED: (grid, res % 2) -> why.
EXPECTED = {('isea7h', 1): 'density-0 fold-crossing conventions differ '
                           'at odd levels (issue #57)'}


def full_resolutions(key, res_list):
    """The subset of `res_list` whose raw files are guaranteed complete:
    fewer rows than the stamped budget means enumerate-all (a file
    exactly AT the budget is ambiguous and skipped), capped at MAX_N."""
    out = set()
    for res in res_list:
        f = pq.ParquetFile(runner.raw_path(key, res))
        rows = f.metadata.num_rows
        budget = int(f.metadata.metadata[b'n_cells'].decode())
        if rows < budget and rows <= MAX_N:
            out.add(res)
    return out


def load(key, res):
    """([latlng ring, ...], [unit-vec3 ring, ...]) for one raw file."""
    t = pq.read_table(runner.raw_path(key, res), columns=['verts'])
    rings = [np.asarray(open_ring(v), float) for v in t['verts'].to_pylist()]
    return rings, [csar.to_vec3(r, geo='latlng_deg') for r in rings]


def centroids(vecs):
    c = np.array([v.mean(0) for v in vecs])
    return c / np.linalg.norm(c, axis=1, keepdims=True)


def nearest(A, B):
    """Index of the nearest row of B for each row of A (unit vectors),
    in blocks to bound the distance matrix."""
    out = np.empty(len(A), dtype=int)
    for lo in range(0, len(A), 1000):
        out[lo:lo + 1000] = np.argmax(A[lo:lo + 1000] @ B.T, axis=1)
    return out


def kabsch(A, B):
    """Rotation minimizing ||R A - B|| over matched unit vectors."""
    U, _, Vt = np.linalg.svd(A.T @ B)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return Vt.T @ np.diag([1, 1, d]) @ U.T


def compare(key_a, key_b, res):
    ra, va = load(key_a, res)
    rb, vb = load(key_b, res)
    nn = nearest(centroids(va), centroids(vb))
    haus = 0.0          # worst vertex mismatch over matched cells (deg)
    d_ar = d_area = 0.0
    pairs_a, pairs_b = [], []
    for i, j in enumerate(nn):
        a, b = va[i], vb[j]
        d = np.arccos(np.clip(a @ b.T, -1, 1))
        haus = max(haus, max(d.min(1).max(), d.min(0).max()))
        ar_a, area_a = stats.cell_stats(ra[i])
        ar_b, area_b = stats.cell_stats(rb[j])
        if not (np.isnan(ar_a) or np.isnan(ar_b)):
            d_ar = max(d_ar, abs(ar_a - ar_b))
        d_area = max(d_area, abs(area_a - area_b) / area_b)
        # Kabsch wants 1:1 vertex pairs: same count, unambiguous match.
        k = d.argmin(axis=1)
        if len(a) == len(b) and len(set(k.tolist())) == len(k):
            pairs_a.append(a)
            pairs_b.append(b[k])
    A, B = np.vstack(pairs_a), np.vstack(pairs_b)
    R = kabsch(A, B)
    angle = np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1)))
    resid = np.degrees(np.arccos(np.clip(((A @ R.T) * B).sum(1), -1, 1)))
    note = EXPECTED.get((cache.parse_key(key_a)[0], res % 2))
    print(f'  r{res:<2} {len(ra):>6} vs {len(rb):>6} cells | vertex max '
          f'{np.degrees(haus) * KM:9.4f} km | max |dAR| {d_ar:.1e} | max '
          f'rel dArea {d_area:.1e} | fit: rot {angle:.6f} deg, residual '
          f'{resid.max() * KM:.4f} km'
          + (f'  [expected: {note}]' if note else ''))


def main():
    if not GRID:
        raise SystemExit('set the grid: just cross-impl <grid>')
    # available_raw is also the stray-file police: an unrecognized file
    # in data/raw/ fails loudly here, same as in the metrics stage.
    raw = dict(metrics.available_raw())
    keys = sorted(k for k in raw if cache.parse_key(k)[0] == GRID)
    if len(keys) < 2:
        raise SystemExit(f'{GRID}: need raw from >=2 implementations, '
                         f'found {keys or "none"} (run `just gen` first)')
    fulls = {k: full_resolutions(k, raw[k]) for k in keys}
    for key_a, key_b in combinations(keys, 2):
        shared = sorted(fulls[key_a] & fulls[key_b])
        print(f'{key_a} vs {key_b} — full-enumeration resolutions '
              f'{shared or "NONE (nothing guaranteed-complete on both sides)"}')
        for res in shared:
            compare(key_a, key_b, res)


if __name__ == '__main__':
    main()
