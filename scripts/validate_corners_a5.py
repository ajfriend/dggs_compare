"""Validate corners-only AR for A5: enclosing-cone AR from geometrically
detected corner vertices vs the dense boundary a5_fast returns (which is
ground truth: any edge bowing outside the corner hull would show up there).
Also reports the corner-count distribution per resolution. Run when
bumping a5_fast.

Run with:  uv run scripts/validate_corners_a5.py
No CLI args (project convention).
"""
import a5_fast as a5
import numpy as np

import skar

from dggs_compare import stats

SEED = 0xC0FFEE
K = 200
RES_LIST = [0, 1, 2, 5, 9, 14, 20, 25, 30]
TURN_DEG = 5.0            # exterior-angle threshold marking a corner


def ar(latlng):
    r = skar.solve(skar.to_vec3(latlng, geo='latlng_deg'), geo='vec3')
    return r.aspect_ratio if isinstance(r, skar.Converged) else np.nan


def corner_idx(ring_latlng):
    """Indices of high-turning points of an open (lat,lng) ring."""
    v = np.asarray(skar.to_vec3(ring_latlng, geo='latlng_deg'), dtype=float)
    e = np.roll(v, -1, axis=0) - v                    # chord to next vertex
    e /= np.linalg.norm(e, axis=1, keepdims=True)
    prev = np.roll(e, 1, axis=0)
    cosang = np.clip(np.einsum('ij,ij->i', prev, e), -1, 1)
    turn = np.degrees(np.arccos(cosang))              # exterior angle at each pt
    return np.nonzero(turn > TURN_DEG)[0], turn


rng = np.random.default_rng(SEED)
print(f'{"res":>4} {"cells":>6} {"pts":>5} {"corners":>9} {"max|dAR|":>10} {"max_rel":>10}')
worst = 0.0
for res in RES_LIST:
    if a5.get_num_cells(res) <= K:
        zones = a5.get_res0_cells() if res == 0 else list(
            c for c0 in a5.get_res0_cells() for c in a5.cell_to_children(c0, res))
    else:
        pts = stats.sample_uniform_lnglat(K * 2, rng)
        zones = list({a5.lonlat_to_cell(lng, lat, res) for lng, lat in pts})[:K]
    max_abs = max_rel = 0.0
    ncorners = set()
    npts = None
    for z in zones:
        ring = a5.cell_to_boundary(z)                 # closed (lng, lat) ring
        npts = len(ring)
        dense = [(lat, lng) for lng, lat in ring[:-1]]
        idx, _ = corner_idx(dense)
        ncorners.add(len(idx))
        five = [dense[i] for i in idx]
        a_d, a_c = ar(dense), ar(five)
        if np.isnan(a_d) or np.isnan(a_c):
            continue
        d = abs(a_d - a_c)
        max_abs = max(max_abs, d)
        max_rel = max(max_rel, d / a_d)
    worst = max(worst, max_abs)
    print(f'{res:>4} {len(zones):>6} {npts:>5} {str(sorted(ncorners)):>9} '
          f'{max_abs:>10.2e} {max_rel:>10.2e}')

print(f'\noverall max |dAR| = {worst:.3e}')
print('corners-only CONFIRMED for A5' if worst < 1e-3 else 'NOT negligible — keep dense boundaries')
