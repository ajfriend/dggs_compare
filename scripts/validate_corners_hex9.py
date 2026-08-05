"""Validate corners-only AR for hex9: hex9 edges are straight in the
octahedral face plane, not sphere geodesics, so in principle the 6-corner
ring could under-sweep the true boundary. Ground truth is a densify=5 ring
(1458 points); for each layer this reports the enclosing-cone AR error of
densify 0..3 rings against it. Measured: the geometric error peaks at ~3e-9
(L2) and the fine-layer residuals are csar solver noise (identical across
densify levels) — corners-only is faithful, so systems/hex9.py has no
stats_ring. Run when bumping the hex9 wheel.

Run with:  uv run scripts/validate_corners_hex9.py
No CLI args (project convention).
"""
import hex9
import numpy as np

import csar

from dggs_compare import stats

SEED = 0xC0FFEE
K = 200
RES_LIST = [0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15]
TRUTH_D = 5            # ground-truth densification (6*3^5 = 1458 points)
CAND_D = (0, 1, 2, 3)  # candidate solver rings


def ar(latlng):
    r = csar.solve(csar.to_vec3(latlng, geo='latlng_deg'), geo='vec3')
    return r.aspect_ratio if isinstance(r, csar.Converged) else np.nan


def ring(u, res, d):
    pts = hex9.cell(u, res, d)                 # closed (lng, lat)
    return [(float(la), float(lo)) for lo, la in pts[:-1]]


def zones_at(res, rng):
    if 12 * 9 ** res <= K:
        base = hex9.grid(-180.0, -90.0, 180.0, 90.0, 0)[0]
        return [b for u0 in base for b in hex9.owned_cells(u0, res)[0]]
    pts = stats.sample_uniform_lnglat(K * 2, rng)
    full = hex9.encode(pts[:, 0].copy(), pts[:, 1].copy())
    bins = hex9.bin(full, res)
    return list({b.tobytes(): b for b in bins}.values())[:K]


rng = np.random.default_rng(SEED)
hdr = ' '.join(f'{f"max|dAR| d{d}":>12}' for d in CAND_D)
print(f'{"layer":>5} {"cells":>6} {hdr}')
worst = {d: 0.0 for d in CAND_D}
for res in RES_LIST:
    zs = zones_at(res, rng)
    errs = {d: 0.0 for d in CAND_D}
    for u in zs:
        truth = ar(ring(u, res, TRUTH_D))
        if np.isnan(truth):
            continue
        for d in CAND_D:
            a = ar(ring(u, res, d))
            if not np.isnan(a):
                errs[d] = max(errs[d], abs(a - truth))
    for d in CAND_D:
        worst[d] = max(worst[d], errs[d])
    cols = ' '.join(f'{errs[d]:>12.2e}' for d in CAND_D)
    print(f'{res:>5} {len(zs):>6} {cols}')

print('\noverall worst per densify: ' +
      ', '.join(f'd{d}={worst[d]:.2e}' for d in CAND_D))
print('corners-only CONFIRMED for hex9' if worst[0] < 1e-6
      else 'NOT negligible — add a stats_ring to systems/hex9.py')
