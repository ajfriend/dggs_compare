"""Build the static site's globe data from the Parquet tables.

A cheap column reshape — the stats were computed at generation time, so
nothing is solved here. The site's distribution plots are the matplotlib PNGs
from `scripts/survey.py`; this module only emits the globe layer, into
web/out/ (gitignored):

- globe/{sys}_r{res}_{pos.f32,idx.u32,ar.f32,ids.json} — ajglobe's native
  flat-binary polygon format. One globe per system, all at a common cell size:
  the resolution whose cell count is closest to H3's at config.GLOBE_H3_RES.
  Equal-area grids -> matching cell count matches average cell area, so this is
  an area match computed from the closed-form counts (no table reads).
- manifest.json — the data-release tag, per-system web colors/labels, the
  chosen globe resolution per system, the solve tolerance, and the shared globe
  AR max (so every globe colors on one comparable scale).
"""

import json
import math
import os
from pathlib import Path

import numpy as np
from matplotlib.colors import to_hex

from . import cache, config

WEB_OUT = Path(__file__).resolve().parents[2] / 'web' / 'out'
GLOBE_DIR = WEB_OUT / 'globe'


def systems():
    """Config-ordered systems that have tables on disk."""
    on_disk = set(cache.available_systems())
    return [s for s in config.TARGET_RES if s in on_disk]


def globe_res_for(sys, anchor_n):
    """The resolution whose (true) cell count is closest — in log-ratio — to
    `anchor_n`, over the system's available resolutions. Area-match: these are
    equal-area grids, so avg cell area = 4*pi*R^2 / N(res), and matching N
    matches cell size. Uses the closed-form counts (config.CELLS_PER_RES), so a
    *sampled* fine-resolution table can never masquerade as a coarse one."""
    n_of = config.CELLS_PER_RES[sys]
    return min(cache.available_resolutions(sys),
               key=lambda r: abs(math.log(n_of(r) / anchor_n)))


def build_globe():
    """Write ajglobe's flat binaries for each system's area-matched globe.
    Returns {system: res} and the shared AR max over those globes' cells."""
    GLOBE_DIR.mkdir(parents=True, exist_ok=True)
    anchor_n = config.CELLS_PER_RES['h3'](config.GLOBE_H3_RES)
    chosen, globe_max = {}, 1.0
    for s in systems():
        res = globe_res_for(s, anchor_n)
        chosen[s] = res
        cols = cache.load_columns(s, res, ['cid', 'verts', 'ar'])
        cell_ars = cols['ar']
        finite = cell_ars[~np.isnan(cell_ars)]
        if finite.size:
            globe_max = max(globe_max, float(finite.max()))
        pos, starts = [], [0]
        for latlng in cols['verts']:
            pos.append(np.asarray(latlng, dtype='<f4')[:, ::-1])   # -> [lng, lat]
            starts.append(starts[-1] + len(latlng))
        stem = GLOBE_DIR / f'{s}_r{res}'
        np.concatenate(pos).tofile(f'{stem}_pos.f32')
        np.asarray(starts, dtype='<u4').tofile(f'{stem}_idx.u32')
        cell_ars.astype('<f4').tofile(f'{stem}_ar.f32')
        Path(f'{stem}_ids.json').write_text(json.dumps(cols['cid']))
        print(f'  globe {s} r{res}: {len(cols["cid"]):,} cells '
              f'(anchor {anchor_n:,}) -> {stem.name}_*')
    return chosen, globe_max


def build_all():
    """globe binaries + manifest.json under web/out/."""
    WEB_OUT.mkdir(parents=True, exist_ok=True)

    print(f'building globe binaries (H3 r{config.GLOBE_H3_RES} anchor)...')
    globe_res, globe_max = build_globe()

    manifest = {
        # Which data release these artifacts were built from (set by pages.yml
        # from the workflow's `tag` input); empty for a local build.
        'tag': os.environ.get('DGGS_DATA_TAG', ''),
        'systems': systems(),
        'colors': {s: to_hex(config.SYS_COLOR[s]) for s in systems()},
        'labels': {s: s.upper() for s in systems()},
        'res_prefix': {s: config.RES_PREFIX[s] for s in systems()},
        'target_res': {s: config.TARGET_RES[s] for s in systems()},
        'globe_res': globe_res,              # {system: resolution} — area-matched
        'globe_h3_res': config.GLOBE_H3_RES,
        'globe_ar_max': globe_max,
        'gap_tol': config.GAP_TOL,
    }
    (WEB_OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    print(f'wrote manifest.json (globe AR range 1..{globe_max:.4f})')
