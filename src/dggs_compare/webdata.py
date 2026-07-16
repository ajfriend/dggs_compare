"""Build the static site's globe data from the Parquet tables.

A cheap column reshape — the stats were computed at generation time, so
nothing is solved here. The site's distribution plots are the matplotlib PNGs
from `scripts/survey.py`; this module only emits the globe layer, into
web/out/ (gitignored):

- globe/{sys}_r{res}_{pos.f32,idx.u32,ar.f32,ids.json} — ajglobe's native
  flat-binary polygon format, for the coarse resolutions only (largest res
  per system with <= GLOBE_MAX_CELLS cells, and everything below). The site
  draws one globe per system at that system's largest such resolution.
- manifest.json — the data-release tag, per-system web colors/labels, the
  solve tolerance, the globe resolutions available, and the shared globe AR
  max (so every globe colors on one comparable scale).
"""

import json
import os
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
from matplotlib.colors import to_hex

from . import cache, config

WEB_OUT = Path(__file__).resolve().parents[2] / 'web' / 'out'
GLOBE_DIR = WEB_OUT / 'globe'

GLOBE_MAX_CELLS = 80_000  # a system's globe shows the largest res at/under this


def systems():
    """Config-ordered systems that have tables on disk."""
    on_disk = set(cache.available_systems())
    return [s for s in config.TARGET_RES if s in on_disk]


def globe_resolutions(s):
    """The coarse resolutions to render for system `s`: the contiguous prefix
    up to the last res with <= GLOBE_MAX_CELLS cells. Stops at the FIRST
    over-cap res — past it the tables hold N_CELLS *sampled* cells, which also
    fall under the cap but are sparse scatter, not globe coverage."""
    out = []
    for res in cache.available_resolutions(s):
        nrows = pq.ParquetFile(cache.table_path(s, res)).metadata.num_rows
        if nrows > GLOBE_MAX_CELLS:
            break
        out.append(res)
    return out


def build_globe():
    """Write ajglobe's flat binaries per coarse (system, res). Returns
    {system: [res, ...]} and the shared AR max over all globe cells."""
    GLOBE_DIR.mkdir(parents=True, exist_ok=True)
    avail, globe_max = {}, 1.0
    for s in systems():
        res_list = globe_resolutions(s)
        avail[s] = res_list
        for res in res_list:
            cols = cache.load_columns(s, res, ['cid', 'verts', 'ar'])
            cell_ars = cols['ar']
            finite = cell_ars[~np.isnan(cell_ars)]
            if finite.size:
                globe_max = max(globe_max, float(finite.max()))
            pos, starts = [], [0]
            for latlng in cols['verts']:
                pos.append(np.asarray(latlng, dtype='<f4')[:, ::-1])  # -> [lng, lat]
                starts.append(starts[-1] + len(latlng))
            stem = GLOBE_DIR / f'{s}_r{res}'
            np.concatenate(pos).tofile(f'{stem}_pos.f32')
            np.asarray(starts, dtype='<u4').tofile(f'{stem}_idx.u32')
            cell_ars.astype('<f4').tofile(f'{stem}_ar.f32')
            Path(f'{stem}_ids.json').write_text(json.dumps(cols['cid']))
            print(f'  globe {s} r{res}: {len(cols["cid"])} cells -> {stem.name}_*')
    return avail, globe_max


def build_all():
    """globe binaries + manifest.json under web/out/."""
    WEB_OUT.mkdir(parents=True, exist_ok=True)

    print('building globe binaries...')
    globe_avail, globe_max = build_globe()

    manifest = {
        # Which data release these artifacts were built from (set by pages.yml
        # from the workflow's `tag` input); empty for a local build.
        'tag': os.environ.get('DGGS_DATA_TAG', ''),
        'systems': systems(),
        'colors': {s: to_hex(config.SYS_COLOR[s]) for s in systems()},
        'labels': {s: s.upper() for s in systems()},
        'res_prefix': {s: config.RES_PREFIX[s] for s in systems()},
        'target_res': {s: config.TARGET_RES[s] for s in systems()},
        'globe_res': globe_avail,
        'globe_ar_max': globe_max,
        'gap_tol': config.GAP_TOL,
    }
    (WEB_OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    print(f'wrote manifest.json (globe AR range 1..{globe_max:.4f})')
