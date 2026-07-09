"""Build the web viewer's static data from the Parquet tables.

A cheap column reshape — the stats were computed at generation time, so
nothing is solved here. Emits into web/out/ (gitignored):

- histograms.json — for every (system, resolution): a FIXED fine-bin
  histogram (NBINS bins over [1, AMAX] + one overflow bin) plus summary
  stats. The fixed grid lets the page re-aggregate to any coarser bin width
  in the browser; AMAX is a high global percentile so the working
  resolutions never reach overflow.
- globe/{sys}_r{res}_{pos.f32,idx.u32,ar.f32,ids.json} — ajglobe's native
  flat-binary polygon format, for the coarse resolutions only (largest res
  per system with <= GLOBE_MAX_CELLS cells, and everything below).
- manifest.json — what exists, per-system web colors/labels, the solve
  tolerance, and the shared globe AR max.
"""

import json
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
from matplotlib.colors import to_hex

from . import cache, config

WEB_OUT = Path(__file__).resolve().parents[2] / 'web' / 'out'
GLOBE_DIR = WEB_OUT / 'globe'

NBINS = 256               # fixed fine bins for the stored histograms
AMAX_PCT = 99.99          # global percentile that sets the fixed grid's top edge
GLOBE_MAX_CELLS = 80_000  # a system's globe shows the largest res at/under this


def systems():
    """Config-ordered systems that have tables on disk."""
    on_disk = set(cache.available_systems())
    return [s for s in config.TARGET_RES if s in on_disk]


def build_histograms():
    """Fixed-grid histograms + stats per (system, res). The grid is shared
    across every system so the page can overlay and re-bin on one axis."""
    ars = {s: {res: cache.load_columns(s, res, ['ar'])['ar']
               for res in cache.available_resolutions(s)}
           for s in systems()}
    allars = np.concatenate([a[~np.isnan(a)] for byres in ars.values()
                             for a in byres.values() if a.size])
    amax = float(np.percentile(allars, AMAX_PCT))
    edges = np.linspace(1.0, amax, NBINS + 1)

    data = {}
    for s, byres in ars.items():
        data[s] = {}
        for res, cell_ars in byres.items():
            a = cell_ars[~np.isnan(cell_ars)]
            if not a.size:
                continue
            data[s][str(res)] = {
                'counts': np.histogram(a, bins=edges)[0].astype(int).tolist(),
                'n': int(a.size),
                'dnc': int(np.isnan(cell_ars).sum()),
                'min': float(a.min()),
                'median': float(np.median(a)),
                'p99': float(np.percentile(a, 99)),
                'max': float(a.max()),
            }
    return {
        'edges': [float(x) for x in edges],  # counts[i] in [edges[i], edges[i+1])
        'nbins': NBINS,
        'amax': amax,
        'data': data,
    }


def globe_resolutions(s):
    """The coarse resolutions to render for system `s`: the contiguous prefix
    up to the last res with <= GLOBE_MAX_CELLS cells. Stops at the FIRST
    over-cap res — past the target resolution the tables hold N_SMALL
    *sampled* cells, which also fall under the cap but are sparse scatter,
    not globe coverage."""
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
    """histograms.json + globe binaries + manifest.json under web/out/."""
    WEB_OUT.mkdir(parents=True, exist_ok=True)

    print('building histograms...')
    hist = build_histograms()
    (WEB_OUT / 'histograms.json').write_text(json.dumps(hist))
    print(f'wrote histograms.json (amax={hist["amax"]:.4f}, {NBINS} bins)')

    print('building globe binaries...')
    globe_avail, globe_max = build_globe()

    manifest = {
        'systems': systems(),
        'colors': {s: to_hex(config.SYS_COLOR[s]) for s in systems()},
        'labels': {s: s.upper() for s in systems()},
        'res_prefix': {s: config.RES_PREFIX[s] for s in systems()},
        'target_res': {s: config.TARGET_RES[s] for s in systems()},
        'hist_res': {s: sorted(int(r) for r in hist['data'].get(s, {}))
                     for s in systems()},
        'globe_res': globe_avail,
        'globe_ar_max': globe_max,
        'gap_tol': config.GAP_TOL,
    }
    (WEB_OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    print(f'wrote manifest.json (globe AR range 1..{globe_max:.4f})')
