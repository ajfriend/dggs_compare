"""Build the static site's globe data from the Parquet tables.

A cheap column reshape — the stats were computed at generation time, so
nothing is measured here. The site's distribution plots are the matplotlib PNGs
from `scripts/survey.py`; this module only emits the globe layer, into
web/out/ (gitignored):

- globe/{sys}_r{res}_{pos.f32,idx.u32,ar.f32,area.f32,ids.json} — ajglobe's
  native flat-binary polygon format plus one f32 per cell per metric (AR, and
  area normalized by the resolution's exact mean 4*pi/N; the viewer derives
  each metric's shared color domain from these). One globe per system per
  ANCHOR: the viewer's Resolution menu offers every H3 resolution 0..
  config.GLOBE_H3_RES as a common cell size, and each system serves the
  resolution whose cell count is closest to that anchor's. Equal-area
  grids -> matching cell count matches average cell area, so this is an
  area match computed from the closed-form counts (no table reads).
  Anchor payloads shrink geometrically toward r0, so all the coarser
  levels together cost ~1/6 of the finest one.
- manifest.json — the data-release tag, per-system web colors/labels, the
  per-anchor globe resolutions per system (+ anchor labels for the menu),
  csar's gap tolerance, and the shared globe AR max (the viewer's
  provisional color domain while binaries load).
"""

import json
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
    """The available-on-disk resolution area-matched to `anchor_n` cells
    (config.count_match_res over the closed-form counts, so a *sampled*
    fine-resolution table can never masquerade as a coarse one)."""
    return config.count_match_res(sys, anchor_n, cache.available_resolutions(sys))


def _emit_globe(s, res):
    """Write one system's flat binaries at `res`; returns its finite AR max."""
    cols = cache.load_columns(s, res, ['cid', 'verts', 'ar', 'area'])
    cell_ars = cols['ar']
    finite = cell_ars[~np.isnan(cell_ars)]
    rel = cols['area'] / config.mean_cell_area(s, res)
    pos, starts = [], [0]
    for latlng in cols['verts']:
        pos.append(np.asarray(latlng, dtype='<f4')[:, ::-1])   # -> [lng, lat]
        starts.append(starts[-1] + len(latlng))
    stem = GLOBE_DIR / f'{s}_r{res}'
    np.concatenate(pos).tofile(f'{stem}_pos.f32')
    np.asarray(starts, dtype='<u4').tofile(f'{stem}_idx.u32')
    cell_ars.astype('<f4').tofile(f'{stem}_ar.f32')
    rel.astype('<f4').tofile(f'{stem}_area.f32')
    Path(f'{stem}_ids.json').write_text(json.dumps(cols['cid']))
    print(f'  globe {s} r{res}: {len(cols["cid"]):,} cells -> {stem.name}_*')
    return float(finite.max()) if finite.size else 1.0


def build_globe():
    """Write flat binaries for each system at every anchor level (H3
    r0..GLOBE_H3_RES, coarse to fine). Returns ({system: [res per
    anchor]}, [anchor label per anchor], shared AR max over all cells)."""
    GLOBE_DIR.mkdir(parents=True, exist_ok=True)
    anchors = range(config.GLOBE_H3_RES + 1)
    chosen = {s: [] for s in systems()}
    labels, globe_max, written = [], 1.0, set()
    for a in anchors:
        anchor_n = config.CELLS_PER_RES['h3'](a)
        km2 = config.mean_cell_area('h3', a) * config.SR2KM2
        labels.append(f'~{float(f"{km2:.3g}"):,.0f} km²/cell (≈H3 r{a})')
        for s in systems():
            res = globe_res_for(s, anchor_n)
            chosen[s].append(res)
            # neighboring anchors can pick the same resolution for a
            # slow-growing grid; the files are keyed by res, so emit once
            if (s, res) not in written:
                written.add((s, res))
                globe_max = max(globe_max, _emit_globe(s, res))
    return chosen, labels, globe_max


def build_all():
    """globe binaries + manifest.json under web/out/."""
    WEB_OUT.mkdir(parents=True, exist_ok=True)

    print(f'building globe binaries (anchors H3 r0..r{config.GLOBE_H3_RES})...')
    globe_res, anchor_labels, globe_max = build_globe()

    manifest = {
        # Which data release these artifacts were built from (set by pages.yml
        # from the workflow's `tag` input); empty for a local build.
        'tag': os.environ.get('DGGS_DATA_TAG', ''),
        'systems': systems(),
        'colors': {s: to_hex(config.SYS_COLOR[s]) for s in systems()},
        'labels': {s: s.upper() for s in systems()},
        'res_prefix': {s: config.RES_PREFIX[s] for s in systems()},
        'target_res': {s: config.TARGET_RES[s] for s in systems()},
        # {system: [resolution per anchor]} — area-matched, coarse to fine;
        # the viewer's Resolution menu indexes these (default: the last)
        'globe_res': globe_res,
        'globe_anchor_labels': anchor_labels,
        'globe_h3_res': config.GLOBE_H3_RES,
        'globe_ar_max': globe_max,
        'gap_tol': config.GAP_TOL,
    }
    (WEB_OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    print(f'wrote manifest.json (globe AR range 1..{globe_max:.4f})')
