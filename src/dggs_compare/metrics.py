"""Metrics stage: raw geometry parquet -> the published tables.

Binding-free by construction: this module reads what the runner wrote
(data/raw/) and never imports a DGGS binding. It computes the per-cell
metrics with the shared solvers (csar AR, sparea area), so every published
number carries ONE solver provenance regardless of which env generated the
geometry.

Also computes each implementation's convergence residuals from the
runner's density-0/dense vertex pairs (max |dAR| per gate level) — the
measured evidence that density-0 vertex lists are faithful solver inputs —
and records them in every final table's metadata.
"""

import json
import time
from importlib.metadata import version as _pkg_version

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from . import config, stats
from .cache import (BATCH, DATA_DIR, SCHEMA, VERTS_TYPE, open_ring,
                    open_writer)
from .runner import RAW_DIR


def _solver_metadata():
    meta = {'gap_tol': repr(config.GAP_TOL),
            'csar_method': config.CSAR_METHOD}
    for pkg in ('csar', 'sparea'):
        try:
            meta[f'version_{pkg}'] = _pkg_version(pkg)
        except Exception:
            pass
    return {k.encode(): v.encode() for k, v in meta.items()}


def _ar(verts):
    import csar
    r = csar.solve(csar.to_vec3(verts, geo='latlng_deg'), geo='vec3')
    return r.aspect_ratio if isinstance(r, csar.Converged) else None


def convergence_residuals(key):
    """{res: max |dAR|} from `key`'s density-0/dense pairs."""
    table = pq.read_table(RAW_DIR / f'{key}_convergence.parquet')
    out = {}
    for res, v0, vd in zip(table['res'].to_pylist(),
                           table['verts'].to_pylist(),
                           table['verts_dense'].to_pylist()):
        a0, ad = _ar(open_ring(v0)), _ar(open_ring(vd))
        if a0 is None or ad is None:
            continue
        out[res] = max(out.get(res, 0.0), abs(a0 - ad))
    return out


def build(key, res, extra_meta=None, out_dir=None):
    """One final table for `key` = '{grid}-{impl}' at `res`."""
    t0 = time.perf_counter()
    out_dir = DATA_DIR if out_dir is None else out_dir
    raw_path = RAW_DIR / f'{key}_r{res}.parquet'
    raw = pq.ParquetFile(raw_path)
    grid = raw.metadata.metadata[b'grid'].decode()

    # Carry the raw provenance forward, minus parquet's automatic
    # ARROW:schema key — that one describes the RAW schema, and embedding
    # it in the final file would break type reconstruction on read
    # (readers would lose the fixed_size_list vertex type).
    meta = {k: v for k, v in (raw.metadata.metadata or {}).items()
            if k != b'ARROW:schema'}
    meta.update(_solver_metadata())
    meta.update(extra_meta or {})

    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f'{key}_r{res}.parquet'
    writer = open_writer(path, SCHEMA, meta)
    dnc = 0
    n = 0
    try:
        for batch in raw.iter_batches(batch_size=BATCH):
            cids = batch['cid'].to_pylist()
            verts, ars, areas = [], [], []
            for vlist in batch['verts'].to_pylist():
                latlng = open_ring(vlist)
                ar, area = stats.cell_stats(latlng)
                verts.append(latlng)
                ars.append(ar)
                areas.append(area)
            dnc += int(np.isnan(ars).sum())
            n += len(cids)
            writer.write_table(pa.table({
                'dggs': pa.array([grid] * len(cids), pa.string()),
                'res': pa.array([res] * len(cids), pa.int32()),
                'cid': pa.array(cids, pa.string()),
                'verts': pa.array(verts, VERTS_TYPE),
                'ar': pa.array(ars, pa.float64()),
                'area': pa.array(areas, pa.float64()),
            }, schema=SCHEMA))
    finally:
        writer.close()
    kb = path.stat().st_size / 1024
    print(f'[{key} r{res:<2}] {n:>8} cells (DNC {dnc}) -> {path.name} '
          f'({kb:.0f} KiB) [{time.perf_counter() - t0:.0f}s]', flush=True)


def available_raw():
    """Sorted (key, [res, ...]) pairs present in data/raw/."""
    import re
    pat = re.compile(r'^(.+)_r(\d+)\.parquet$')
    found = {}
    for p in RAW_DIR.glob('*_r*.parquet'):
        if (m := pat.match(p.name)):
            found.setdefault(m.group(1), []).append(int(m.group(2)))
    return sorted((k, sorted(rs)) for k, rs in found.items())


def build_all(out_dir=None):
    """Finals for everything in data/raw/, with convergence residuals
    computed first and stamped into each table's metadata."""
    for key, res_list in available_raw():
        residuals = convergence_residuals(key)
        print(f'[{key}] convergence max |dAR| per level: '
              + '  '.join(f'r{r}={v:.1e}' for r, v in sorted(residuals.items())))
        extra = {b'convergence_max_dar': json.dumps(residuals).encode()}
        for res in res_list:
            build(key, res, extra_meta=extra, out_dir=out_dir)
