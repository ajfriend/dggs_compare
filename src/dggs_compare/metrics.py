"""Metrics stage: raw geometry parquet -> the published tables.

Binding-free by construction: this module reads what the runner wrote
(data/raw/) and never imports a DGGS binding. It computes the per-cell
metrics with the shared solvers (csar AR, sparea area), so every published
number carries ONE solver provenance regardless of which env generated the
geometry.

Also the admission gate: each implementation's convergence residuals
(max |dAR|, density 0 vs dense, from the runner's vertex pairs) are
measured here at the published solver settings, stamped into every final
table's metadata, and checked against config.CONV_TOL — a grid whose
density-0 vertex lists are not faithful solver inputs fails loudly unless
it is in config.CONV_EXPECTED_RED with a documented reason.
"""

import json
import time
from importlib.metadata import version as _pkg_version

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from . import cache, config, runner, stats
from .cache import BATCH, DATA_DIR, SCHEMA, open_ring, open_writer


def _solver_metadata():
    meta = {'gap_tol': repr(config.GAP_TOL),
            'csar_method': config.CSAR_METHOD}
    for pkg in ('csar', 'sparea'):
        try:
            meta[f'{runner.META_VERSION_PREFIX}{pkg}'] = _pkg_version(pkg)
        except Exception:
            pass
    return {k.encode(): v.encode() for k, v in meta.items()}


def _forwarded(raw_meta):
    """The runner-vocabulary subset of a raw file's metadata (META_KEYS +
    version_*); parquet's automatic keys and anything else stay behind."""
    keep = {k.encode() for k in runner.META_KEYS}
    prefix = runner.META_VERSION_PREFIX.encode()
    return {k: v for k, v in raw_meta.items()
            if k in keep or k.startswith(prefix)}


def convergence_residuals(key):
    """{res: max |dAR|} from `key`'s density-0/dense pairs, at the
    published solver settings (stats.ar)."""
    table = pq.read_table(runner.conv_path(key))
    out = {}
    for res, v0, vd in zip(table['res'].to_pylist(),
                           table['verts'].to_pylist(),
                           table['verts_dense'].to_pylist()):
        a0, ad = stats.ar(open_ring(v0)), stats.ar(open_ring(vd))
        if a0 is None or ad is None:
            continue
        out[res] = max(out.get(res, 0.0), abs(a0 - ad))
    return out


def build(key, res, extra_meta=None, out_dir=None):
    """One final table for `key` = '{grid}-{impl}' at `res`."""
    t0 = time.perf_counter()
    out_dir = DATA_DIR if out_dir is None else out_dir
    raw = pq.ParquetFile(runner.raw_path(key, res))
    raw_meta = raw.metadata.metadata
    grid = raw_meta[b'grid'].decode()
    # The filename and the metadata must agree on the identity — a copied
    # or hand-renamed raw file must not publish under the wrong key.
    assert key == cache.key(grid, raw_meta[b'impl'].decode()), key

    meta = _forwarded(raw_meta)
    meta.update(_solver_metadata())
    meta.update(extra_meta or {})

    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / cache.table_name(key, res)
    writer = open_writer(path, SCHEMA, meta)
    dnc = 0
    n = 0
    t_solve = 0.0
    try:
        for batch in raw.iter_batches(batch_size=BATCH):
            ts = time.perf_counter()
            ars, areas = [], []
            for vlist in batch['verts'].to_pylist():
                ar, area = stats.cell_stats(open_ring(vlist))
                ars.append(ar)
                areas.append(area)
            t_solve += time.perf_counter() - ts
            dnc += int(np.isnan(ars).sum())
            rows = len(batch)
            n += rows
            # cid/verts pass through untouched (the contract guarantees
            # open vertex lists, so the solver input above IS the stored
            # geometry); only the metric columns are new.
            writer.write_table(pa.table({
                'dggs': pa.array([grid] * rows, pa.string()),
                'res': pa.array([res] * rows, pa.int32()),
                'cid': batch['cid'],
                'verts': batch['verts'],
                'ar': pa.array(ars, pa.float64()),
                'area': pa.array(areas, pa.float64()),
            }, schema=SCHEMA))
    finally:
        writer.close()
    kb = path.stat().st_size / 1024
    t = time.perf_counter() - t0
    # solve = the per-cell csar/sparea loop (incl. the arrow->python vlist
    # decode); the residual is raw read + arrow build + zstd-19 write.
    print(f'[{key} r{res:<2}] {n:>8} cells (DNC {dnc}) -> {path.name} '
          f'({kb:.0f} KiB) [{t:.1f}s: solve {t_solve:.1f} = '
          f'{runner._us_per_cell(t_solve, n)}]', flush=True)


def available_raw():
    """Sorted (key, [res, ...]) pairs present in data/raw/. A file that
    matches neither the table scheme nor '{key}_convergence.parquet' is a
    loud error — a stray file must not silently join (or be silently
    excluded from) the published artifact."""
    found = {}
    bad = []
    for p in runner.RAW_DIR.glob('*.parquet'):
        if (parsed := cache.parse_table_name(p.name)):
            grid, impl, res = parsed
            found.setdefault(cache.key(grid, impl), []).append(res)
        elif not p.name.endswith('_convergence.parquet'):
            bad.append(p.name)
    if bad:
        raise RuntimeError(f'unrecognized files in data/raw/: {sorted(bad)}')
    return sorted((k, sorted(rs)) for k, rs in found.items())


def build_all():
    """Finals for everything in data/raw/: the admission gate first (a
    grid over CONV_TOL fails unless carved out), then the tables, with the
    residuals stamped into each table's metadata."""
    for key, res_list in available_raw():
        residuals = convergence_residuals(key)
        worst = max(residuals.values(), default=0.0)
        grid, _ = cache.parse_key(key)
        print(f'[{key}] convergence max |dAR| per level: '
              + '  '.join(f'r{r}={v:.1e}'
                          for r, v in sorted(residuals.items())))
        if worst >= config.CONV_TOL:
            reason = config.CONV_EXPECTED_RED.get(grid)
            if reason is None:
                raise RuntimeError(
                    f'{key}: convergence residual {worst:.1e} >= '
                    f'{config.CONV_TOL:g} — density-0 vertex lists are not '
                    f'faithful solver inputs for this grid')
            print(f'[{key}] over tolerance, expected: {reason}')
        extra = {b'convergence_max_dar': json.dumps(residuals).encode()}
        for res in res_list:
            build(key, res, extra_meta=extra)
