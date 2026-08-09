"""Generation runner: takes a GridImpl, writes raw geometry parquet.

This is the half of the pipeline that touches a DGGS binding, and it runs
inside the implementation script's own env (see scripts/systems/). It
writes, per resolution:

    data/raw/{grid}-{impl}_r{res}.parquet     cid + verts (density 0),
                                              cid-sorted

plus one small convergence artifact per implementation: density-0 and
dense vertex pairs for K sampled cells at the gate levels, so the metrics
stage can measure the sampling-density residual binding-free.

Selection regimes, seeds, and row order are identical to the classic
cache.build_table path (same seeds -> same cells), so ports verify by
column comparison.
"""

import os
import time
from importlib.metadata import version as _pkg_version
from pathlib import Path

import numpy as np
import pyarrow as pa

from . import config, stats
from .cache import BATCH, VERTS_TYPE, open_writer

RAW_DIR = Path(__file__).resolve().parents[2] / 'data' / 'raw'
RAW_COMPRESSION_LEVEL = 3   # written once, read once — don't pay zstd 19
MAX_DRAW_FACTOR = 60        # safety cap on point draws in the sample regime

RAW_SCHEMA = pa.schema([('cid', pa.string()), ('verts', VERTS_TYPE)])
CONV_SCHEMA = pa.schema([('res', pa.int32()), ('cid', pa.string()),
                         ('verts', VERTS_TYPE), ('verts_dense', VERTS_TYPE)])

# The raw-metadata vocabulary: exactly these keys (plus version_*) are
# written by _metadata and forwarded into the published tables by the
# metrics stage. Anything else in a raw file's metadata stays behind.
META_KEYS = ('grid', 'impl', 'seed', 'per_res_seed', 'n_cells',
             'resolutions', 'conv_seed', 'conv_levels', 'conv_k',
             'conv_samples')
META_VERSION_PREFIX = 'version_'


# ----- artifact naming: the one home -------------------------------------
def key(grid, impl_name):
    """The artifact key '{grid}-{impl}'. Invertible only if neither half
    contains '-', so that is enforced here."""
    assert '-' not in grid and '-' not in impl_name, (grid, impl_name)
    return f'{grid}-{impl_name}'


def raw_path(key_, res):
    return RAW_DIR / f'{key_}_r{res}.parquet'


def conv_path(key_):
    return RAW_DIR / f'{key_}_convergence.parquet'


def _metadata(impl):
    """Provenance recorded on every raw file (all-str, parquet-style)."""
    meta = {
        'grid': impl.grid,
        'impl': impl.impl,
        'seed': hex(config.SEED),
        'per_res_seed': str(config.PER_RES_SEED),
        'n_cells': str(config.N_CELLS),
        'resolutions': ','.join(str(r) for r in impl.resolutions()),
        'conv_seed': hex(config.CONV_SEED),
        'conv_levels': ','.join(str(r) for r in config.CONV_LEVELS),
        'conv_k': str(config.CONV_K),
        'conv_samples': str(config.CONV_SAMPLES),
    }
    for pkg in ('dggs_compare', *impl.packages):
        try:
            meta[f'{META_VERSION_PREFIX}{pkg}'] = _pkg_version(pkg)
        except Exception:
            pass
    return {k.encode(): v.encode() for k, v in meta.items()}


def _sample_distinct(impl, res, n, rng, *, batch, log_skips):
    """`n` distinct cells via the three selection regimes (see the
    config.N_CELLS comment): enumerate-all, enumerate+subsample, or
    draw-until-n in `batch`-point draws. Returns (cells, mode)."""
    total = impl.num_cells(res)
    if total <= n:
        return list(impl.enumerate_cells(res)), 'all'

    if total <= config.SUBSAMPLE_MAX_RATIO * n:
        cells = list(impl.enumerate_cells(res))
        idx = rng.choice(len(cells), n, replace=False)
        return [cells[i] for i in idx], 'subsam'

    seen, cells = set(), []
    drawn = 0
    while len(cells) < n:
        if drawn >= MAX_DRAW_FACTOR * n:
            raise RuntimeError(
                f'{impl.grid}-{impl.impl} r{res}: {drawn:,} draws yielded '
                f'only {len(cells):,}/{n:,} distinct cells')
        k = min(batch, MAX_DRAW_FACTOR * n - drawn)
        pts = stats.sample_uniform_latlng(k, rng).tolist()
        hits = impl.cells_at(res, pts)
        for (lat, lng), c in zip(pts, hits):
            if c is None:
                # The engine couldn't resolve the point (rare deep-level
                # singular points) — draw again; each logged specimen is a
                # concrete example for an upstream report.
                if log_skips:
                    print(f'    unresolved point skip: ({lat:.6f}, '
                          f'{lng:.6f}) [{impl.grid}-{impl.impl} r{res}]',
                          flush=True)
                continue
            if c not in seen:
                seen.add(c)
                cells.append(c)
                if len(cells) == n:
                    break
        drawn += k
    return cells, 'sample'


def _select_cells(impl, res):
    """The resolution's cell set: exactly N_CELLS cells (or every cell
    where fewer exist / config.FULL_RES demands complete coverage)."""
    if res in config.FULL_RES.get(impl.grid, ()):
        return list(impl.enumerate_cells(res)), 'all'
    rng = np.random.default_rng(
        [config.SEED, res] if config.PER_RES_SEED else config.SEED)
    return _sample_distinct(impl, res, config.N_CELLS, rng,
                            batch=100_000, log_skips=True)


def _write_res(impl, res, meta):
    """One resolution's raw table."""
    t0 = time.perf_counter()
    cells, mode = _select_cells(impl, res)
    cells = sorted(zip(impl.cid_strs(cells), cells), key=lambda cc: cc[0])
    path = raw_path(key(impl.grid, impl.impl), res)
    writer = open_writer(path, RAW_SCHEMA, meta,
                         compression_level=RAW_COMPRESSION_LEVEL)
    try:
        for lo in range(0, len(cells), BATCH):
            chunk = cells[lo:lo + BATCH]
            cids, clist = zip(*chunk)
            writer.write_table(pa.table(
                {'cid': pa.array(cids, pa.string()),
                 'verts': pa.array(impl.boundaries(res, clist), VERTS_TYPE)},
                schema=RAW_SCHEMA))
    finally:
        writer.close()
    kb = path.stat().st_size / 1024
    print(f'[{path.stem} ] {mode:>6} {len(cells):>8} cells '
          f'({kb:.0f} KiB) [{time.perf_counter() - t0:.0f}s]', flush=True)


def _write_convergence(impl, meta):
    """Density-0 + dense vertex pairs for the gate levels."""
    rng = np.random.default_rng(config.CONV_SEED)
    path = conv_path(key(impl.grid, impl.impl))
    writer = open_writer(path, CONV_SCHEMA, meta,
                         compression_level=RAW_COMPRESSION_LEVEL)
    try:
        for level in config.CONV_LEVELS:
            if level not in impl.resolutions():
                continue
            cells, _ = _sample_distinct(impl, level, config.CONV_K, rng,
                                        batch=config.CONV_K * 4,
                                        log_skips=False)
            writer.write_table(pa.table(
                {'res': pa.array([level] * len(cells), pa.int32()),
                 'cid': pa.array(impl.cid_strs(cells), pa.string()),
                 'verts': pa.array(impl.boundaries(level, cells),
                                   VERTS_TYPE),
                 'verts_dense': pa.array(
                     impl.boundaries(level, cells, config.CONV_SAMPLES),
                     VERTS_TYPE)},
                schema=CONV_SCHEMA))
    finally:
        writer.close()
    print(f'[{path.stem}] convergence pairs written', flush=True)


def generate(impl):
    """Write the raw geometry artifacts for `impl` — all resolutions, or
    the subset named in the DGGS_COMPARE_RES env var (comma-separated; for
    quick local runs — CI always generates everything)."""
    t0 = time.perf_counter()
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    meta = _metadata(impl)
    res_list = list(impl.resolutions())
    if (only := os.environ.get('DGGS_COMPARE_RES')):
        wanted = {int(r) for r in only.split(',')}
        res_list = [r for r in res_list if r in wanted]
    # Live ETA weighted by CELL counts (coarse resolutions are nearly free;
    # every deep one is a full N_CELLS build).
    full = config.FULL_RES.get(impl.grid, ())
    cells = {r: impl.num_cells(r) if r in full
             else min(impl.num_cells(r), config.N_CELLS) for r in res_list}
    total_cells = sum(cells.values())
    done_cells = 0
    for i, res in enumerate(res_list):
        _write_res(impl, res, meta)
        done_cells += cells[res]
        done = time.perf_counter() - t0
        if i + 1 < len(res_list):
            eta = done / done_cells * (total_cells - done_cells)
            print(f'    {i + 1}/{len(res_list)} resolutions '
                  f'({done_cells:,}/{total_cells:,} cells) in {done:.0f}s '
                  f'(~{eta:.0f}s to go)', flush=True)
    _write_convergence(impl, meta)
    print(f'[{key(impl.grid, impl.impl)}] {len(res_list)} resolutions in '
          f'{time.perf_counter() - t0:.0f}s', flush=True)
