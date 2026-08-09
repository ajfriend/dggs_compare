"""Generation runner: takes a GridImpl, writes raw geometry parquet.

This is the half of the pipeline that touches a DGGS binding, and it runs
inside the implementation script's own env (see scripts/systems/). It
writes, per resolution:

    data/raw/{grid}-{impl}_r{res}.parquet     cid + verts (density 0),
                                              cid-sorted

plus one small convergence artifact per implementation:

    data/raw/{grid}-{impl}_convergence.parquet
        (res, cid, verts, verts_dense) for K sampled cells at the gate
        levels — density 0 next to density CONV_SAMPLES, so the metrics
        stage can measure the sampling-density residual binding-free.

Selection regimes, seeds, and row order are identical to the classic
cache.build_table path (same seeds -> same cells), so ports verify by
column comparison. Every file's metadata records grid, impl, package
versions, the seed policy, and the resolutions covered.
"""

import time
from importlib.metadata import version as _pkg_version
from pathlib import Path

import numpy as np
import pyarrow as pa

from . import config, stats
from .cache import BATCH, MAX_DRAW_FACTOR, VERTS_TYPE, open_writer

RAW_DIR = Path(__file__).resolve().parents[2] / 'data' / 'raw'

RAW_SCHEMA = pa.schema([('cid', pa.string()), ('verts', VERTS_TYPE)])
CONV_SCHEMA = pa.schema([('res', pa.int32()), ('cid', pa.string()),
                         ('verts', VERTS_TYPE), ('verts_dense', VERTS_TYPE)])

# Convergence sampling (mirrors the classic scripts/convergence.py knobs):
# K cells per gate level, density 0 vs CONV_SAMPLES.
CONV_SEED = 0xC0FFEE
CONV_LEVELS = [0, 1, 2, 3, 5, 8, 11]
CONV_K = 300
CONV_SAMPLES = 40


def _metadata(impl):
    """Provenance recorded on every raw file (all-str, parquet-style)."""
    meta = {
        'grid': impl.grid,
        'impl': impl.impl,
        'seed': hex(config.SEED),
        'per_res_seed': str(config.PER_RES_SEED),
        'n_cells': str(config.N_CELLS),
        'resolutions': ','.join(str(r) for r in impl.resolutions()),
    }
    for pkg in ('dggs_compare', *getattr(impl, 'packages', ())):
        try:
            meta[f'version_{pkg}'] = _pkg_version(pkg)
        except Exception:
            pass
    return {k.encode(): v.encode() for k, v in meta.items()}


def _select_cells(impl, res):
    """The resolution's cell set: exactly N_CELLS cells (or every cell where
    fewer exist / config.FULL_RES demands complete coverage). Three regimes —
    see the config.N_CELLS comment. Returns (cells, mode)."""
    n = config.N_CELLS
    total = impl.num_cells(res)
    rng = np.random.default_rng(
        [config.SEED, res] if config.PER_RES_SEED else config.SEED)

    if total <= n or res in config.FULL_RES.get(impl.grid, ()):
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
        k = min(100_000, MAX_DRAW_FACTOR * n - drawn)
        pts = stats.sample_uniform_latlng(k, rng).tolist()
        hits = impl.cells_at(res, pts)
        for (lat, lng), c in zip(pts, hits):
            if c is None:
                # The engine couldn't resolve the point (rare deep-level
                # singular points) — draw again, and log the specimen.
                print(f'    unresolved point skip: ({lat:.6f}, {lng:.6f}) '
                      f'[{impl.grid}-{impl.impl} r{res}]', flush=True)
                continue
            if c not in seen:
                seen.add(c)
                cells.append(c)
                if len(cells) == n:
                    break
        drawn += k
    return cells, 'sample'


def _verts_arrays(vlists):
    return [[[float(la), float(ln)] for la, ln in v] for v in vlists]


def _write_res(impl, res, meta):
    """One resolution's raw table. Returns the cell count."""
    t0 = time.perf_counter()
    cells, mode = _select_cells(impl, res)
    cells = sorted(zip(impl.cid_strs(cells), cells), key=lambda cc: cc[0])
    key = f'{impl.grid}-{impl.impl}'
    path = RAW_DIR / f'{key}_r{res}.parquet'
    writer = open_writer(path, RAW_SCHEMA, meta)
    try:
        for lo in range(0, len(cells), BATCH):
            chunk = cells[lo:lo + BATCH]
            cids, clist = zip(*chunk)
            verts = _verts_arrays(impl.boundaries(res, clist))
            writer.write_table(pa.table(
                {'cid': pa.array(cids, pa.string()),
                 'verts': pa.array(verts, VERTS_TYPE)}, schema=RAW_SCHEMA))
    finally:
        writer.close()
    kb = path.stat().st_size / 1024
    print(f'[{key} r{res:<2}] {mode:>6} {len(cells):>8} cells '
          f'-> {path.name} ({kb:.0f} KiB) '
          f'[{time.perf_counter() - t0:.0f}s]', flush=True)
    return len(cells)


def _write_convergence(impl, meta):
    """Density-0 + dense vertex pairs for the gate levels, so the metrics
    stage can measure the sampling-density residual without the binding."""
    rng = np.random.default_rng(CONV_SEED)
    key = f'{impl.grid}-{impl.impl}'
    path = RAW_DIR / f'{key}_convergence.parquet'
    writer = open_writer(path, CONV_SCHEMA, meta)
    try:
        for level in CONV_LEVELS:
            if level not in impl.resolutions():
                continue
            cells = _conv_cells(impl, level, rng)
            cids = impl.cid_strs(cells)
            base = _verts_arrays(impl.boundaries(level, cells))
            dense = _verts_arrays(
                impl.boundaries(level, cells, CONV_SAMPLES))
            writer.write_table(pa.table(
                {'res': pa.array([level] * len(cells), pa.int32()),
                 'cid': pa.array(cids, pa.string()),
                 'verts': pa.array(base, VERTS_TYPE),
                 'verts_dense': pa.array(dense, VERTS_TYPE)},
                schema=CONV_SCHEMA))
    finally:
        writer.close()
    print(f'[{key}] convergence pairs -> {path.name}', flush=True)


def _conv_cells(impl, level, rng):
    total = impl.num_cells(level)
    if total <= CONV_K:
        return list(impl.enumerate_cells(level))
    if total <= 4 * CONV_K:
        cells = list(impl.enumerate_cells(level))
        idx = rng.choice(len(cells), CONV_K, replace=False)
        return [cells[i] for i in idx]
    seen, out = set(), []
    drawn = 0
    while len(out) < CONV_K:
        if drawn >= 60 * CONV_K:
            raise RuntimeError(
                f'{drawn:,} draws yielded only {len(out)}/{CONV_K} cells')
        pts = stats.sample_uniform_latlng(CONV_K * 4, rng).tolist()
        drawn += len(pts)
        for c in impl.cells_at(level, pts):
            if c is not None and c not in seen:
                seen.add(c)
                out.append(c)
                if len(out) >= CONV_K:
                    break
    return out


def generate(impl, only=None):
    """Write the raw geometry artifacts for `impl` (all resolutions, or
    the subset in `only`)."""
    t0 = time.perf_counter()
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    meta = _metadata(impl)
    res_list = [r for r in impl.resolutions() if only is None or r in only]
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
    print(f'[{impl.grid}-{impl.impl}] {len(res_list)} resolutions in '
          f'{time.perf_counter() - t0:.0f}s', flush=True)
