"""The data artifact: one Parquet table per (system, resolution), one row per
cell, with the computed stats bundled alongside the geometry.

Schema:

    dggs   string                            constant per file, e.g. 'isea7h'
    res    int32                             constant per file
    cid    string                            cell id text
    verts  list<fixed_size_list<double, 2>>  ring of [lat, lng] degrees (open)
    ar     float64                           enclosing-cone aspect ratio (csar);
                                             NaN = did-not-converge at GAP_TOL
    area   float64                           spherical area, steradians (sparea)

Downstream consumers (survey, dnc-check, web data, DuckDB/pandas users) read
columns — no solver in the read path. Provenance (seed policy, budgets,
solver settings, library versions) travels in each file's Parquet metadata.

Files: data/cells/{dggs}_r{res}.parquet at the repo root (gitignored;
published as GitHub data releases, one asset per table). Each holds exactly
N_CELLS cells, enumerated in full where the resolution has fewer.
"""

import os
import time
from importlib.metadata import version as _pkg_version
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from . import config, registry, stats

# Repo-internal library: src/dggs_compare/cache.py -> repo root is parents[2].
DATA_DIR = Path(__file__).resolve().parents[2] / 'data' / 'cells'

# fixed_size_list(2): each vertex is exactly [lat, lng] deg; the outer list is
# the variable-length ring (6 for hexagons, 5 for the pentagons, 4 for quads).
VERTS_TYPE = pa.list_(pa.list_(pa.float64(), 2))
# Leaf column path of the nested verts doubles (for per-column encoding).
VERTS_LEAF = 'verts.list.element.list.element'
SCHEMA = pa.schema([
    ('dggs', pa.string()),
    ('res', pa.int32()),
    ('cid', pa.string()),
    ('verts', VERTS_TYPE),
    ('ar', pa.float64()),
    ('area', pa.float64()),
])


def _provenance():
    """Config + versions recorded in every table's file metadata."""
    meta = {
        'seed': hex(config.SEED),
        'per_res_seed': str(config.PER_RES_SEED),
        'n_cells': str(config.N_CELLS),
        'gap_tol': repr(config.GAP_TOL),
        'csar_method': config.CSAR_METHOD,
    }
    for pkg in ('csar', 'sparea', 'h3', 's2sphere', 'a5_fast', 'dggal',
                'hex9'):
        try:
            meta[f'version_{pkg}'] = _pkg_version(pkg)
        except Exception:
            pass
    return {k.encode(): v.encode() for k, v in meta.items()}


def open_ring(ring):
    """Drop a closing vertex if the ring repeats its first point."""
    if len(ring) >= 2 and tuple(ring[0]) == tuple(ring[-1]):
        return ring[:-1]
    return ring


def table_path(dggs, res):
    """Canonical Parquet path for a (system, resolution) table."""
    return DATA_DIR / f'{dggs}_r{res}.parquet'


def _existing_path(dggs, res):
    path = table_path(dggs, res)
    if not path.exists():
        raise FileNotFoundError(
            f'{path} not found — build the tables first with `just gen`, or '
            f'download a data release with `just fetch-data <tag>`.')
    return path


BATCH = 50_000        # cells per streamed row group (bounds build memory)
MAX_DRAW_FACTOR = 60  # safety cap on point draws in the sample regime

# Per-column encodings, applied to whichever of these columns a schema has.
# BYTE_STREAM_SPLIT packs the float64s ~28% smaller, losslessly (shared
# sign/exponent bytes compress; random mantissa bytes stay out of the way).
# DELTA_BYTE_ARRAY prefix-compresses the cid-sorted ids.
_ENCODINGS = {'cid': 'DELTA_BYTE_ARRAY',
              VERTS_LEAF: 'BYTE_STREAM_SPLIT',
              'ar': 'BYTE_STREAM_SPLIT',
              'area': 'BYTE_STREAM_SPLIT'}


def open_writer(path, schema, metadata, compression_level=19):
    """A ParquetWriter with the house encodings for `schema`.

    zstd level 19 (the default, for PUBLISHED tables) is ~free to read
    (decode speed is level-independent; measured 0.04s vs 0.05s on the
    1.18M-row table) and pays only at write time (0.4s -> 9s there). It
    matters most on full-enumeration tables, where cid-sorted neighbors
    have near-identical vertex bytes that long-range matching exploits:
    -23% there, ~-5% on sampled tables. Intermediates that are written
    once and read once (data/raw/) should pass a low level instead.
    """
    names = set(schema.names)
    return pq.ParquetWriter(
        path, schema.with_metadata(metadata),
        compression='zstd',
        compression_level=compression_level,
        use_dictionary=[c for c in ('dggs', 'res') if c in names],
        column_encoding={k: v for k, v in _ENCODINGS.items()
                         if k.split('.')[0] in names},
    )


def _select_cells(sysmod, dggs, res):
    """The resolution's cell set: exactly N_CELLS cells (or every cell where
    fewer exist / config.FULL_RES demands complete coverage). Three regimes —
    see the config.N_CELLS comment. Returns (cells, mode)."""
    n = config.N_CELLS
    total = sysmod.num_cells(res)
    rng = np.random.default_rng(
        [config.SEED, res] if config.PER_RES_SEED else config.SEED)

    if total <= n or res in config.FULL_RES.get(dggs, ()):
        return list(sysmod.enumerate_cells(res)), 'all'

    if total <= config.SUBSAMPLE_MAX_RATIO * n:
        cells = list(sysmod.enumerate_cells(res))
        idx = rng.choice(len(cells), n, replace=False)
        return [cells[i] for i in idx], 'subsam'

    seen, cells = set(), []
    drawn = 0
    while len(cells) < n:
        if drawn >= MAX_DRAW_FACTOR * n:
            raise RuntimeError(
                f'{dggs} r{res}: {drawn:,} draws yielded only '
                f'{len(cells):,}/{n:,} distinct cells')
        k = min(100_000, MAX_DRAW_FACTOR * n - drawn)
        pts = stats.sample_uniform_latlng(k, rng).tolist()
        hits = sysmod.cells_at(res, pts)
        for (lat, lng), c in zip(pts, hits):
            if c is None:
                # The engine couldn't resolve the point (DGGAL nullZone at
                # rare singular points) — draw again, and log the specimen:
                # each one is a concrete example for the upstream report.
                print(f'    nullZone skip: ({lat:.6f}, {lng:.6f}) '
                      f'[{dggs} r{res}]', flush=True)
                continue
            if c not in seen:
                seen.add(c)
                cells.append(c)
                if len(cells) == n:
                    break
        drawn += k
    return cells, 'sample'


def build_table(dggs, res):
    """Build the `(dggs, res)` table — geometry + stats in one pass — and
    stream it to Parquet in BATCH-cell row groups, so memory stays flat at
    any budget."""
    t0 = time.perf_counter()
    sysmod = registry.get(dggs)
    cells, mode = _select_cells(sysmod, dggs, res)

    # Sort by cid for a canonical, deterministic row order (independent of
    # sampling order): enables Parquet cid page-stats / range pushdown, and
    # lets DELTA_BYTE_ARRAY prefix-compress the sorted ids.
    cells = sorted(zip(sysmod.cid_strs(cells), cells), key=lambda cc: cc[0])

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    path = table_path(dggs, res)
    writer = open_writer(path, SCHEMA, _provenance())
    dnc = 0
    try:
        for lo in range(0, len(cells), BATCH):
            chunk = cells[lo:lo + BATCH]
            cids, clist = zip(*chunk)
            boundaries = sysmod.boundaries(res, clist)
            verts, ars, areas = [], [], []
            for vlist in boundaries:
                latlng = [[float(la), float(ln)] for la, ln in open_ring(vlist)]
                ar, area = stats.cell_stats(latlng)
                verts.append(latlng)
                ars.append(ar)
                areas.append(area)
            dnc += int(np.isnan(ars).sum())
            writer.write_table(pa.table({
                'dggs': pa.array([dggs] * len(cids), pa.string()),
                'res': pa.array([res] * len(cids), pa.int32()),
                'cid': pa.array(cids, pa.string()),
                'verts': pa.array(verts, VERTS_TYPE),
                'ar': pa.array(ars, pa.float64()),
                'area': pa.array(areas, pa.float64()),
            }, schema=SCHEMA))
    finally:
        writer.close()

    kb = os.path.getsize(path) / 1024
    print(f'[{dggs} r{res:<2}] {mode:>6} {len(cells):>8} cells '
          f'(DNC {dnc}) -> {path.name} ({kb:.0f} KiB) '
          f'[{time.perf_counter() - t0:.0f}s]', flush=True)
    return path


def build_system(dggs):
    """Build every resolution's table for one system."""
    t0 = time.perf_counter()
    sysmod = registry.get(dggs)
    res_list = list(sysmod.resolutions())
    # Live ETA for whoever is watching the CI log, weighted by CELL counts
    # (known exactly up front from the three-regime budget) — a per-
    # resolution average is useless here because the coarse resolutions are
    # nearly free and every deep one is a full N_CELLS build. Per-cell cost
    # is ~constant across a system's resolutions (boundary FFI + solve are
    # both per-cell flat; for DGGAL grids the boundary call is the bigger
    # share), so remaining-cells x observed rate is honest from the first
    # resolutions on. Known wobble: isea3h's odd levels cost ~4x its even
    # ones (refined stats rings).
    full = config.FULL_RES.get(dggs, ())
    cells = {r: sysmod.num_cells(r) if r in full
             else min(sysmod.num_cells(r), config.N_CELLS) for r in res_list}
    total_cells = sum(cells.values())
    done_cells = 0
    for i, res in enumerate(res_list):
        build_table(dggs, res)
        done = time.perf_counter() - t0
        done_cells += cells[res]
        if i + 1 < len(res_list):
            eta = done / done_cells * (total_cells - done_cells)
            print(f'    {i + 1}/{len(res_list)} resolutions '
                  f'({done_cells:,}/{total_cells:,} cells) in {done:.0f}s '
                  f'(~{eta:.0f}s to go)', flush=True)
    print(f'[{dggs}] all {len(res_list)} resolutions in '
          f'{time.perf_counter() - t0:.0f}s', flush=True)


def build_all():
    """Build every table for every system in the registry."""
    for name in registry.names():
        build_system(name)


# ----- readers -------------------------------------------------------------

def available_systems():
    """Sorted distinct DGGS names that have tables in data/cells/."""
    import re
    pat = re.compile(r'^(.+)_r\d+\.parquet$')
    names = {m.group(1) for p in DATA_DIR.glob('*_r*.parquet')
             if (m := pat.match(p.name))}
    return sorted(names)


def available_resolutions(dggs):
    """Sorted resolutions that have a `dggs` table on disk."""
    import re
    pat = re.compile(rf'^{re.escape(dggs)}_r(\d+)\.parquet$')
    res = [int(m.group(1)) for p in DATA_DIR.glob(f'{dggs}_r*.parquet')
           if (m := pat.match(p.name))]
    return sorted(res)


def load_columns(dggs, res, columns):
    """The requested columns of one table as a dict of numpy/python arrays.

    Float columns come back as numpy arrays; `cid` as a list of str; `verts`
    as a list of (M, 2) float arrays.
    """
    table = pq.read_table(_existing_path(dggs, res), columns=list(columns))
    out = {}
    for col in columns:
        if col == 'cid':
            out[col] = table['cid'].to_pylist()
        elif col == 'verts':
            out[col] = [np.asarray(v, dtype=float)
                        for v in table['verts'].to_pylist()]
        else:
            out[col] = table[col].to_numpy(zero_copy_only=False)
    return out


def load_cells(dggs, res):
    """Yield (cid, (M, 2) lat/lng array) — streaming, memory-flat."""
    path = _existing_path(dggs, res)
    for batch in pq.ParquetFile(path).iter_batches(columns=['cid', 'verts']):
        for row in batch.to_pylist():
            yield row['cid'], np.asarray(row['verts'], dtype=float)


def load_cells_sample(dggs, res, n):
    """Yield (cid, (M, 2) array) for up to `n` cells drawn at random.

    Rows are stored sorted by cid (spatially coherent), so a prefix would be
    a clustered patch of the globe; random indices give a representative
    sample. Reproducible via config.SEED.
    """
    path = _existing_path(dggs, res)
    table = pq.read_table(path, columns=['cid', 'verts'])
    total = table.num_rows
    if total > n:
        idx = np.random.default_rng(config.SEED).choice(total, n, replace=False)
        table = table.take(idx)
    for cid, verts in zip(table['cid'].to_pylist(), table['verts'].to_pylist()):
        yield cid, np.asarray(verts, dtype=float)
