"""The data artifact: one Parquet table per (system, resolution), one row per
cell, with the computed stats bundled alongside the geometry.

Schema:

    dggs   string                            constant per file, e.g. 'isea7h'
    res    int32                             constant per file
    cid    string                            cell id text
    verts  list<fixed_size_list<double, 2>>  ring of [lat, lng] degrees (open)
    ar     float64                           enclosing-cone aspect ratio (skar);
                                             NaN = did-not-converge at GAP_TOL
    area   float64                           spherical area, steradians (sparea)

Downstream consumers (survey, dnc-check, web data, DuckDB/pandas users) read
columns — no solver in the read path. Provenance (seed policy, budgets,
solver settings, library versions) travels in each file's Parquet metadata.

Files: data/cells/{dggs}_r{res}.parquet at the repo root (gitignored;
published as GitHub data releases). Each holds up to N_BIG cells at/below the
system's target resolution and N_SMALL above it, enumerated in full where the
resolution has fewer cells.
"""

import os
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
        'n_big': str(config.N_BIG),
        'n_small': str(config.N_SMALL),
        'gap_tol': repr(config.GAP_TOL),
        'skar_method': config.SKAR_METHOD,
    }
    for pkg in ('skar', 'sparea', 'h3', 's2sphere', 'a5_fast', 'dggal'):
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


def build_table(dggs, res):
    """Build the `(dggs, res)` table — geometry + stats in one pass — and
    write it to Parquet.

    Draws up to `n` cells (N_BIG at/below the system's target resolution,
    N_SMALL above). If the resolution has `num_cells(res) <= n`, every cell
    is enumerated (exact, complete); otherwise `n` uniform-on-sphere points
    are drawn (module SEED policy) and deduped to distinct cells.
    """
    sysmod = registry.get(dggs)
    n = config.N_BIG if res <= config.TARGET_RES[dggs] else config.N_SMALL
    if sysmod.num_cells(res) <= n:
        zones = list(sysmod.enumerate_cells(res))
        mode = 'all'
    else:
        rng = np.random.default_rng(
            [config.SEED, res] if config.PER_RES_SEED else config.SEED)
        seen, zones = set(), []
        for lng, lat in stats.sample_uniform_lnglat(n, rng):
            z = sysmod.cell_at(res, float(lat), float(lng))
            if z not in seen:
                seen.add(z)
                zones.append(z)
        mode = 'sample'

    # Sort by cid for a canonical, deterministic row order (independent of
    # sampling order): enables Parquet cid page-stats / range pushdown, and
    # lets DELTA_BYTE_ARRAY prefix-compress the sorted ids.
    rows = sorted(((sysmod.cid_str(z), open_ring(sysmod.cell_boundary(z)))
                   for z in zones), key=lambda cr: cr[0])
    cids, verts, ars, areas = [], [], [], []
    for cid, ring in rows:
        latlng = [[float(la), float(ln)] for la, ln in ring]
        ar, area = stats.cell_stats(latlng)
        cids.append(cid)
        verts.append(latlng)
        ars.append(ar)
        areas.append(area)

    table = pa.table({
        'dggs': pa.array([dggs] * len(cids), pa.string()),
        'res': pa.array([res] * len(cids), pa.int32()),
        'cid': pa.array(cids, pa.string()),
        'verts': pa.array(verts, VERTS_TYPE),
        'ar': pa.array(ars, pa.float64()),
        'area': pa.array(areas, pa.float64()),
    }, schema=SCHEMA).replace_schema_metadata(_provenance())

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    path = table_path(dggs, res)
    # BYTE_STREAM_SPLIT packs the float64s ~28% smaller, losslessly (shared
    # sign/exponent bytes compress; random mantissa bytes stay out of the way).
    # DELTA_BYTE_ARRAY prefix-compresses the sorted cids.
    pq.write_table(
        table, path,
        compression='zstd',
        use_dictionary=['dggs', 'res'],
        column_encoding={'cid': 'DELTA_BYTE_ARRAY',
                         VERTS_LEAF: 'BYTE_STREAM_SPLIT',
                         'ar': 'BYTE_STREAM_SPLIT',
                         'area': 'BYTE_STREAM_SPLIT'},
    )
    kb = os.path.getsize(path) / 1024
    dnc = int(np.isnan(ars).sum())
    print(f'[{dggs} r{res:<2}] {mode:>6} {len(cids):>7} cells '
          f'(DNC {dnc}) -> {path.name} ({kb:.0f} KiB)', flush=True)
    return path


def build_system(dggs):
    """Build every resolution's table for one system."""
    for res in registry.get(dggs).resolutions():
        build_table(dggs, res)


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
