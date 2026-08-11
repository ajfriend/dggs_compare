"""The data artifact: one Parquet table per (grid, implementation,
resolution), one row per cell, with the computed stats bundled alongside
the geometry.

Schema:

    dggs   string                            the grid, e.g. 'isea7h'
    res    int32                             constant per file
    cid    string                            cell id text
    verts  list<fixed_size_list<double, 2>>  vertex list, [lat, lng] deg, open
    ar     float64                           enclosing-cone aspect ratio (csar);
                                             NaN = did-not-converge at GAP_TOL
    area   float64                           spherical area, steradians (sparea)
    irregular  bool                         implementation-DECLARED exceptional
                                            cells (a hex grid's 12 pentagons);
                                            see interface.py

Files: data/cells/{grid}-{impl}_r{res}.parquet (gitignored; published as
GitHub data releases, one asset per table). Written by the metrics stage
(see metrics.py); provenance travels in each file's Parquet metadata.

Downstream consumers (survey, dnc-check, web data, DuckDB/pandas users)
read columns — no measurement code in the read path. The readers here
are keyed by GRID; a grid with several implementations resolves through
config.PRIMARY_IMPL.
"""

import re
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from . import config

# Repo-internal library: src/dggs_compare/cache.py -> repo root is parents[2].
DATA_DIR = Path(__file__).resolve().parents[2] / 'data' / 'cells'

# fixed_size_list(2): each vertex is exactly [lat, lng] deg; the outer list is
# the variable-length vertex list (6 for hexagons, 5 for pentagons, 4 quads).
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
    ('irregular', pa.bool_()),
])

BATCH = 50_000        # cells per streamed row group (bounds build memory)

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


def open_ring(ring):
    """Drop a closing vertex if the vertex list repeats its first point."""
    if len(ring) >= 2 and tuple(ring[0]) == tuple(ring[-1]):
        return ring[:-1]
    return ring


# ----- artifact naming: THE one home ---------------------------------------
# An implementation is named by the key '{grid}-{impl}'; its tables are
# '{key}_r{res}.parquet' and its script is 'scripts/systems/{key}.py'. The
# key is invertible only because neither half may contain '-': composition
# asserts it, and every parse in the pipeline routes through KEY_RE.

KEY_RE = r'([^-]+)-([^-]+)'                # the two halves of a key
_TABLE_PAT = re.compile(rf'^{KEY_RE}_r(\d+)\.parquet$')


def key(grid, impl):
    """The artifact key '{grid}-{impl}'."""
    assert '-' not in grid and '-' not in impl, (grid, impl)
    return f'{grid}-{impl}'


def parse_key(key_):
    """(grid, impl) from a key; raises on a malformed one."""
    grid, impl = key_.split('-')
    return grid, impl


def table_name(key_, res):
    return f'{key_}_r{res}.parquet'


def parse_table_name(name):
    """(grid, impl, res) from a table filename, or None if it doesn't
    match the scheme."""
    m = _TABLE_PAT.match(name)
    return (m.group(1), m.group(2), int(m.group(3))) if m else None


def impl_table_path(grid, impl, res):
    return DATA_DIR / table_name(key(grid, impl), res)


def table_path(dggs, res):
    """Canonical Parquet path for `dggs` (a grid name) at `res`, resolved
    through the grid's primary implementation."""
    return impl_table_path(dggs, config.PRIMARY_IMPL[dggs], res)


def _existing_path(dggs, res, impl=None):
    path = impl_table_path(dggs, impl or config.PRIMARY_IMPL[dggs], res)
    if not path.exists():
        raise FileNotFoundError(
            f'{path} not found — build the tables first (`just gen` then '
            f'`just metrics`), or download a data release with '
            f'`just fetch-data <tag>`.')
    return path


# ----- readers -------------------------------------------------------------

def available_tables():
    """{(grid, impl): sorted [res, ...]} for every table in data/cells/.
    Files that don't match the {grid}-{impl}_r{res} scheme are ignored
    (not part of the artifact)."""
    found = {}
    for p in DATA_DIR.glob('*.parquet'):
        if (parsed := parse_table_name(p.name)):
            grid, impl, res = parsed
            found.setdefault((grid, impl), []).append(res)
    return {k: sorted(v) for k, v in found.items()}


def available_systems():
    """Sorted grids whose PRIMARY_IMPL tables are present."""
    tables = available_tables()
    return sorted(g for (g, i) in tables
                  if config.PRIMARY_IMPL.get(g) == i)


def available_resolutions(dggs, impl=None):
    """Sorted resolutions of `dggs` with a table on disk (the primary
    implementation's, unless `impl` says otherwise)."""
    return available_tables().get(
        (dggs, impl or config.PRIMARY_IMPL[dggs]), [])


def table_metadata(grid, impl, res):
    """One table's Parquet metadata, decoded to a str dict."""
    meta = pq.ParquetFile(impl_table_path(grid, impl, res)).metadata.metadata
    return {k.decode(): v.decode() for k, v in (meta or {}).items()}


def load_columns(dggs, res, columns, impl=None):
    """The requested columns of one table as a dict of numpy/python arrays.

    Float columns come back as numpy arrays; `cid` as a list of str; `verts`
    as a list of (M, 2) float arrays.
    """
    table = pq.read_table(_existing_path(dggs, res, impl), columns=list(columns))
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


def load_cells(dggs, res, impl=None):
    """Yield (cid, (M, 2) lat/lng array) — streaming, memory-flat."""
    path = _existing_path(dggs, res, impl)
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
