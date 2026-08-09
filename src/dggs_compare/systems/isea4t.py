"""ISEA4T — Snyder equal-area icosahedral aperture-4 triangle (via DGGRID).

The first DGGRID-backed system: batch subprocess calls instead of
per-cell FFI (see dggrid_engine) — the batch-first registry contract
means the module is stateless: each call maps directly onto one or two
DGGRID invocations. Cells are plain SEQNUM ints; they are only unique
within a resolution, which is why the contract passes `res` alongside.

MAX_RES 28: DGGRID accepts finer, but 20*4^r overflows uint64 seqnums
past r30 (observed: r31 ids come out smaller than r30's), and r28 cells
are already ~6 cm^2.

Density-0 stats are exact here (spot-checked to ~1e-12 against densely
sampled boundaries): the aperture-4 triangle lattice keeps icosahedron
edges ON cell edges at every level, so there are no distortion vertices.
"""

from dggs_compare.dggrid_engine import Engine

_engine = Engine('ISEA4T')

MAX_RES = 28


def resolutions():
    return range(MAX_RES + 1)


def num_cells(res):
    return 20 * 4 ** res


def cells_at(res, points):
    return _engine.cells_at(res, points)


# Zero-padded to the max seqnum's width: cache.build_table sorts rows by
# cid TEXT for spatially coherent Parquet pages, and variable-width
# decimals would sort '10' < '9'.
_CID_WIDTH = len(str(num_cells(MAX_RES)))


def cid_strs(cells):
    return [f'{c:0{_CID_WIDTH}d}' for c in cells]


def boundaries(res, cells, samples_per_edge=0):
    return _engine.boundaries(res, cells, samples_per_edge)


def enumerate_cells(res):
    yield from _engine.enumerate_ids(res)
