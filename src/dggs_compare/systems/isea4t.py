"""ISEA4T — Snyder equal-area icosahedral aperture-4 triangle (via DGGRID).

The first DGGRID-backed system: batch subprocess calls instead of
per-cell FFI (see dggrid_engine) — the batch-first registry contract
means the module is stateless: each call maps directly onto one or two
DGGRID invocations. Zones are plain SEQNUM ints; they are only unique
within a resolution, which is why the contract passes `res` alongside.

MAX_RES 28: DGGRID accepts finer, but 20*4^r overflows uint64 seqnums
past r30 (observed: r31 ids come out smaller than r30's), and r28 cells
are already ~6 cm^2.

Corners-only stats are exact here (spot-checked to ~1e-12 against
densified boundaries): the aperture-4 triangle lattice keeps icosahedron
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


def cid_strs(zones):
    # Zero-padded to fixed width (max seqnum at MAX_RES is 19 digits):
    # cache.build_table sorts rows by cid TEXT for spatially coherent
    # Parquet pages, and variable-width decimals would sort '10' < '9'.
    return [f'{z:019d}' for z in zones]


def boundaries(res, zones):
    rings = _engine.boundaries(res, zones)
    return [rings[z] for z in zones]


def refined_boundaries(res, zones, refine):
    rings = _engine.boundaries(res, zones, refine)
    return [rings[z] for z in zones]


def enumerate_cells(res):
    yield from _engine.enumerate_ids(res)
