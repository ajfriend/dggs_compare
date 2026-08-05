"""ISEA4T — Snyder equal-area icosahedral aperture-4 triangle (via DGGRID).

The first DGGRID-backed system: batch subprocess calls instead of
per-cell FFI (see dggrid_engine) — the batch-first registry contract
means the module is stateless: each call maps directly onto one or two
DGGRID invocations. Zone handles are (res, seqnum) tuples since DGGRID
SEQNUM addresses are only unique within a resolution.

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
    return [None if s is None else (res, s)
            for s in _engine.cells_at(res, points)]


def cid_str(z):
    res, seq = z
    return str(seq)


def boundaries(zones):
    """One clipped generation per resolution present in `zones`."""
    by_res = {}
    for res, seq in zones:
        by_res.setdefault(res, []).append(seq)
    rings = {}
    for res, seqs in by_res.items():
        for seq, ring in _engine.boundaries(res, seqs).items():
            rings[(res, seq)] = ring
    return [rings[z] for z in zones]


def enumerate_cells(res):
    for seq in _engine.enumerate_ids(res):
        yield (res, seq)
