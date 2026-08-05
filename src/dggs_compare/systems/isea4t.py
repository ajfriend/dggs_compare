"""ISEA4T — Snyder equal-area icosahedral aperture-4 triangle (via DGGRID).

The first DGGRID-backed system: batch subprocess calls instead of
per-cell FFI (see dggrid_engine). Zone handles are (res, seqnum) tuples.
Implements the batch registry hooks (cells_at_batch, boundaries_batch) —
per-cell calls work too, but the pipeline prefers the batch forms, and
per-cell cell_at costs a whole subprocess.

MAX_RES 28: DGGRID accepts finer, but 20*4^r overflows uint64 seqnums
past r30 (observed: r31 ids come out smaller than r30's), and r28 cells
are already ~6 cm^2.
"""

from dggs_compare.dggrid_engine import Engine

_engine = Engine('ISEA4T')

MAX_RES = 28
_BOUNDARY_CACHE = {}      # (res, seqnum) -> ring, filled by boundaries_batch


def resolutions():
    return range(MAX_RES + 1)


def num_cells(res):
    return 20 * 4 ** res


def cell_at(res, lat, lng):
    return (res, _engine.cells_at(res, [(lat, lng)])[0])


def cells_at_batch(res, latlngs):
    """[(res, seqnum) or None, ...] for a batch of (lat, lng) points."""
    return [None if s is None else (res, s)
            for s in _engine.cells_at(res, latlngs)]


def cid_str(z):
    res, seq = z
    return str(seq)


def cell_boundary(z):
    ring = _BOUNDARY_CACHE.get(z)
    if ring is None:
        boundaries_batch([z])
        ring = _BOUNDARY_CACHE[z]
    return ring


def boundaries_batch(zones):
    """Rings for a batch of zones (one clipped generation per resolution
    present in the batch). The module-level cache holds only the latest
    batch — a whole resolution's rings would be gigabytes."""
    _BOUNDARY_CACHE.clear()
    by_res = {}
    for z in zones:
        by_res.setdefault(z[0], []).append(z[1])
    for res, seqs in by_res.items():
        for seq, ring in _engine.boundaries(res, seqs).items():
            _BOUNDARY_CACHE[(res, seq)] = ring
    return {z: _BOUNDARY_CACHE[z] for z in zones}


def enumerate_cells(res):
    for seq in _engine.enumerate_ids(res):
        yield (res, seq)
