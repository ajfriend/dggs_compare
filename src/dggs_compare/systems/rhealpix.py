"""rHEALPix — equal-area HEALPix-derived quad DGGS (via DGGAL)."""

from dggs_compare.dggal_engine import Adapter

_adapter = Adapter('rHEALPix')


def adapter():
    """The live-engine Adapter (edge refinement, point->cell streams, ...)."""
    return _adapter


def resolutions():
    return range(_adapter.max_level() + 1)


def num_cells(res):
    return _adapter.count(res)


def cells_at(res, points):
    return [_adapter.zone_at(res, lng, lat) for lat, lng in points]


def cid_str(z):
    return _adapter.cid_str(z)


def boundaries(zones):
    return [_adapter.cell_boundary(z) for z in zones]


def enumerate_cells(res):
    yield from _adapter.enumerate(res)
