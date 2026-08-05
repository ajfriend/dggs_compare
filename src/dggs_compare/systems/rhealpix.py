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


def cid_strs(zones):
    return [_adapter.cid_str(z) for z in zones]


def boundaries(res, zones):
    return [_adapter.cell_boundary(z) for z in zones]


def refined_boundaries(res, zones, refine):
    return [_adapter.refined_boundary(z, refine) for z in zones]


def enumerate_cells(res):
    yield from _adapter.enumerate(res)
