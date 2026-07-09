"""ISEA7H — Snyder equal-area icosahedral aperture-7 hex (via DGGAL)."""

from dggs_compare.dggal_engine import Adapter

_adapter = Adapter('ISEA7H')


def adapter():
    """The live-engine Adapter (edge refinement, point->cell streams, ...)."""
    return _adapter


def resolutions():
    return range(_adapter.max_level() + 1)


def num_cells(res):
    return _adapter.count(res)


def cell_at(res, lat, lng):
    return _adapter.zone_at(res, lng, lat)


def cid_str(z):
    return _adapter.cid_str(z)


def cell_boundary(z):
    return _adapter.cell_boundary(z)


def enumerate_cells(res):
    yield from _adapter.enumerate(res)
