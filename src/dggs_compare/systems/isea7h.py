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


def cells_at(res, points):
    return _adapter.cells_at(res, points)


def cid_strs(cells):
    return _adapter.cid_strs(cells)


def boundaries(res, cells, samples_per_edge=0):
    return _adapter.boundaries(cells, samples_per_edge)


def enumerate_cells(res):
    yield from _adapter.enumerate(res)
