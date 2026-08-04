"""IVEA3H — icosahedral vertex-oriented equal-area aperture-3 hex (via DGGAL).

Same layout as ISEA3H but with distortion spread smoothly instead of
concentrated in seams (the same ISEA/IVEA relationship as the 7H pair).
"""

from dggs_compare.dggal_engine import Adapter

_adapter = Adapter('IVEA3H')


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
