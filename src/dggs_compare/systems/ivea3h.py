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


def cells_at(res, points):
    return _adapter.cells_at(res, points)


def cid_strs(zones):
    return _adapter.cid_strs(zones)


def boundaries(res, zones):
    return _adapter.boundaries(zones)


def refined_boundaries(res, zones, refine):
    return _adapter.refined_boundaries(zones, refine)


def enumerate_cells(res):
    yield from _adapter.enumerate(res)
