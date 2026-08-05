"""ISEA3H — Snyder equal-area icosahedral aperture-3 hex (via DGGAL).

KNOWN APPROXIMATION (issue #25): stats are corners-only, like every other
grid, but isea3h's odd-level cells straddling an icosahedron edge kink
there — the real boundary bulges past the corner hexagon by 100s of
meters (ground-truthed via point->zone), so corners-only AR is off by up
to ~4e-3 for those cells, and `just validate-corners` reports this line
red by design. The honest fix (a `stats_rings` override of edge-refined
vertices, see the registry contract) is disabled because dggal 0.0.6's
runtime degrades under that path's per-cell array churn — revisit per
issue #25 (upstream distortion-vertices API / dggal 0.0.7 / direct-C
corners). Even levels are corner-exact (<2e-8), as is all of IVEA3H.
"""

from dggs_compare.dggal_engine import Adapter

_adapter = Adapter('ISEA3H')


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
