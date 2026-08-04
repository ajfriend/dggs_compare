"""ISEA3H — Snyder equal-area icosahedral aperture-3 hex (via DGGAL)."""

from dggs_compare.dggal_engine import Adapter

_adapter = Adapter('ISEA3H')

# Even refinement counts land a vertex on the icosahedron-edge kink (odd
# counts miss it and reproduce the corners-only answer). 6 matches
# refine-40 AR to ~1e-8 (validated over 4.5k cells at levels 9/11/13)
# at half the solve cost and a third the FFI-array size of 20.
STATS_REFINE = 6


def stats_ring(z):
    """Solver-ring override for odd levels; None where corners suffice.

    Odd-level cells straddling an icosahedron edge kink there: the real
    boundary bulges past the corner hexagon by 100s of meters (points in
    the bulge classify to this cell via point->zone), putting corners-only
    AR off by up to ~4e-3. Even levels are corner-exact (<2e-8), as are
    all IVEA3H levels — vertex-oriented distortion doesn't kink.
    """
    if _adapter.level(z) % 2:
        return _adapter.refined_boundary(z, STATS_REFINE)
    return None


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
