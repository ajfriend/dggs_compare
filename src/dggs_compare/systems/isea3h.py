"""ISEA3H — Snyder equal-area icosahedral aperture-3 hex (via DGGAL).

KNOWN APPROXIMATION (issue #25): stats use density-0 vertices, like every
other grid, but isea3h's odd-level cells straddling an icosahedron edge
kink there — the real boundary bulges past the density-0 hexagon by 100s
of meters (ground-truthed via point->cell), so density-0 AR is off by up
to ~4e-3 for those cells, and the convergence check reports this line red
by design. The honest fix (solving on densely sampled vertices) is
disabled because dggal 0.0.6's runtime degrades under that path's
per-cell array churn — revisit per issue #25 (upstream
distortion-vertices API / dggal 0.0.7 / direct-C vertices). Even levels
are exact at density 0 (<2e-8), as is all of IVEA3H.
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
    return _adapter.cells_at(res, points)


def cid_strs(cells):
    return _adapter.cid_strs(cells)


def boundaries(res, cells, samples_per_edge=0):
    return _adapter.boundaries(cells, samples_per_edge)


def enumerate_cells(res):
    yield from _adapter.enumerate(res)
