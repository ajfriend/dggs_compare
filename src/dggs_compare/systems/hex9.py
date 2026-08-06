"""hex9 — octahedral aperture-9 equal-area hexagonal DGGS (via the hex9 wheel).

hex9 addresses the WGS84 authalic sphere (geodetic->authalic conversion and an
equal-area warp are built into the binding; both are always on). Cells are
hexagons at every layer — the 12 topological pentagons (two per octahedral
vertex) are represented as 6-vertex rings with a reflected half-hex across the
seam, so no pentagon special-casing is needed here.

Zone handles are `(layer, uuid_bytes)` tuples: a hex9 bin uuid does not
encode its layer, and `cid_strs` (which gets no `res`) needs it for
labeling. The cid is the canonical keyed label ("<digits>.<key>") — its
hierarchical digits make the tables' text sort spatially coherent, it is
fixed-width within a resolution, and it round-trips exactly via
hex9.parse_label (bare labels are ambiguous at split-hex bodies).

Cell edges are straight in the octahedral face plane, not sphere geodesics —
but the corner ring is still faithful for the solvers (the edges never bow
outside the corner-determined enclosing cone), so no `stats_rings` override
is needed. `refined_boundaries` traces the true hex9 edges via the
binding's native densification, so `just validate-corners` checks
corners-vs-refined for hex9 like every other system.
"""

import math

import hex9
import numpy as np


def resolutions():
    return range(hex9.lmax() + 1)   # layers 0..30


def num_cells(res):
    return 12 * 9 ** res            # octahedral aperture-9: 12, 108, 972, ...


def _u(z):
    """The (layer, bytes) zone handle's uuid as the uint8[16] hex9 expects."""
    return np.frombuffer(z[1], dtype=np.uint8)


def cells_at(res, points):
    pts = np.asarray(points, dtype=float).reshape(-1, 2)   # (lat, lng) rows
    full = hex9.encode(np.ascontiguousarray(pts[:, 1]),    # lng
                       np.ascontiguousarray(pts[:, 0]))    # lat
    return [(res, b.tobytes()) for b in hex9.bin(full, res)]


def cid_strs(zones):
    return [hex9.label(_u(z), z[0], True) for z in zones]


def _ring(z, d=0):
    # hex9 winds its rings clockwise seen from outside; sparea needs CCW
    # (a CW ring reads as the 4pi-complement polygon), hence the reversal.
    ring = hex9.cell(_u(z), z[0], d)   # closed (lng, lat) ring
    return [(float(la), float(lo)) for lo, la in ring[-2::-1]]


def boundaries(res, zones):
    return [_ring(z) for z in zones]


def refined_boundaries(res, zones, refine):
    # densify=d inserts 3^d - 1 points per edge along the true hex9 edge
    # (straight in the octahedral face plane — hex9's own edge model).
    d = max(1, math.ceil(math.log(refine + 1) / math.log(3)))
    return [_ring(z, d) for z in zones]


def enumerate_cells(res):
    # The 12 layer-0 base cells, then each one's canonical 9^res-cell
    # partition of the layer (owned_cells streams bins base-cell by base-cell,
    # so memory stays flat at the full-enumeration resolutions).
    base = hex9.grid(-180.0, -90.0, 180.0, 90.0, 0)[0]
    for u0 in base:
        bins, _curves = hex9.owned_cells(u0, res)
        for b in bins:
            yield res, b.tobytes()
