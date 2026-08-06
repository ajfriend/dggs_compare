"""hex9 — octahedral aperture-9 equal-area hexagonal DGGS (via the hex9 wheel).

hex9 addresses the WGS84 authalic sphere (geodetic->authalic conversion and an
equal-area warp are built into the binding; both are always on). Cells are
hexagons at every layer — the 12 topological pentagons (two per octahedral
vertex) are represented as 6-vertex rings with a reflected half-hex across the
seam, so no pentagon special-casing is needed here.

Zone handles are `(layer, uuid_bytes)` tuples: a hex9 bin uuid does not
encode its layer, and `cid_strs` (which gets no `res`) needs it for
labeling. The cid is the canonical keyed label ("<digits>.<key>"), which
round-trips exactly via hex9.parse_label — bare labels are ambiguous at
split-hex bodies.

Cell edges are straight in the octahedral face plane, not sphere geodesics —
but the corner ring is still faithful for the solvers (the edges never bow
outside the corner-determined enclosing cone), so no `stats_rings` override
is needed. `refined_boundaries` traces the true hex9 edges via the
binding's native densification, so `just validate-corners` checks
corners-vs-refined for hex9 like every other system.
"""

import hex9
import numpy as np

# The binding supports layers 0..30, but aperture 9 descends so fast that
# r19 cells (~0.3 cm^2) are already smaller than ANY other system's finest
# (s2 L30 ~ 0.74 cm^2), and past r19 the f64 solver floor takes over: r20 is
# ~96% DNC and exposes a csar edge-case crash. Nothing beyond r19 is
# numerically meaningful, so the pipeline stops there.
MAX_RES = 19


def resolutions():
    return range(MAX_RES + 1)


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
    # One batch label() call (a scalar call per zone costs ~2x). The batch
    # form takes a single layer; callers pass single-resolution batches.
    layer = zones[0][0]
    assert all(z[0] == layer for z in zones)
    arr = np.frombuffer(b''.join(b for _, b in zones), np.uint8).reshape(-1, 16)
    return hex9.label(arr, layer, True)


def _ring(z, d=0):
    # hex9 winds its rings clockwise seen from outside; sparea needs CCW
    # (a CW ring reads as the 4pi-complement polygon): the row slice drops
    # the closing vertex and reverses, the column flip swaps (lng, lat) ->
    # (lat, lng), all before the one C-speed tolist().
    return hex9.cell(_u(z), z[0], d)[-2::-1, ::-1].tolist()


def boundaries(res, zones):
    return [_ring(z) for z in zones]


def refined_boundaries(res, zones, refine):
    # smallest d with 3^d - 1 densification points per edge >= refine
    d = 1
    while 3 ** d - 1 < refine:
        d += 1
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
