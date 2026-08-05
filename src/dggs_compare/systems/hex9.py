"""hex9 — octahedral aperture-9 equal-area hexagonal DGGS (via the hex9 wheel).

hex9 addresses the WGS84 authalic sphere (geodetic->authalic conversion and an
equal-area warp are built into the binding; both are always on). Cells are
hexagons at every layer — the 12 topological pentagons (two per octahedral
vertex) are represented as 6-vertex rings with a reflected half-hex across the
seam, so no pentagon special-casing is needed here.

Cell ids: a hex9 bin uuid does not encode its layer, so the zone handle is the
tuple `(layer, uuid_bytes)` — hashable, and every consumer below unpacks it.
`cid_str` is the canonical keyed label ("<digits>.<key>"), which round-trips
exactly via hex9.parse_label; bare labels are ambiguous at split-hex bodies.

Cell edges are straight in the octahedral face plane, not sphere geodesics —
but the corner ring is still faithful for the solvers: the edges never bow
outside the corner-determined enclosing cone, so no `stats_ring` is needed.
Validated by scripts/validate_corners_hex9.py (corners vs densify<=3 rings
against a densify=5 ground truth): geometric max |dAR| is ~3e-9 at L2 and
falls off ~9x per layer; the ~1e-7 residuals at fine layers are csar solver
noise (identical across densify levels), inside the 1e-6 gap tolerance.
"""

import hex9
import numpy as np


def resolutions():
    return range(hex9.lmax() + 1)   # layers 0..30


def num_cells(res):
    return 12 * 9 ** res            # octahedral aperture-9: 12, 108, 972, ...


def _u(z):
    """The (layer, bytes) zone handle's uuid as the uint8[16] hex9 expects."""
    return np.frombuffer(z[1], dtype=np.uint8)


def cell_at(res, lat, lng):
    full = hex9.encode(np.array([lng], float), np.array([lat], float))
    return res, hex9.bin(full, res)[0].tobytes()


def cid_str(z):
    return hex9.label(_u(z), z[0], True)   # keyed label: exact round-trip


def cell_boundary(z):
    # hex9 winds its rings clockwise seen from outside; sparea needs CCW
    # (a CW ring reads as the 4pi-complement polygon), hence the reversal.
    ring = hex9.cell(_u(z), z[0])   # closed (lng, lat) ring
    return [(float(la), float(lo)) for lo, la in ring[-2::-1]]


def enumerate_cells(res):
    # The 12 layer-0 base cells, then each one's canonical 9^res-cell
    # partition of the layer (owned_cells streams bins base-cell by base-cell,
    # so memory stays flat at the full-enumeration resolutions).
    base = hex9.grid(-180.0, -90.0, 180.0, 90.0, 0)[0]
    for u0 in base:
        bins, _curves = hex9.owned_cells(u0, res)
        for b in bins:
            yield res, b.tobytes()
