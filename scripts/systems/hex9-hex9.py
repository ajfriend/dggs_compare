# /// script
# requires-python = ">=3.11,<3.14"
# dependencies = [
#     "hex9>=2.3.1",
#     "dggs-compare",
# ]
#
# [tool.uv.sources]
# dggs-compare = { path = "../..", editable = true }
# ///
"""hex9 — octahedral aperture-9 equal-area hexagonal DGGS (via the hex9
wheel).

hex9 addresses the WGS84 authalic sphere (geodetic->authalic conversion
and an equal-area warp are built into the binding; both are always on).
Cells are hexagons at every layer — the 12 topological pentagons (two per
octahedral vertex) are represented as 6-vertex lists with a reflected
half-hex across the seam, so no pentagon special-casing is needed here.

Cell handles are `(layer, uuid_bytes)` tuples: a hex9 bin uuid does not
encode its layer, and `cid_strs` (which gets no `res`) needs it for
labeling. The cid is the canonical keyed label ("<digits>.<key>"), which
round-trips exactly via hex9.parse_label — bare labels are ambiguous at
split-hex bodies.

MAX_RES 19: the binding supports layers 0..30, but aperture 9 descends so
fast that r19 cells (~0.3 cm^2) are already smaller than ANY other
system's finest (s2 L30 ~ 0.74 cm^2), and past r19 the f64 solver floor
takes over: r20 is ~96% DNC and exposes a csar edge-case crash. Nothing
beyond r19 is numerically meaningful, so the pipeline stops there.

Run with:  uv run scripts/systems/hex9-hex9.py
"""

import hex9
import numpy as np

from dggs_compare import runner, stats

MAX_RES = 19


class Hex9:
    grid = 'hex9'
    impl = 'hex9'
    packages = ('hex9',)

    def resolutions(self):
        return range(MAX_RES + 1)

    def num_cells(self, res):
        return 12 * 9 ** res        # octahedral aperture-9: 12, 108, ...

    def _u(self, c):
        return np.frombuffer(c[1], dtype=np.uint8)

    def cells_at(self, res, points):
        pts = np.asarray(points, dtype=float).reshape(-1, 2)  # (lat, lng)
        full = hex9.encode(np.ascontiguousarray(pts[:, 1]),   # lng
                           np.ascontiguousarray(pts[:, 0]))   # lat
        return [(res, b.tobytes()) for b in hex9.bin(full, res)]

    def cid_strs(self, cells):
        # One batch label() call (a scalar call per cell costs ~2x). The
        # batch form takes a single layer; callers pass single-resolution
        # batches.
        layer = cells[0][0]
        assert all(c[0] == layer for c in cells)
        arr = np.frombuffer(b''.join(b for _, b in cells),
                            np.uint8).reshape(-1, 16)
        return hex9.label(arr, layer, True)

    def boundaries(self, res, cells, samples_per_edge=0):
        # hex9's densify=d samples 3^d - 1 vertices per edge along the
        # true hex9 edge (straight in the octahedral face plane): the
        # smallest d meeting the requested density. hex9 winds its vertex
        # lists clockwise seen from outside; sparea needs CCW (a CW list
        # reads as the 4pi-complement polygon): the row slice drops the
        # closing vertex and reverses, the column flip swaps (lng, lat) ->
        # (lat, lng), all before the one C-speed tolist().
        d = 0
        while 3 ** d - 1 < samples_per_edge:
            d += 1
        # hex9 addresses WGS84 geodetic coordinates; the authalic map is
        # how this system gets to the sphere.
        return stats.authalic_rings(
            [hex9.cell(self._u(c), c[0], d)[-2::-1, ::-1].tolist()
             for c in cells])

    def enumerate_cells(self, res):
        # The 12 layer-0 base cells, then each one's canonical 9^res-cell
        # partition of the layer (owned_cells streams bins base-cell by
        # base-cell, so memory stays flat at full-enumeration resolutions).
        base = hex9.grid(-180.0, -90.0, 180.0, 90.0, 0)[0]
        for u0 in base:
            bins, _curves = hex9.owned_cells(u0, res)
            for b in bins:
                yield res, b.tobytes()


if __name__ == '__main__':
    runner.generate(Hex9())
