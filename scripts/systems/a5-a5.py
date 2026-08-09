# /// script
# requires-python = ">=3.11,<3.14"
# dependencies = [
#     "a5_fast",
#     "dggs-compare",
# ]
#
# [tool.uv.sources]
# dggs-compare = { path = "../..", editable = true }
# ///
"""A5 — pentagonal equal-area DGGS (via the Rust/PyO3 a5_fast binding).

a5_fast's `cell_to_boundary` returns an adaptively densified vertex list
(321 points at res 0 down to 6 from res ~9 up). The extra points don't
change the enclosing-cone AR, so density 0 reduces the list to its corner
vertices (turning-angle peaks): 5 for the pentagons, 3 for the res-1
"quintant" triangles. Higher density slerps between a5's own sampled
vertices.

Run with:  uv run scripts/systems/a5-a5.py
"""

import a5_fast as a5
import numpy as np

from dggs_compare import runner, stats

# Exterior-angle threshold marking a corner (degrees). Densified edge
# points turn by small fractions of a degree; true corners by tens.
_TURN_DEG = 5.0


def _corners(latlng):
    """Reduce an open (lat, lng)-degree vertex list to its corners."""
    la, lo = np.radians(np.asarray(latlng, dtype=float)).T
    v = np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo),
                         np.sin(la)])
    e = np.roll(v, -1, axis=0) - v                    # chord to next vertex
    e /= np.linalg.norm(e, axis=1, keepdims=True)
    cosang = np.clip(np.einsum('ij,ij->i', np.roll(e, 1, axis=0), e), -1, 1)
    idx = np.nonzero(np.degrees(np.arccos(cosang)) > _TURN_DEG)[0]
    return [latlng[i] for i in idx] if len(idx) >= 3 else latlng


class A5:
    grid = 'a5'
    impl = 'a5'
    packages = ('a5_fast',)

    def resolutions(self):
        return range(31)            # a5 supports 0..30 (a5.MAX_RESOLUTION)

    def num_cells(self, res):
        return a5.get_num_cells(res)

    def cells_at(self, res, points):
        return [a5.lonlat_to_cell(lng, lat, res)   # int, hashable
                for lat, lng in points]

    def cid_strs(self, cells):
        return [a5.u64_to_hex(c) for c in cells]

    def _boundary(self, c):
        verts = a5.cell_to_boundary(c)   # closed (lng, lat) list, densified
        return _corners([(lat, lng) for lng, lat in verts[:-1]])

    def boundaries(self, res, cells, samples_per_edge=0):
        if samples_per_edge:
            return [stats.refine_geodesic(
                        [(lat, lng)
                         for lng, lat in a5.cell_to_boundary(c)[:-1]],
                        samples_per_edge)
                    for c in cells]
        return [self._boundary(c) for c in cells]

    def enumerate_cells(self, res):
        for c0 in a5.get_res0_cells():
            if res == 0:
                yield c0
            else:
                yield from a5.cell_to_children(c0, res)


if __name__ == '__main__':
    runner.generate(A5())
