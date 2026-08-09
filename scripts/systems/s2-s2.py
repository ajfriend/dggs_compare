# /// script
# requires-python = ">=3.11,<3.14"
# dependencies = [
#     "s2sphere",
#     "dggs-compare",
# ]
#
# [tool.uv.sources]
# dggs-compare = { path = "../..", editable = true }
# ///
"""S2 — Google's quadtree-on-a-cube DGGS (via s2sphere).

Already on the sphere: s2 cell edges lie in planes through the origin
(constant-u/v lines on the cube), i.e. great circles, so higher sampling
density is slerp between the vertices — definitional, since s2sphere's
own containment uses these great circles.

Run with:  uv run scripts/systems/s2-s2.py
"""

import s2sphere

from dggs_compare import runner, stats


class S2:
    grid = 's2'
    impl = 's2'
    packages = ('s2sphere',)

    def resolutions(self):
        return range(31)            # s2sphere supports levels 0..30

    def num_cells(self, res):
        return 6 * 4 ** res

    def _cell_at(self, res, lat, lng):
        # from_lat_lng yields a level-30 leaf; walk up to the cell at `res`.
        leaf = s2sphere.CellId.from_lat_lng(
            s2sphere.LatLng.from_degrees(lat, lng))
        return leaf.parent(res).id()   # int, hashable

    def cells_at(self, res, points):
        return [self._cell_at(res, lat, lng) for lat, lng in points]

    def cid_strs(self, cells):
        return [format(c, '016x') for c in cells]

    def _boundary(self, c):
        cell = s2sphere.Cell(s2sphere.CellId(c))
        verts = []
        for i in range(4):
            ll = s2sphere.LatLng.from_point(cell.get_vertex(i))
            verts.append((ll.lat().degrees, ll.lng().degrees))
        return verts

    def boundaries(self, res, cells, samples_per_edge=0):
        if samples_per_edge:
            return [stats.refine_geodesic(self._boundary(c), samples_per_edge)
                    for c in cells]
        return [self._boundary(c) for c in cells]

    def enumerate_cells(self, res):
        for cid in s2sphere.CellId.walk(res):
            yield cid.id()


if __name__ == '__main__':
    runner.generate(S2())
