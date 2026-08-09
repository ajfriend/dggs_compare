# /// script
# requires-python = ">=3.11,<3.14"
# dependencies = [
#     "h3>=4",
#     "dggs-compare",
# ]
#
# [tool.uv.sources]
# dggs-compare = { path = "../..", editable = true }
# ///
"""H3 — Uber's hexagonal hierarchical DGGS (gnomonic icosahedral).

Already on the sphere: h3 treats lat/lng as spherical coordinates, and its
edges are geodesics between the boundary vertices (which include the
icosahedron-edge distortion vertices), so higher sampling density is slerp
between them.

Run with:  uv run scripts/systems/h3-h3.py
"""

import h3

from dggs_compare import runner, stats


class H3:
    grid = 'h3'
    impl = 'h3'
    packages = ('h3',)

    def resolutions(self):
        return range(16)            # h3 supports 0..15

    def num_cells(self, res):
        return h3.get_num_cells(res)

    def cells_at(self, res, points):
        return [h3.latlng_to_cell(lat, lng, res) for lat, lng in points]

    def cid_strs(self, cells):
        return [str(c) for c in cells]

    def boundaries(self, res, cells, samples_per_edge=0):
        if samples_per_edge:
            return [stats.refine_geodesic(h3.cell_to_boundary(c),
                                          samples_per_edge)
                    for c in cells]
        return [h3.cell_to_boundary(c) for c in cells]   # (lat, lng) degrees

    def enumerate_cells(self, res):
        for c0 in h3.get_res0_cells():
            if res == 0:
                yield c0
            else:
                yield from h3.cell_to_children(c0, res)


if __name__ == '__main__':
    runner.generate(H3())
