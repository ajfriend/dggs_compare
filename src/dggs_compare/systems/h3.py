"""H3 — Uber's hexagonal hierarchical DGGS (gnomonic icosahedral)."""

import h3

from dggs_compare import stats


def resolutions():
    return range(16)            # h3 supports 0..15


def num_cells(res):
    return h3.get_num_cells(res)


def cells_at(res, points):
    return [h3.latlng_to_cell(lat, lng, res) for lat, lng in points]


def cid_strs(zones):
    return [str(z) for z in zones]


def boundaries(res, zones):
    return [h3.cell_to_boundary(z) for z in zones]   # (lat, lng) deg corners


def refined_boundaries(res, zones, refine):
    # h3 edges ARE geodesics between the boundary vertices (which already
    # include icosahedron-edge distortion vertices), so refinement is slerp.
    return [stats.refine_geodesic(h3.cell_to_boundary(z), refine)
            for z in zones]


def enumerate_cells(res):
    for c0 in h3.get_res0_cells():
        if res == 0:
            yield c0
        else:
            yield from h3.cell_to_children(c0, res)
