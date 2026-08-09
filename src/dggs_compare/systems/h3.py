"""H3 — Uber's hexagonal hierarchical DGGS (gnomonic icosahedral)."""

import h3

from dggs_compare import stats


def resolutions():
    return range(16)            # h3 supports 0..15


def num_cells(res):
    return h3.get_num_cells(res)


def cells_at(res, points):
    return [h3.latlng_to_cell(lat, lng, res) for lat, lng in points]


def cid_strs(cells):
    return [str(c) for c in cells]


def boundaries(res, cells, samples_per_edge=0):
    # h3 edges ARE geodesics between the boundary vertices (which already
    # include icosahedron-edge distortion vertices), so higher sampling
    # density is slerp between them. Note: since denser vertices come from
    # the same edge model, the convergence check reduces to solver numerics
    # here (agreement is a convexity theorem) — it cannot falsify the
    # geodesic-edge model itself; only a point-classification ground-truth
    # test could.
    if samples_per_edge:
        return [stats.refine_geodesic(h3.cell_to_boundary(c), samples_per_edge)
                for c in cells]
    return [h3.cell_to_boundary(c) for c in cells]   # (lat, lng) degrees


def enumerate_cells(res):
    for c0 in h3.get_res0_cells():
        if res == 0:
            yield c0
        else:
            yield from h3.cell_to_children(c0, res)
