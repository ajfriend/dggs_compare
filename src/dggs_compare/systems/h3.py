"""H3 — Uber's hexagonal hierarchical DGGS (gnomonic icosahedral)."""

import h3


def resolutions():
    return range(16)            # h3 supports 0..15


def num_cells(res):
    return h3.get_num_cells(res)


def cells_at(res, points):
    return [h3.latlng_to_cell(lat, lng, res) for lat, lng in points]


def cid_str(z):
    return str(z)


def boundaries(zones):
    return [h3.cell_to_boundary(z) for z in zones]   # (lat, lng) deg corners


def enumerate_cells(res):
    for c0 in h3.get_res0_cells():
        if res == 0:
            yield c0
        else:
            yield from h3.cell_to_children(c0, res)
