"""H3 — Uber's hexagonal hierarchical DGGS (gnomonic icosahedral)."""

import h3


def resolutions():
    return range(16)            # h3 supports 0..15


def num_cells(res):
    return h3.get_num_cells(res)


def cell_at(res, lat, lng):
    return h3.latlng_to_cell(lat, lng, res)


def cid_str(z):
    return str(z)


def cell_boundary(z):
    return h3.cell_to_boundary(z)   # [(lat, lng), ...] deg, corners only


def enumerate_cells(res):
    for c0 in h3.get_res0_cells():
        if res == 0:
            yield c0
        else:
            yield from h3.cell_to_children(c0, res)
