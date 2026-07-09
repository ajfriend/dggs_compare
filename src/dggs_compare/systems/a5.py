"""A5 — pentagonal equal-area DGGS (via the Rust/PyO3 a5_fast binding)."""

import a5_fast as a5


def resolutions():
    return range(31)            # a5 supports 0..30 (a5.MAX_RESOLUTION)


def num_cells(res):
    return a5.get_num_cells(res)


def cell_at(res, lat, lng):
    return a5.lonlat_to_cell(lng, lat, res)   # int, hashable


def cid_str(z):
    return a5.u64_to_hex(z)


def cell_boundary(z):
    ring = a5.cell_to_boundary(z)   # closed ring of (lng, lat)
    return [(lat, lng) for lng, lat in ring]   # -> (lat, lng); closing repeat handled downstream


def enumerate_cells(res):
    for c0 in a5.get_res0_cells():
        if res == 0:
            yield c0
        else:
            yield from a5.cell_to_children(c0, res)
