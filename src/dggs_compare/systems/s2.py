"""S2 — Google's quadtree-on-a-cube DGGS."""

import s2sphere


def resolutions():
    return range(31)            # s2sphere supports levels 0..30


def num_cells(res):
    return 6 * 4 ** res


def cell_at(res, lat, lng):
    # from_lat_lng yields a level-30 leaf; walk up to the cell at `res`.
    leaf = s2sphere.CellId.from_lat_lng(s2sphere.LatLng.from_degrees(lat, lng))
    return leaf.parent(res).id()   # int, hashable


def cid_str(z):
    return format(z, '016x')


def cell_boundary(z):
    cell = s2sphere.Cell(s2sphere.CellId(z))
    ring = []
    for i in range(4):
        ll = s2sphere.LatLng.from_point(cell.get_vertex(i))
        ring.append((ll.lat().degrees, ll.lng().degrees))
    return ring


def enumerate_cells(res):
    for cid in s2sphere.CellId.walk(res):
        yield cid.id()
