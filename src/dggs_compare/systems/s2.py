"""S2 — Google's quadtree-on-a-cube DGGS."""

import s2sphere

from dggs_compare import stats


def resolutions():
    return range(31)            # s2sphere supports levels 0..30


def num_cells(res):
    return 6 * 4 ** res


def _cell_at(res, lat, lng):
    # from_lat_lng yields a level-30 leaf; walk up to the cell at `res`.
    leaf = s2sphere.CellId.from_lat_lng(s2sphere.LatLng.from_degrees(lat, lng))
    return leaf.parent(res).id()   # int, hashable


def cells_at(res, points):
    return [_cell_at(res, lat, lng) for lat, lng in points]


def cid_strs(zones):
    return [format(z, '016x') for z in zones]


def _boundary(z):
    cell = s2sphere.Cell(s2sphere.CellId(z))
    ring = []
    for i in range(4):
        ll = s2sphere.LatLng.from_point(cell.get_vertex(i))
        ring.append((ll.lat().degrees, ll.lng().degrees))
    return ring


def boundaries(res, zones):
    return [_boundary(z) for z in zones]


def refined_boundaries(res, zones, refine):
    # s2 cell edges lie in planes through the origin (constant-u/v lines on
    # the cube), i.e. great circles — refinement is slerp between corners.
    # Here that's definitional (s2sphere's own containment uses these great
    # circles), so validate-corners reduces to a solver-numerics check:
    # agreement is a convexity theorem, not an empirical finding.
    return [stats.refine_geodesic(_boundary(z), refine) for z in zones]


def enumerate_cells(res):
    for cid in s2sphere.CellId.walk(res):
        yield cid.id()
