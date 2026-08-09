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


def cid_strs(cells):
    return [format(c, '016x') for c in cells]


def _boundary(c):
    cell = s2sphere.Cell(s2sphere.CellId(c))
    verts = []
    for i in range(4):
        ll = s2sphere.LatLng.from_point(cell.get_vertex(i))
        verts.append((ll.lat().degrees, ll.lng().degrees))
    return verts


def boundaries(res, cells, samples_per_edge=0):
    # s2 cell edges lie in planes through the origin (constant-u/v lines on
    # the cube), i.e. great circles — higher sampling density is slerp
    # between the vertices. Here that's definitional (s2sphere's own
    # containment uses these great circles), so the convergence check
    # reduces to solver numerics: agreement is a convexity theorem, not an
    # empirical finding.
    if samples_per_edge:
        return [stats.refine_geodesic(_boundary(c), samples_per_edge)
                for c in cells]
    return [_boundary(c) for c in cells]


def enumerate_cells(res):
    for cid in s2sphere.CellId.walk(res):
        yield cid.id()
