"""A5 — pentagonal equal-area DGGS (via the Rust/PyO3 a5_fast binding).

a5_fast's `cell_to_boundary` returns an adaptively densified ring (321
points at res 0 down to 6 from res ~9 up). The extra points don't change the
enclosing-cone AR — validated to max |dAR| = 3.5e-9 across resolutions by
scripts/validate_corners_a5.py — so `cell_boundary` reduces the ring to its
corner vertices (turning-angle peaks): 5 for the pentagons, 3 for the res-1
"quintant" triangles.
"""

import a5_fast as a5
import numpy as np

# Exterior-angle threshold marking a corner (degrees). Densified edge points
# turn by small fractions of a degree; true corners by tens of degrees.
_TURN_DEG = 5.0


def resolutions():
    return range(31)            # a5 supports 0..30 (a5.MAX_RESOLUTION)


def num_cells(res):
    return a5.get_num_cells(res)


def cells_at(res, points):
    return [a5.lonlat_to_cell(lng, lat, res)   # int, hashable
            for lat, lng in points]


def cid_str(z):
    return a5.u64_to_hex(z)


def _corners(latlng):
    """Reduce an open (lat, lng)-degree ring to its corner vertices."""
    la, lo = np.radians(np.asarray(latlng, dtype=float)).T
    v = np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo),
                         np.sin(la)])
    e = np.roll(v, -1, axis=0) - v                    # chord to next vertex
    e /= np.linalg.norm(e, axis=1, keepdims=True)
    cosang = np.clip(np.einsum('ij,ij->i', np.roll(e, 1, axis=0), e), -1, 1)
    idx = np.nonzero(np.degrees(np.arccos(cosang)) > _TURN_DEG)[0]
    return [latlng[i] for i in idx] if len(idx) >= 3 else latlng


def _boundary(z):
    ring = a5.cell_to_boundary(z)   # closed ring of (lng, lat), densified
    return _corners([(lat, lng) for lng, lat in ring[:-1]])


def boundaries(zones):
    return [_boundary(z) for z in zones]


def enumerate_cells(res):
    for c0 in a5.get_res0_cells():
        if res == 0:
            yield c0
        else:
            yield from a5.cell_to_children(c0, res)
