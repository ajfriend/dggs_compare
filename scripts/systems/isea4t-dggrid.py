# /// script
# requires-python = ">=3.11,<3.14"
# dependencies = [
#     "dggs-compare",
# ]
#
# [tool.uv.sources]
# dggs-compare = { path = "../..", editable = true }
# ///
"""ISEA4T — Snyder equal-area icosahedral aperture-4 triangle (via DGGRID).

Batch subprocess calls instead of per-cell FFI (see dggrid_engine; the
dggrid binary comes from `just install-dggrid` or DGGS_COMPARE_DGGRID —
it is not a python package, so `packages` is empty). Cells are plain
SEQNUM ints; they are only unique within a resolution, which is why the
contract passes `res` alongside.

MAX_RES 28: DGGRID accepts finer, but 20*4^r overflows uint64 seqnums
past r30 (observed: r31 ids come out smaller than r30's), and r28 cells
are already ~6 cm^2.

Density-0 stats are exact here (spot-checked to ~1e-12 against densely
sampled boundaries): the aperture-4 triangle lattice keeps icosahedron
edges ON cell edges at every level, so there are no distortion vertices.

Run with:  uv run scripts/systems/isea4t-dggrid.py
"""

from _dggrid_engine import Engine
from dggs_compare import runner

MAX_RES = 28


class ISEA4T:
    grid = 'isea4t'
    impl = 'dggrid'
    packages = ()

    def __init__(self):
        self._engine = Engine('ISEA4T')
        # Zero-padded to the max seqnum's width: rows are sorted by cid
        # TEXT for spatially coherent Parquet pages, and variable-width
        # decimals would sort '10' < '9'.
        self._cid_width = len(str(self.num_cells(MAX_RES)))

    def resolutions(self):
        return range(MAX_RES + 1)

    def num_cells(self, res):
        return 20 * 4 ** res

    def cells_at(self, res, points):
        return self._engine.cells_at(res, points)

    def cid_strs(self, cells):
        return [f'{c:0{self._cid_width}d}' for c in cells]

    def boundaries(self, res, cells, samples_per_edge=0):
        return self._engine.boundaries(res, cells, samples_per_edge)

    def enumerate_cells(self, res):
        yield from self._engine.enumerate_ids(res)


if __name__ == '__main__':
    runner.generate(ISEA4T())
