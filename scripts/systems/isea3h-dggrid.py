# /// script
# requires-python = ">=3.11,<3.14"
# dependencies = [
#     "dggs-compare",
# ]
#
# [tool.uv.sources]
# dggs-compare = { path = "../..", editable = true }
# ///
"""ISEA3H — Snyder equal-area icosahedral aperture-3 hex (via DGGRID).

The second isea3h implementation (issue #57): independent geometry to
check against isea3h-dggal and to compare runtimes; PRIMARY_IMPL keeps
consumers on the dggal tables. Same batch-subprocess engine as isea4t;
DGGRID computes on the authalic sphere, so its coordinates are used as-is
(no authalic transform — that is the dggal script's record, for a binding
that declares WGS84 geodetic output).

MAX_RES 33 matches the dggal implementation's range so the tables align
resolution-for-resolution (DGGRID itself accepts finer; 10*3^33+2 seqnums
sit well inside uint64).

The odd-level icosa-edge kink (issue #25) is a property of the GRID, not
of dggal: DGGRID also emits 6 density-0 vertices for those cells, so the
same convergence carve-out (config.CONV_EXPECTED_RED, keyed by grid)
covers this implementation. Whether its residuals match dggal's is one of
the things the cross-implementation comparison is for.

Run with:  uv run scripts/systems/isea3h-dggrid.py
"""

from _dggrid_engine import Engine
from dggs_compare import runner

MAX_RES = 33


class ISEA3H:
    grid = 'isea3h'
    impl = 'dggrid'
    packages = ()

    def __init__(self):
        self._engine = Engine('ISEA3H')
        # Zero-padded to the max seqnum's width: rows are sorted by cid
        # TEXT for spatially coherent Parquet pages, and variable-width
        # decimals would sort '10' < '9'.
        self._cid_width = len(str(self.num_cells(MAX_RES)))

    def resolutions(self):
        return range(MAX_RES + 1)

    def num_cells(self, res):
        return 10 * 3 ** res + 2

    def cells_at(self, res, points):
        return self._engine.cells_at(res, points)

    def cid_strs(self, cells):
        return [f'{c:0{self._cid_width}d}' for c in cells]

    def boundaries(self, res, cells, samples_per_edge=0):
        return self._engine.boundaries(res, cells, samples_per_edge)

    def enumerate_cells(self, res):
        yield from self._engine.enumerate_ids(res)


if __name__ == '__main__':
    runner.generate(ISEA3H())
