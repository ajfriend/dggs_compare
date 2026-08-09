# /// script
# requires-python = ">=3.11,<3.14"
# dependencies = [
#     "dggal>=0.0.6",
#     "dggs-compare",
# ]
#
# [tool.uv.sources]
# dggs-compare = { path = "../..", editable = true }
# ///
"""ISEA7H — Snyder equal-area icosahedral aperture-7 hex (via DGGAL).

Run with:  uv run scripts/systems/isea7h-dggal.py
"""

from dggs_compare.dggal_engine import Adapter
from dggs_compare import runner


class Impl:
    grid = 'isea7h'
    impl = 'dggal'
    packages = ('dggal',)

    def __init__(self):
        self._a = Adapter('ISEA7H')

    def resolutions(self):
        return range(self._a.max_level() + 1)

    def num_cells(self, res):
        return self._a.count(res)

    def cells_at(self, res, points):
        return self._a.cells_at(res, points)

    def cid_strs(self, cells):
        return self._a.cid_strs(cells)

    def boundaries(self, res, cells, samples_per_edge=0):
        return self._a.boundaries(cells, samples_per_edge)

    def enumerate_cells(self, res):
        yield from self._a.enumerate(res)


if __name__ == '__main__':
    runner.generate(Impl())
