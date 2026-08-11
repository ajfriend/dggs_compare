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
consumers on the dggal tables. DGGRID computes on the authalic sphere,
so its coordinates are used as-is.

DGGAL_ORIENTATION (see _dggrid_engine, its one home) pins dggal's
documented icosahedron placement, so the two implementations' grids
coincide cell-for-cell (measured: max vertex discrepancy ~20 cm, the
AIGEN 12-decimal text floor; per-cell AR and area agree to ~1e-10) and
cross-implementation comparison is direct rather than up-to-a-rotation.

max_res 33 matches the dggal implementation's range so the tables align
resolution-for-resolution (DGGRID itself accepts finer; 10*3^33+2 seqnums
sit well inside uint64).

The odd-level icosa-edge kink (issue #25) is a property of the GRID, not
of dggal: DGGRID also emits 6 density-0 vertices for those cells, so the
same convergence carve-out (config.CONV_EXPECTED_RED, keyed by grid)
covers this implementation.

Run with:  uv run scripts/systems/isea3h-dggrid.py
"""

from _dggrid_engine import DGGAL_ORIENTATION, Adapter
from dggs_compare import runner

if __name__ == '__main__':
    runner.generate(Adapter('ISEA3H', max_res=33,
                            num_cells=lambda r: 10 * 3 ** r + 2,
                            params=DGGAL_ORIENTATION, pentagons=True))
