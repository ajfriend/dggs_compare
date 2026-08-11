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

Batch subprocess calls instead of per-cell FFI (see _dggrid_engine; the
dggrid binary comes from `just install-dggrid` or DGGS_COMPARE_DGGRID).
Cells are plain SEQNUM ints; they are only unique within a resolution,
which is why the contract passes `res` alongside. Uses DGGRID's default
orientation (vert0 11.25 E) — no dggal counterpart to align with.

max_res 28: DGGRID accepts finer, but 20*4^r overflows uint64 seqnums
past r30 (observed: r31 ids come out smaller than r30's), and r28 cells
are already ~6 cm^2.

KNOWN APPROXIMATION (issue #53): density-0 AR is exact here
(spot-checked to ~1e-12; the aperture-4 lattice keeps icosahedron edges
ON cell edges at every level, so no cell straddles a fold), but
density-0 AREA is not — ISEA's angular distortion curves the
projected-straight edges, so the 3-corner chord is off by a
scale-invariant ±16% for the icosa-edge-adjacent class at every level,
and broadly at coarse levels. The true cells are exactly equal-area;
densified boundaries restore that (measurements in issue #53).

Run with:  uv run scripts/systems/isea4t-dggrid.py
"""

from _dggrid_engine import Adapter
from dggs_compare import runner

if __name__ == '__main__':
    runner.generate(Adapter('ISEA4T', max_res=28,
                            num_cells=lambda r: 20 * 4 ** r))
