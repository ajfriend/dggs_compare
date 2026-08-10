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

Density-0 stats are exact here (spot-checked to ~1e-12 against densely
sampled boundaries): the aperture-4 triangle lattice keeps icosahedron
edges ON cell edges at every level, so there are no distortion vertices.

Run with:  uv run scripts/systems/isea4t-dggrid.py
"""

from _dggrid_engine import Adapter
from dggs_compare import runner

if __name__ == '__main__':
    runner.generate(Adapter('ISEA4T', max_res=28,
                            num_cells=lambda r: 20 * 4 ** r))
