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
"""ISEA3H — Snyder equal-area icosahedral aperture-3 hex (via DGGAL).

KNOWN APPROXIMATION (issue #25): stats use density-0 vertices, like every
other grid, but isea3h's odd-level cells straddling an icosahedron edge
kink there — the real boundary bulges past the density-0 hexagon by 100s
of meters (ground-truthed via point->cell), so density-0 AR is off by up
to ~4e-3 for those cells, and the convergence gate carves this grid out
(config.CONV_EXPECTED_RED). Revisit per issue #25.

Run with:  uv run scripts/systems/isea3h-dggal.py
"""

from _dggal_engine import Adapter
from dggs_compare import runner, stats

if __name__ == '__main__':
    runner.generate(Adapter('ISEA3H', to_sphere=stats.authalic_rings,
                            pentagons=True))
