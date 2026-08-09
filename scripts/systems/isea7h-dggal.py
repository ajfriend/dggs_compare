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
from dggs_compare import runner, stats

if __name__ == '__main__':
    runner.generate(Adapter('ISEA7H', to_sphere=stats.authalic_rings))
