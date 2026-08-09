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
"""rHEALPix — equal-area HEALPix-derived quad DGGS (via DGGAL).

Run with:  uv run scripts/systems/rhealpix-dggal.py
"""

from dggs_compare.dggal_engine import GridImplAdapter
from dggs_compare import runner

if __name__ == '__main__':
    runner.generate(GridImplAdapter('rhealpix', 'rHEALPix'))
