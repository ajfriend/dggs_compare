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
"""IVEA7H — icosahedral vertex-oriented equal-area aperture-7 hex (via
DGGAL). Same layout as ISEA7H but with distortion spread smoothly instead
of concentrated in seams.

Run with:  uv run scripts/systems/ivea7h-dggal.py
"""

from dggs_compare.dggal_engine import GridImplAdapter
from dggs_compare import runner

if __name__ == '__main__':
    runner.generate(GridImplAdapter('ivea7h', 'IVEA7H'))
