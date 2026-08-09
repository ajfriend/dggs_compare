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
"""IVEA3H — icosahedral vertex-oriented equal-area aperture-3 hex (via
DGGAL). Same layout as ISEA3H but with distortion spread smoothly instead
of concentrated in seams (the same ISEA/IVEA relationship as the 7H pair).

Run with:  uv run scripts/systems/ivea3h-dggal.py
"""

from _dggal_engine import Adapter
from dggs_compare import runner, stats

if __name__ == '__main__':
    runner.generate(Adapter('IVEA3H', to_sphere=stats.authalic_rings))
