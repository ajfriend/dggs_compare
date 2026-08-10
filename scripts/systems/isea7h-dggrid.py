# /// script
# requires-python = ">=3.11,<3.14"
# dependencies = [
#     "dggs-compare",
# ]
#
# [tool.uv.sources]
# dggs-compare = { path = "../..", editable = true }
# ///
"""ISEA7H — Snyder equal-area icosahedral aperture-7 hex (via DGGRID).

The second isea7h implementation (issue #57): independent geometry to
check against isea7h-dggal and to compare runtimes; PRIMARY_IMPL keeps
consumers on the dggal tables. DGGRID computes on the authalic sphere,
so its coordinates are used as-is.

ORIENTATION pins the same icosahedron placement as isea3h-dggrid —
dggal's documented vertex 0 at authalic arctan(golden-ratio) N, 11.20 E,
azimuth 0.

KNOWN DIVERGENCE (`just cross-impl isea7h`, issue #57): EVEN levels
agree with dggal at the numerics floor (~20 cm), but at ODD levels the
two implementations use different density-0 conventions for cells
crossing icosahedron folds — same cell centers, boundary vertices apart
by km (r1: 20 cells at 11.0 km; r3: max |dAR| 1.2e-2). Each is
self-consistent under its own densification (both pass the convergence
gate at ~1e-9); dggal's odd-level polygons tile the sphere exactly,
DGGRID's over-cover slightly (+2.3e-3 sr at r1). The odd-level kink
of #25, surfacing as a cross-implementation convention difference.

max_res 19 matches the dggal implementation's range so the tables align
resolution-for-resolution (DGGRID itself accepts finer; 10*7^19+2
seqnums sit well inside uint64).

Run with:  uv run scripts/systems/isea7h-dggrid.py
"""

from _dggrid_engine import Adapter
from dggs_compare import runner

ORIENTATION = {
    'dggs_vert0_lon': '11.20',
    'dggs_vert0_lat': '58.282525588538994675786',
    'dggs_vert0_azimuth': '0.0',
}

if __name__ == '__main__':
    runner.generate(Adapter('ISEA7H', max_res=19,
                            num_cells=lambda r: 10 * 7 ** r + 2,
                            params=ORIENTATION))
