"""System discovery: the `systems/` folder IS the registry.

`systems/` contains exactly one module per DGGS and nothing else (it's an
implicit namespace package — no __init__.py), so adding a grid to the whole
pipeline is dropping in one file that implements the contract below. `get()`
imports a system's module on first use — the ONLY laziness in the pipeline —
so `import dggs_compare` and the table-reading consumers never load a DGGS
binding (or spin up the DGGAL engine); inside a system module, its binding
is a plain top-level import.

The contract is batch-first: subprocess-backed systems (isea4t/DGGRID)
resolve a whole chunk per call, and in-process backends just loop — the
reverse adaptation (batch on top of per-cell) is what stateful hacks are
made of. Every systems/ module implements:

    resolutions()          -> range of generatable resolutions (0..finest)
    num_cells(res)         -> int, total cells at `res`
    cells_at(res, points)  -> [zone or None, ...] for [(lat, lng), ...];
                              zones are opaque hashables; None = the engine
                              couldn't resolve the point (rare DGGAL
                              deep-level singular points) — samplers skip
                              and draw again
    cid_str(z)             -> str text id stored in the tables; must sort
                              (as text) in a spatially coherent order —
                              cache.build_table orders rows by it for
                              Parquet page locality (fixed-width ids)
    boundaries(zones)      -> [[(lat, lng), ...] open corner ring, ...]
                              aligned with `zones`, degrees
    enumerate_cells(res)   -> iterator of every zone at `res`

DGGAL-backed modules additionally expose `adapter()` — the live-engine
`dggal_engine.Adapter` (neighbors, edge-refined vertices, …) for analyses
that need more than a bag of cells.

Optional: `stats_rings(zones)` -> [ring or None, ...] fed to the AR/area
solvers INSTEAD of the corner rings where not None — for systems whose
corners do not always bound the cell (isea3h: odd-level cells kink at
icosahedron edges). The tables' verts column always stores corner rings.
"""

import importlib
from functools import cache
from pathlib import Path

SYSTEMS_DIR = Path(__file__).parent / 'systems'


@cache
def names():
    """Sorted system names — the systems/ folder listing."""
    return sorted(p.stem for p in SYSTEMS_DIR.glob('*.py'))


@cache
def get(name):
    """The system module for `name` (lazy import; cached)."""
    if name not in names():
        raise KeyError(f'unknown DGGS {name!r}; available: {names()}')
    return importlib.import_module(f'dggs_compare.systems.{name}')
