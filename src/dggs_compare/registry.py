"""System discovery: the `systems/` folder IS the registry.

`systems/` contains exactly one module per DGGS and nothing else (it's an
implicit namespace package — no __init__.py), so adding a grid to the whole
pipeline is dropping in one file that implements the contract below. `get()`
imports a system's module on first use — the ONLY laziness in the pipeline —
so `import dggs_compare` and the table-reading consumers never load a DGGS
binding (or spin up the DGGAL engine); inside a system module, its binding
is a plain top-level import.

Every systems/ module implements:

    resolutions()            -> range of generatable resolutions (0..finest)
    num_cells(res)           -> int, total cells at `res`
    cell_at(res, lat, lng)   -> native cell id (hashable) at the point, or
                                None when the engine can't resolve it (rare
                                DGGAL deep-level singular points) — samplers
                                skip and draw again
    cid_str(z)               -> str text id stored in the tables
    cell_boundary(z)         -> [(lat, lng), ...] degrees, open corner ring

DGGAL-backed modules additionally expose `adapter()` — the live-engine
`dggal_engine.Adapter` (neighbors, edge-refined vertices, …) for analyses
that need more than a bag of cells.

Optional: `stats_ring(z)` -> [(lat, lng), ...] open ring fed to the AR/area
solvers INSTEAD of cell_boundary(z), for systems whose corner ring does not
faithfully bound the cell (isea3h: odd-level cells kink at icosahedron
edges). The tables' verts column always stores the corner ring.
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
