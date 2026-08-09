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
made of. Cell handles are plain hashables whose meaning may be scoped to
a resolution (DGGRID seqnums), so every geometry call carries `res`.

A cell's boundary is always the same kind of object — a list of vertices,
in order, joined by great-circle edges — and the only knob is the SAMPLING
DENSITY: how many extra vertices are sampled along each native edge.
Density 0 is the minimal vertex set. For systems whose native edges are
not great circles, higher density traces the true edge more closely;
scripts/convergence.py verifies the metrics are converged in this knob.
Every systems/ module implements:

    resolutions()          -> range of generatable resolutions (0..finest)
    num_cells(res)         -> int, total cells at `res`
    cells_at(res, points)  -> [cell or None, ...] for [(lat, lng), ...];
                              None = the engine couldn't resolve the point
                              (rare DGGAL deep-level singular points) —
                              samplers skip and draw again
    cid_strs(cells)        -> [str, ...] text ids stored in the tables;
                              must sort (as text) in a spatially coherent
                              order — cache.build_table orders rows by
                              them for Parquet page locality (fixed-width)
    boundaries(res, cells, samples_per_edge=0)
                           -> [[(lat, lng), ...] open vertex list, ...]
                              aligned with `cells`, degrees, with AT LEAST
                              `samples_per_edge` extra vertices sampled
                              along each native edge (binding APIs
                              quantize, so modules may return more)
    enumerate_cells(res)   -> iterator of every cell at `res`

DGGAL-backed modules additionally expose `adapter()` — the live-engine
`dggal_engine.Adapter` (neighbors, point->cell, …) for analyses that need
more than a bag of cells.
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
