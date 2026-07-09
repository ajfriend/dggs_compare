"""Build the (system, resolution) Parquet tables: geometry + per-cell stats
(skar aspect ratio, sparea area) in one pass, one env.

Systems come from the library registry (the systems/ folder); config
(SEED, budgets, target resolutions, solver settings) from
dggs_compare.config. Output -> data/cells/ (gitignored; published as GitHub
data releases).

Run with:  just gen              (every system)
           just gen isea7h       (one system — how CI parallelizes, one
                                  runner per system, artifacts merged after)

The selector arrives via the DGGS_COMPARE_GEN env var (set by the justfile
recipe) rather than CLI args, per project convention.
"""

import os

from dggs_compare import cache

if __name__ == '__main__':
    which = os.environ.get('DGGS_COMPARE_GEN', 'all')
    if which == 'all':
        cache.build_all()
    else:
        cache.build_system(which)
