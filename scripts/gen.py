"""Build every (system, resolution) Parquet table: geometry + per-cell stats
(skar aspect ratio, sparea area) in one pass, one native env.

Systems come from the library registry (the systems/ folder); config
(SEED, budgets, target resolutions, solver settings) from
dggs_compare.config. Output -> data/cells/ (gitignored; published as GitHub
data releases).

Run with:  just gen
No CLI args (project convention).
"""

from dggs_compare import cache

if __name__ == '__main__':
    cache.build_all()
