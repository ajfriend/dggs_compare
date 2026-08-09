"""dggs_compare — internal library for the DGGS comparison pipeline.

Code organization only: this package is consumed exclusively by the scripts
and workflows in this repo and is never published to PyPI. The repo's actual
products are the data artifacts (GitHub data releases of the per-cell Parquet
tables) and the web pages built from them.

Modules:
    config        pipeline constants — the single source of truth
    interface     the implementation contract (GridImpl)
    runner        stage 1: GridImpl -> raw geometry parquet (data/raw/)
    metrics       stage 2: raw -> published tables (binding-free)
    dggal_engine  shared DGGAL glue + the live-engine Adapter
    dggrid_engine DGGRID batch-subprocess engine
    stats         solvers (csar AR, sparea area) + sphere helpers
    cache         the published tables: IO + readers (data/cells/)
    checks        DNC invariants + artifact/config coherence
    webdata       web-viewer artifacts derived from the tables

The implementation registry is the scripts/systems/ file listing; each
script resolves its own env. `import dggs_compare` stays light — the
library core never loads a DGGS binding or a solver.
"""

from . import config  # noqa: F401  (light, dependency-free)
