"""dggs_compare — internal library for the DGGS comparison pipeline.

Code organization only: this package is consumed exclusively by the scripts
and workflows in this repo and is never published to PyPI. The repo's actual
products are the data artifacts (GitHub data releases of the per-cell Parquet
tables) and the web pages built from them.

Modules:
    config        pipeline constants — the single source of truth
    registry      system discovery (the systems/ folder IS the registry)
    systems/      one module per DGGS, nothing else
    dggal_engine  shared DGGAL glue + the live-engine Adapter
    stats         per-cell aspect ratio (skar) + area (sparea)
    cache         the Parquet tables: build (geometry + stats) and read
    checks        DNC invariants + corners-only validation
    webdata       web-viewer artifacts derived from the tables

The registry imports a system's module on first use (the only laziness in
the pipeline), so `import dggs_compare` and the table-reading consumers stay
light — they never load a DGGS binding.
"""

from . import config, registry  # noqa: F401  (light, dependency-free)
