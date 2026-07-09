"""Print the release notes for a data release to stdout — the provenance
block that travels with the published artifacts (the same facts ride in each
table's Parquet metadata). Used by `just data-publish`."""

from datetime import date
from importlib.metadata import version

from dggs_compare import config, registry

if __name__ == '__main__':
    versions = ', '.join(
        f'{pkg} {version(pkg)}'
        for pkg in ('skar', 'sparea', 'h3', 's2sphere', 'a5_fast', 'dggal'))
    print(f"""\
Data release — the per-cell Parquet tables (`cells-parquet.tar`, one table
per system x resolution: geometry + `ar` + `area` columns) and the web
viewer's derived files (flat-named: `globe--*`, `full--*`,
`histograms.json`, `manifest.json`; the viewer loads them via
`?data=<this release's download URL>`).

- date: {date.today().isoformat()}
- systems: {', '.join(registry.names())}
- sampling: seed {config.SEED:#x}, per-res-seed {config.PER_RES_SEED}, \
budget {config.N_CELLS:,}
- solver: gap_tol {config.GAP_TOL:g}, method {config.SKAR_METHOD}
- versions: {versions}

Fetch the tables: `just fetch-data <tag>` (or read any asset directly with
pyarrow/DuckDB).""")
