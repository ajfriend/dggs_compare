"""Print the release notes for a data release to stdout — the provenance
block that travels with the published artifacts (fuller per-implementation
facts, including binding versions, ride in each table's Parquet metadata).
Used by `just data-publish`."""

from datetime import date
from importlib.metadata import version

from dggs_compare import checks, config

if __name__ == '__main__':
    impls = ', '.join(f'{g}-{i}' for g, i in checks.implementations())
    versions = ', '.join(
        f'{pkg} {version(pkg)}' for pkg in ('csar', 'sparea'))
    print(f"""\
Data release — the per-cell Parquet tables (one table per
grid-implementation x resolution: geometry + `ar` + `area` columns) and
the web viewer's derived files (flat-named: `globe--*`, `full--*`,
`histograms.json`, `manifest.json`).

- date: {date.today().isoformat()}
- implementations: {impls}
- sampling: seed {config.SEED:#x}, per-res-seed {config.PER_RES_SEED}, \
budget {config.N_CELLS:,}
- solver: gap_tol {config.GAP_TOL:g}, method {config.CSAR_METHOD}
- solver versions: {versions} (binding versions ride in each table's
  Parquet metadata)

Fetch the tables: `just fetch-data <tag>` (or read any asset directly with
pyarrow/DuckDB).""")
