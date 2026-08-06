# dggs_compare

Comparing discrete global grid systems — **H3, S2, A5, ISEA7H, IVEA7H,
**ISEA3H, IVEA3H, ISEA4T, rHEALPix, hex9** — by per-cell statistics: shape (aspect
ratio via [csar](https://github.com/ajfriend/csar_py)'s enclosing-cone
solver) and area (via [sparea](https://pypi.org/project/sparea/)), at
area-matched resolutions.

**The products of this repo are data artifacts and web pages.** The Parquet
tables — one row per cell, with the computed `ar` and `area` bundled next to
the geometry — are published as GitHub data releases; downstream consumers
(the analyses here, the web viewer, or anyone with pyarrow/DuckDB) just read
columns. The code is organized as an internal library (`src/dggs_compare/`)
with thin scripts on top; it is never published to PyPI.

**Aspect ratio (AR)** throughout means the major/minor semi-axis ratio (a/b,
a≥b) of a cell's enclosing-cone ellipse — the discrete per-cell analogue of
Tissot's indicatrix. AR = 1 is isotropic; AR > 1 is anisotropic ("squished").
Since these grids are ~equal-area, AR is where the unavoidable distortion
surfaces.

Not listed separately: **IGEO7** is geometrically identical to ISEA7H (it is
ISEA7H with Z7 cell indexing), so its shape/area statistics are already
covered by the ISEA7H tables — only the cell ID strings would differ.
**IVEA4T** has no implementation to compare against: DGGRID has no IVEA
projection, and dggal has no triangle grids.

## Layout

```
src/dggs_compare/     the internal library (organization only, not for PyPI)
  config.py             pipeline constants — the single source of truth
  systems/              ONE FILE PER DGGS, nothing else; the folder is the
                        registry (see "Adding a grid")
  registry.py           folder discovery + lazy imports
  dggal_engine.py       shared DGGAL glue + the live-engine Adapter
  dggrid_engine.py      DGGRID batch-subprocess engine (backs isea4t)
  stats.py              per-cell AR (csar) + area (sparea)
  cache.py              the Parquet tables: build + read (data/cells/)
  checks.py             DNC invariants (cached-ar or re-solve modes)
  webdata.py            web-viewer artifacts from the tables
scripts/              thin callers: gen, survey, calibrate, dnc_check,
                      validate_corners, web_*, explorations/
web/                  the static comparison site (survey plots + ajglobe
                      globes colored by AR); data from web/out/
data/cells/           the tables (gitignored; published as data releases)
notebooks/            interactive companions
```

Everything runs in **one env** (`just sync`) — no per-script envs, no re-exec
tricks; the registry imports a system's module on first use, so table-reading
consumers never load a DGGS binding. One platform wrinkle: dggal 0.0.6's
macOS arm64 wheel still bundles x86_64 dylibs, so on Apple Silicon the env
runs x86_64 under Rosetta (one justfile line); Linux is native (CI generates
the canonical data there). The one system outside the env: isea4t shells out
to a DGGRID binary (no wheels exist anywhere) — `just install-dggrid` builds
it into `.tools/` in ~2 min; CI installs it only on the runner that needs it.

## Use

```sh
just gen               # build all tables: geometry + stats, one pass (~min)
just survey            # AR comparison plots -> out/
just calibrate         # area-match resolutions across systems
just dnc-check         # assert the DNC invariants (pass/fail)
just site              # build the static site into web/out/ (survey plots + globes)
just web               # build the site, then serve it at :8000
just validate-corners  # stats-inputs-vs-refined-edges check, every system
```

## Adding a grid

The `systems/` folder is the registry — the data-release gen matrix is
derived from it at run time. A new DGGS touches:

1. `src/dggs_compare/systems/<name>.py` — the module (the batch-first
   contract is `registry.py`'s docstring; mirror any existing system; add
   `stats_rings` only if the corner rings are not faithful).
2. `config.py` — every dict in `PER_SYSTEM`: `CELLS_PER_RES`, `TARGET_RES`
   (count-match, then confirm with `just calibrate` once tables exist),
   `SYS_COLOR`.

Then run `just validate-corners` (the admission gate) and cut a data release — `just dnc-check` fails its publish gate
unless every registry system has tables and `PER_SYSTEM` config.

## Table schema

One Parquet file per `(system, resolution)` at `data/cells/{sys}_r{res}.parquet`:

| column | type | |
|---|---|---|
| `dggs` | string | system name (constant per file) |
| `res` | int32 | resolution/level (constant per file) |
| `cid` | string | cell id text |
| `verts` | list<[lat, lng] f64> | corner ring, degrees, open |
| `ar` | float64 | enclosing-cone aspect ratio; NaN = did-not-certify |
| `area` | float64 | spherical area, steradians |

Provenance (seed policy, budgets, solver settings, library versions) rides in
each file's Parquet metadata. Coarse resolutions are enumerated in full;
finer ones hold exactly 1,000,000 sampled cells (enumerate-and-subsample
near the cap, draw-until-n beyond it — see `config.N_CELLS`).

## Data releases

The tables + the viewer's derived files are published as **decoupled GitHub
releases** (`data-v1`, `data-v2`, …) — cut when the *inputs* change
(seed/budgets, the grid set, calibration, a solver bump worth reflecting),
not per code release. The canonical producer is CI (native linux):

```sh
gh workflow run data-release.yml -f tag=data-vN   # generate + gate + publish
just fetch-data data-vN                           # pull the tables (10s of GB)
```

The optional `-f runner=ubuntu-24.04-arm` input switches every job in the
release to arm64 runners (faster single-core; one runner type per release so
a data artifact never mixes float provenance across architectures).

Each release carries provenance notes (the same facts ride in every table's
Parquet metadata). Release assets support HTTP Range but send **no CORS
headers**, so browsers can't fetch them cross-origin — CLI/API consumers
(`fetch-data`, pyarrow, DuckDB) are unaffected. The hosted site therefore lives
on **GitHub Pages**, rebuilt from a release's tables at deploy time (never
entering git):

```sh
gh workflow run pages.yml -f tag=data-vN                    # full: fetch tables, `just site`, publish
gh workflow run pages.yml -f tag=data-vN -f rebuild=false   # fast: front-end-only change, ~30s
```

Every plot and globe on the published site is generated in the full job from the
release's tables — the site is a pure function of the data artifact. A full run
also caches the built `web/out/` (keyed per tag); `rebuild=false` restores that
cache and only re-copies the static files (HTML/CSS/JS), so a viewer-only change
deploys in ~30s instead of re-fetching all the tables. The first run after a
data release must be a full one, to seed the cache.

## Relationship to csar_py

csar comes from the released tag pinned in `pyproject.toml`
(`[tool.uv.sources]`). This repo is also **csar's pre-release regression
gate**: point the pin at a candidate, `uv sync`, set `RESOLVE=True` in
`scripts/dnc_check.py`, and run it — that re-solves ~5.7M real cells with the
installed csar against invariants that have caught real solver differences.

Extracted from [skar_py](https://github.com/ajfriend/skar_py) (history through
its PR #16); restructured for native dggal 0.0.6 (no more Rosetta) with stats
bundled into the tables.
