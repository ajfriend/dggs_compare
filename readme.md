# dggs_compare

Comparing discrete global grid systems — **H3, S2, A5, ISEA7H, IVEA7H,
ISEA3H, IVEA3H, ISEA4T, rHEALPix, hex9** — by per-cell statistics: shape
(aspect ratio via [csar](https://github.com/ajfriend/csar_py)'s enclosing-cone
solver) and area (via [sparea](https://pypi.org/project/sparea/)), at
area-matched resolutions.

**Live site: <https://ajfriend.com/dggs_compare/>** — interactive globes
and every comparison plot, rebuilt from each data release.

**The products of this repo are data artifacts and web pages.** The Parquet
tables — one row per cell, with the computed `ar` and `area` bundled next to
the geometry — are published as GitHub data releases; downstream consumers
(the analyses here, the web viewer, or anyone with pyarrow/DuckDB) just read
columns. The code is organized as an internal library (`src/dggs_compare/`)
with thin scripts on top; it is never published to PyPI.

**Aspect ratio (AR)** throughout means the major/minor semi-axis ratio (a/b,
a≥b) of a cell's enclosing-cone ellipse — the discrete per-cell analogue of
Tissot's indicatrix. AR = 1 is isotropic; AR > 1 is anisotropic (elongated).
Tiling a sphere forces some elongation somewhere — AR measures each cell's
share of it.

Not listed separately: **IGEO7** is geometrically identical to ISEA7H (it is
ISEA7H with Z7 cell indexing), so its shape/area statistics are already
covered by the ISEA7H tables — only the cell ID strings would differ.
**IVEA4T** has no implementation to compare against: DGGRID has no IVEA
projection, and dggal has no triangle grids.

## Layout

```
src/dggs_compare/     the internal library (organization only, not for PyPI;
                      NO system-specific code — it knows only the contract)
  config.py             pipeline constants — the single source of truth
  interface.py          the implementation contract (GridImpl)
  runner.py             stage 1: GridImpl -> raw geometry (data/raw/)
  metrics.py            stage 2: raw -> published tables, binding-free
  stats.py              measurement code (csar AR, sparea area) + sphere helpers
  cache.py              the published tables: IO + readers (data/cells/)
  checks.py             DNC invariants + artifact/config coherence
  webdata.py            web-viewer artifacts from the tables
scripts/systems/      ONE PEP 723 SCRIPT PER (grid, implementation), named
                      {grid}-{impl}.py; the folder is the registry and the
                      CI matrix (see "Adding a grid"). Shared engine glue
                      lives here too, underscore-prefixed (_dggal_engine.py,
                      _dggrid_engine.py) — scripts import it directly
scripts/              thin callers: metrics, survey, calibrate, dnc_check,
                      convergence, web_*
web/                  the static comparison site (survey plots + ajglobe
                      globes colored by a selectable metric: AR or relative
                      area); data from web/out/
data/raw/             stage-1 geometry (gitignored, intermediate)
data/cells/           the published tables (gitignored; data releases)
```

Conclusions from retired one-shot explorations are archived in
`exploration-findings.md`; their scripts live in git history.

Each implementation script resolves its **own env** from its PEP 723
header — its binding plus the library, never co-resolved with any other
system, so implementations can differ in dependencies, python version,
even OS/arch without conflicting. The project env (`just sync`) holds the
library plus the analysis tools and never depends on a DGGS binding.
Platform wrinkle: dggal 0.0.6's macOS arm64 wheel bundles x86_64 dylibs,
so on Apple Silicon everything runs x86_64 under Rosetta (one justfile
line); Linux is native (CI generates the canonical data there). isea4t
shells out to a DGGRID binary (no wheels exist anywhere) —
`just install-dggrid` builds it into `.tools/` in ~2 min.

## Use

```sh
just gen               # stage 1: raw geometry, one env per system -> data/raw/
just gen hex9-hex9     # ... just one implementation
just metrics           # stage 2: measure + gate -> the tables in data/cells/
just survey            # AR + area comparison plots -> out/
just calibrate         # area-match resolutions across systems (closed-form)
just dnc-check         # assert the DNC invariants (pass/fail)
just convergence       # report the stamped convergence residuals
just cross-impl isea3h # cross-implementation agreement report (from data/raw/)
just site              # build the static site into web/out/ (survey plots + globes)
just web               # build the site, then serve it at :8000
```

## Adding a grid

The `scripts/systems/` listing is the registry — the data-release gen
matrix is derived from it. A new implementation touches:

1. `scripts/systems/{grid}-{impl}.py` — a PEP 723 script declaring its own
   dependencies, defining a `GridImpl` class (the contract is
   `interface.py`'s docstring; mirror any existing script), and handing it
   to `runner.generate`. Bringing cells to the sphere is the
   implementation's job; the script is the record of how. Glue shared
   between scripts lives beside them as underscore-prefixed modules
   (never in the library) — a bare `from _dggal_engine import Adapter`
   works because running a script puts its own directory on `sys.path`.
2. `config.py` — every dict in `PER_SYSTEM`: `CELLS_PER_RES`, `TARGET_RES`
   (`just calibrate` picks it from the closed-form counts — no tables
   needed), `SYS_COLOR`, `PRIMARY_IMPL`, and `EXPECTED_IRREGULAR` — decide
   whether the implementation declares exceptional cells (the optional
   `irregular` contract method; a hex grid's 12 pentagons) and record the
   expected count, which `just dnc-check` asserts against the tables.

A second implementation of an existing grid touches only step 1 — the
grid's `PER_SYSTEM` entries already exist; the one decision is whether
`PRIMARY_IMPL` should point at it.

`just metrics` then applies the convergence admission gate (density-0
vertex lists must be faithful inputs to the AR measurement), and
`just dnc-check` fails its publish gate unless every registry
implementation has tables and every grid has `PER_SYSTEM` config.

## Table schema

One Parquet file per `(grid, implementation, resolution)` at
`data/cells/{grid}-{impl}_r{res}.parquet`:

| column | type | |
|---|---|---|
| `dggs` | string | system name (constant per file) |
| `res` | int32 | resolution/level (constant per file) |
| `cid` | string | cell id text |
| `verts` | list<[lat, lng] f64> | vertex list, degrees, open; on the unit sphere (see below) |
| `ar` | float64 | enclosing-cone aspect ratio; NaN = did-not-certify |
| `area` | float64 | spherical area, steradians |
| `irregular` | bool | the implementation's DECLARED exceptional cells (a hex grid's 12 pentagons); never inferred from geometry |

Every cell is a **spherical polygon**: `verts` are sphere coordinates and
edges are great circles, and all metrics are measured on that object.
Bringing cells to the sphere is each implementation script's job —
systems declared on an ellipsoid map their vertices through the
area-preserving authalic latitude (`stats.authalic_rings`), so a system
that is exactly equal-area on its declared surface measures exactly
equal-area here (the methodology issue: #42).

Provenance (seed policy, budgets, csar settings, library versions) rides in
each file's Parquet metadata. Coarse resolutions are enumerated in full;
finer ones hold exactly 1,000,000 sampled cells (enumerate-and-subsample
near the cap, draw-until-n beyond it — see `config.N_CELLS`). Sample
points are drawn uniformly on the sphere and handed to each
implementation in its own coordinate convention (geodetic for
WGS84-declared systems), so the effective sampling measure deviates from
uniform by at most the authalic map's ~0.1% density distortion — this
shifts only which cells are sampled, never any measured value, and is
why `just cross-impl` compares full enumerations only. Consumers note:
a point-sampled table's ROWS are size-biased (cells appear
~proportionally to their area), so distributions computed from one are
size-biased unless weighted by 1/area — the published plots do exactly
that; per-cell values are unaffected.

## Data releases

The tables + the viewer's derived files are published as **decoupled GitHub
releases** (`data-v1`, `data-v2`, …) — cut when the *inputs* change
(seed/budgets, the grid set, calibration, a measurement-code bump worth reflecting),
not per code release. The canonical producer is CI (native linux):

```sh
gh workflow run data-release.yml -f tag=data-vN   # generate + gate + publish
just fetch-data data-vN                           # pull the tables (10s of GB)
```

The optional `-f runner=ubuntu-24.04-arm` input switches every job in the
release to arm64 runners (faster single-core; one runner type per release so
a data artifact never mixes float provenance across architectures).

Every PR runs this exact workflow as its test (`pr.yml` calls it with a
tiny cell budget and publish off): every per-script env, both engine
paths, the gen→metrics handoff, every gate, and release staging —
stopping just short of the upload.

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
data release must be a full one, to seed the cache — and so must the first
run after a change to what `just site` emits (a stale cache lacks the new
files; the viewer degrades gracefully but incompletely).

## Relationship to csar_py

csar comes from PyPI (`csar>=0.1.1`). This repo is also **csar's
pre-release regression gate**: point a `[tool.uv.sources]` entry at a
csar_py branch or rev, `uv sync`, set `RESOLVE=True` in
`scripts/dnc_check.py`, and run it — that re-solves ~5.7M real cells with the
installed csar against invariants that have caught real solver differences.

Extracted from [skar_py](https://github.com/ajfriend/skar_py) (history through
its PR #16); restructured for dggal 0.0.6 with stats
bundled into the tables.
