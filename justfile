_:
    just --list

# ONE env for everything — no per-script envs, no re-exec tricks. Two
# wrinkles: (1) dggal 0.0.6's macOS arm64 wheel is still half-broken (the
# Python extension is arm64 but the bundled libecrt/libdggal dylibs are
# x86_64), so on Apple Silicon the whole env runs x86_64 under Rosetta;
# Linux wheels are correct — CI generates the canonical data natively on
# ubuntu. (2) The versions differ by platform because of two wheel gaps
# that happen to interlock: hex9's wheels are tagged abi3 but ship a
# cpython-312-ONLY extension (broken on 3.13 — the mac x86_64 path is
# immune since no wheel exists there and the sdist builds correctly),
# while a5_fast ships mac wheels only for cp313. So: linux 3.12 (all
# wheels genuinely work), mac x86_64 3.13 (a5 wheel + hex9 sdist).
python := if os() == "macos" { "cpython-3.13-macos-x86_64-none" } else { "3.12" }

sync:
    uv sync --python {{python}}

# Build the (system, resolution) Parquet tables: geometry + per-cell stats
# (csar AR, sparea area) in one pass. Systems come from the library registry
# (src/dggs_compare/systems/ — one file per DGGS). Output -> data/cells/
# (gitignored; published as GitHub data releases). Pass a system name to
# build just that one (how CI parallelizes: one runner per system).
gen system="all":
    DGGS_COMPARE_GEN={{system}} uv run scripts/gen.py

# Aspect-ratio survey: reads the ar column (no solving) -> out/histograms.png,
# extremes.png, by_res_<system>.png.
survey:
    uv run scripts/survey.py

# Match each system's working resolution to an H3-r9 cell by cell count
# (closed forms, no tables). Bake picks into dggs_compare.config.TARGET_RES.
calibrate:
    uv run scripts/calibrate.py

# DNC invariants: working resolutions clean, DNC only at the finest sub-metre
# levels, monotone. Exits non-zero on a regression. Set RESOLVE=True in the
# script to re-solve with the installed csar (the csar pre-release gate).
dnc-check:
    uv run scripts/dnc_check.py

# Stats-input-vs-refined-boundary AR validation for EVERY registry system —
# the check that admits a new grid to the pipeline.
validate-corners:
    uv run scripts/validate_corners.py

# Build the web viewer's static data (histograms + ajglobe globe binaries +
# the full-globe page's complete-coverage binaries + manifest) from the
# tables -> web/out/ (gitignored). A column reshape — nothing is solved.
web-data:
    uv run scripts/web_data.py

# Build EVERYTHING the static site serves, from the tables, into web/out/:
# the survey PNGs (cross-system + by-resolution histograms, best/worst cells)
# AND the ajglobe globe binaries + manifest. This is the single command the
# published site runs (pages.yml) after fetching a data release, so every plot
# on the site is generated from that release's tables — not from anything on
# disk. Reads the `ar` column; nothing is solved but the two extreme cells the
# survey re-draws per system.
site: survey web-data
    cp out/histograms.png out/extremes.png out/by_res_*.png web/out/

# Serve the static site at http://localhost:8000 (builds it first).
web: site
    uv run -m http.server 8000 -d web

# Refresh the vendored ajglobe bundle (checked in under web/vendor/) from the
# sibling repo's dist (`just build` over there first if src changed).
web-vendor:
    cp ../ajglobe/dist/ajglobe.min.js web/vendor/ajglobe.min.js

# ----- data releases ------------------------------------------------------
# The tables + web data are published as decoupled GitHub releases (data-v1,
# data-v2, ...), cut when the INPUTS change (seed/budgets, the grid set,
# calibration, or a solver bump worth reflecting) — not per code release.
# Normally cut from CI: gh workflow run data-release.yml -f tag=data-vN

# Stage the flat-named globe binaries + manifest for the release's flat asset
# namespace (out/globe/<f> -> globe--<f>). The Parquet tables are NOT staged —
# they upload straight from data/cells/ (their names are already flat, and
# skipping the copy halves the publish job's disk footprint at the 1M-cell
# scale). The site itself is rebuilt from the tables at deploy time (pages.yml
# runs `just site`), so these globe assets are a convenience, not the source.
data-pack:
    rm -rf data-stage && mkdir -p data-stage
    cp web/out/manifest.json data-stage/
    for f in web/out/globe/*; do cp "$f" "data-stage/globe--$(basename "$f")"; done
    ls data-stage | wc -l

# Create the GitHub release (if absent) and upload: one asset per Parquet
# table (each far under the 2GB/asset cap; a single tar would exceed it)
# plus the staged viewer files. --clobber everywhere: idempotent, so a
# failed/retried publish just overwrites.
data-publish tag:
    uv run scripts/data_notes.py > data-stage/NOTES.md
    gh release view {{tag}} > /dev/null 2>&1 || \
        gh release create {{tag}} --title {{tag}} --notes-file data-stage/NOTES.md
    gh release upload {{tag}} data/cells/*.parquet --clobber
    gh release upload {{tag}} data-stage/* --clobber

# Download a data release's tables into data/cells/ — the instant alternative
# to `just gen`. (10s of GB at the 1M-cell budget; grab single tables with
# `gh release download <tag> -p '<sys>_r<res>.parquet'` instead if that's
# all you need.)
fetch-data tag:
    mkdir -p data/cells
    gh release download {{tag}} -p '*.parquet' -D data/cells --clobber

# Open JupyterLab (notebooks/).
lab:
    uv run --group lab jupyter lab

purge:
    just _rm .venv
    just _rm .DS_Store
    just _rm __pycache__
    just _rm uv.lock
    just _rm .ipynb_checkpoints

_rm pattern:
    -@find . -name "{{pattern}}" -prune -exec rm -rf {} +

# Build the DGGRID binary into .tools/dggrid (needed by the isea4t system;
# no wheels exist anywhere — clones + cmake-builds without GDAL, ~2 min)
install-dggrid:
    rm -rf .tools/dggrid-src
    git clone --depth 1 https://github.com/sahrk/dggrid.git .tools/dggrid-src
    cmake -S .tools/dggrid-src -B .tools/dggrid-src/build -DCMAKE_BUILD_TYPE=Release
    cmake --build .tools/dggrid-src/build --target dggrid -j 8
    cp .tools/dggrid-src/build/src/apps/dggrid/dggrid .tools/dggrid
    rm -rf .tools/dggrid-src
