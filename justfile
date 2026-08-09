_:
    just --list

# The project env holds the library, solvers, and plotting; each implementation
# script in scripts/systems/ resolves its OWN env from its PEP 723 header
# (binding + library, never co-resolved with anything else). Platform
# wrinkle: dggal 0.0.6's macOS arm64 wheel bundles x86_64 dylibs, so on
# Apple Silicon everything runs x86_64 under Rosetta; Linux is native (CI
# generates the canonical data there).
python := if os() == "macos" { "cpython-3.13-macos-x86_64-none" } else { "3.13" }

# Exported so EVERY `uv run` in these recipes (including the per-script
# envs) uses the platform pin: a bare `uv run` otherwise consults
# .python-version, which knows the version but not the macOS ARCH — on
# Apple Silicon it would silently pick native arm64 (dggal's broken
# dylibs). Env var beats the file.
export UV_PYTHON := python

sync:
    uv sync --python {{python}}

# The registry, read as a key list ({grid}-{impl} per line): the CI build
# matrix and `gen all` both consume this — the one shell reading of the
# scripts/systems/ listing. Underscore-prefixed files are shared engine
# modules, not registry entries.
systems:
    @basename -a -s .py scripts/systems/[!_]*.py

# Stage 1: raw cell geometry -> data/raw/ (gitignored). One PEP 723 script
# per (grid, implementation) in scripts/systems/, each in its own env; the
# folder listing is the registry and the CI matrix. Pass a key to run one
# (e.g. `just gen isea3h-dggal`).
gen key="all":
    #!/usr/bin/env sh
    if [ "{{key}}" != "all" ]; then
        exec uv run "scripts/systems/{{key}}.py"
    fi
    # Run every implementation even if one fails (mirrors CI's
    # fail-fast: false): one run yields the full failure set AND all the
    # healthy raw artifacts, which `just metrics` can then use per-key.
    fail=0
    for k in $(just systems); do
        uv run "scripts/systems/$k.py" || { echo "FAILED: $k"; fail=1; }
    done
    exit $fail

# Stage 2: raw geometry -> the published tables in data/cells/ (per-cell
# csar AR + sparea area, one solver provenance for every system), applying
# the convergence admission gate.
metrics:
    uv run scripts/metrics.py

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

# Report the convergence residuals stamped into the published tables (the
# gate itself runs in `just metrics`). A pure metadata read.
convergence:
    uv run scripts/convergence.py

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
# `gh release download <tag> -p '<grid>-<impl>_r<res>.parquet'` instead if that's
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
