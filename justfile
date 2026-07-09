_:
    just --list

# ONE env for everything — no per-script envs, no re-exec tricks. The single
# wrinkle: dggal 0.0.6's macOS arm64 wheel is still half-broken (the Python
# extension is arm64 but the bundled libecrt/libdggal dylibs are x86_64), so
# on Apple Silicon the whole env runs x86_64 under Rosetta. Linux wheels are
# correct — CI generates the canonical data natively on ubuntu.
python := if os() == "macos" { "cpython-3.13-macos-x86_64-none" } else { "3.13" }

sync:
    uv sync --python {{python}}

# Build the (system, resolution) Parquet tables: geometry + per-cell stats
# (skar AR, sparea area) in one pass. Systems come from the library registry
# (src/dggs_compare/systems/ — one file per DGGS). Output -> data/cells/
# (gitignored; published as GitHub data releases). Pass a system name to
# build just that one (how CI parallelizes: one runner per system).
gen system="all":
    DGGS_COMPARE_GEN={{system}} uv run scripts/gen.py

# Aspect-ratio survey: reads the ar column (no solving) -> out/histograms.png,
# extremes.png, by_res_<system>.png.
survey:
    uv run scripts/survey.py

# Match each system's resolution to an H3-r9 cell by area (reads the area
# column). Bake picks into dggs_compare.config.TARGET_RES.
calibrate:
    uv run scripts/calibrate.py

# DNC invariants: working resolutions clean, DNC only at the finest sub-metre
# levels, monotone. Exits non-zero on a regression. Set RESOLVE=True in the
# script to re-solve with the installed skar (the skar pre-release gate).
dnc-check:
    uv run scripts/dnc_check.py

# Corners-vs-edge-refined AR validation for the DGGAL grids — run when
# admitting a new DGGAL system to the pipeline.
validate-corners:
    uv run scripts/validate_corners.py

# Build the web viewer's static data (histograms + ajglobe globe binaries +
# manifest) from the tables -> web/out/ (gitignored). A column reshape.
web-data:
    uv run scripts/web_data.py

# Serve the web viewer at http://localhost:8000 (builds the data first).
web: web-data
    uv run -m http.server 8000 -d web

# Refresh the vendored ajglobe bundle (checked in under web/vendor/) from the
# sibling repo's dist (`just build` over there first if src changed).
web-vendor:
    cp ../ajglobe/dist/ajglobe.min.js web/vendor/ajglobe.min.js

# Full-globe experiment (web/globe_full.html): EVERY ivea7h cell at r1-r3 +
# the r5/r6 torture tests. Two passes, both native; each skips levels already
# on disk. Output -> web/out/full/ (gitignored).
web-full-geom:
    uv run scripts/web_full_geom.py

web-full: web-full-geom
    uv run scripts/web_full_ar.py

# ----- data releases ------------------------------------------------------
# The tables + web data are published as decoupled GitHub releases (data-v1,
# data-v2, ...), cut when the INPUTS change (seed/budgets, the grid set,
# calibration, or a solver bump worth reflecting) — not per code release.
# Normally cut from CI: gh workflow run data-release.yml -f tag=data-vN

# Stage everything a data release ships: the tables as one tar + the web
# viewer files flat-named for the release's flat asset namespace
# (out/globe/<f> -> globe--<f>, out/full/<f> -> full--<f>).
data-pack:
    rm -rf data-stage && mkdir -p data-stage
    tar -cf data-stage/cells-parquet.tar -C data cells
    cp web/out/histograms.json web/out/manifest.json data-stage/
    for f in web/out/globe/*; do cp "$f" "data-stage/globe--$(basename "$f")"; done
    for f in web/out/full/*; do cp "$f" "data-stage/full--$(basename "$f")"; done
    ls data-stage | wc -l

# Create the GitHub release and upload the staged assets (run data-pack first).
data-publish tag:
    uv run scripts/data_notes.py > data-stage/NOTES.md
    gh release create {{tag}} --title {{tag}} --notes-file data-stage/NOTES.md
    gh release upload {{tag}} data-stage/cells-parquet.tar data-stage/*.json
    gh release upload {{tag}} data-stage/globe--*
    gh release upload {{tag}} data-stage/full--*

# Download a data release's tables into data/cells/ — the instant alternative
# to `just gen`.
fetch-data tag:
    mkdir -p data
    gh release download {{tag}} -p cells-parquet.tar -O - | tar -xf - -C data

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
