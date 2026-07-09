_:
    just --list

# Generate random-cell Parquet sets for each DGGS (dggs_cache/cells/). Each
# gen_*.py is a standalone PEP 723 / uv-run script carrying its own DGGS
# library + Python, so the libraries never share an env; output goes to
# dggs_cache/cells/out/ (gitignored, cached). Re-run only for fresh/larger
# sets. dggal ships an arch-broken arm64 wheel; gen_dggal self-re-execs under
# x86_64/Rosetta on Apple Silicon (so a plain `uv run` works everywhere).
gen-cells: gen-h3 gen-s2 gen-a5 gen-dggal

gen-h3:
    uv run dggs_cache/cells/gen_h3.py

gen-s2:
    uv run dggs_cache/cells/gen_s2.py

gen-a5:
    uv run dggs_cache/cells/gen_a5.py

gen-dggal:
    uv run dggs_cache/cells/gen_dggal.py

# Run the DGGS aspect-ratio survey at an H3-r9-matched resolution. Reads the
# pre-generated Parquet cell sets (run `just gen-cells` first), solves each
# with skar, writes PNGs to dggs_cache/out/. DGGS-library-free, native.
survey:
    uv run --group cells dggs_cache/survey.py

# Recalibrate the resolutions that match H3 r9 cell area. Reads the small cell
# sets (`just gen-cells` first); skar-free and DGGS-library-free (standalone
# uv script). Bake the picks into the generators' TARGET.
calibrate:
    uv run dggs_cache/calibrate.py

# Check the DNC invariants across every DGGS / resolution: working resolutions
# clean, and DNC only at the finest sub-metre levels (monotone, no islands).
# Exits non-zero on a regression. Doubles as skar_py's pre-release gate: run
# it (after bumping the skar tag in pyproject.toml + `uv sync`) against a
# release candidate before tagging skar.
dnc-check:
    uv run --group cells dggs_cache/dnc_check.py

# Build the static data for the DGGS aspect-ratio web viewer: solves every
# cached cell with skar and emits out/histograms.json + out/globe/* (ajglobe's
# flat binaries) + out/manifest.json under dggs_cache/web/ (gitignored).
web-data:
    uv run --group cells dggs_cache/web/build_data.py

# Serve the DGGS aspect-ratio web viewer at http://localhost:8000 (builds the
# data first). Static page: dynamic histograms (any system/resolution) + two
# synced orthographic globes (ajglobe) with cells colored by aspect ratio.
web: web-data
    uv run -m http.server 8000 -d dggs_cache/web

# Refresh the vendored ajglobe bundle (checked in under web/vendor/) from the
# sibling repo's dist (`just build` over there first if src changed).
web-vendor:
    cp ../ajglobe/dist/ajglobe.min.js dggs_cache/web/vendor/ajglobe.min.js

# Full-globe experiment (globe_full.html): EVERY ivea7h cell at r1-r3 + the
# r5/r6 torture tests (not the sampled cache). Two passes. (1) gen the
# complete geometry to binary via DGGAL under Rosetta; (2) solve aspect ratios
# natively with skar. ~minutes for r6 (1.18M cells); both passes skip levels
# already on disk. Output -> dggs_cache/web/out/full/ (gitignored).
web-full-geom:
    uv run dggs_cache/web/gen_ivea7h_full_geom.py

# Run `just web-full-geom` first (the slow Rosetta pass); this solves the AR.
web-full:
    uv run --group cells dggs_cache/web/build_ivea7h_full_ar.py

# The dggs_old live-engine analyses need dggal, whose macOS arm64 wheel is
# arch-broken — build the dedicated x86_64 (Rosetta) env for them:
#   just dggs-sync
#   UV_PROJECT_ENVIRONMENT=.venv-dggs uv run --no-sync dggs_old/explorations/<script>.py
dggs_python := "cpython-3.13-macos-x86_64"
dggs_env := ".venv-dggs"

dggs-sync:
    UV_PROJECT_ENVIRONMENT={{dggs_env}} uv sync --python {{dggs_python}} --group dggs

# Open JupyterLab (notebooks/).
lab:
    uv run --group lab jupyter lab

purge:
    just _rm .venv
    just _rm .venv-dggs
    just _rm .DS_Store
    just _rm __pycache__
    just _rm uv.lock
    just _rm .ipynb_checkpoints

_rm pattern:
    -@find . -name "{{pattern}}" -prune -exec rm -rf {} +
