# dggs_compare

Comparing discrete global grid systems — **H3, S2, A5, ISEA7H, IVEA7H,
rHEALPix** — by per-cell statistics: shape (aspect ratio via
[skar](https://github.com/ajfriend/skar_py)'s enclosing-cone solver) and area
(via [sparea](https://pypi.org/project/sparea/)), at area-matched resolutions.

Extracted from [skar_py](https://github.com/ajfriend/skar_py) (see its history
through PR #16 for how this pipeline evolved); split out so the library repo
stays a library and the investigation can grow freely.

**Aspect ratio (AR)** throughout means the major/minor semi-axis ratio (a/b,
a≥b) of a cell's enclosing-cone ellipse — the discrete per-cell analogue of
Tissot's indicatrix. AR = 1 is isotropic; AR > 1 is anisotropic ("squished").
Since these grids are ~equal-area, AR is where the unavoidable distortion
surfaces.

## Layout

```
dggs_cache/     the Parquet-cache pipeline: generate random cells per DGGS
                once, then run every analysis off those files — natively,
                with no DGGS library in the analysis env. Survey plots,
                resolution calibration, DNC invariants, explorations, and an
                interactive web viewer (ajglobe globes + dynamic histograms).
dggs_old/       live-engine analyses that can't run off a cell snapshot
                (point->cell scans, neighbors, edge refinement); dggal under
                a dedicated x86_64/Rosetta env on Apple Silicon.
notebooks/      interactive companions (dggs_survey.ipynb).
```

## Use

```sh
just gen-cells   # generate the Parquet cell sets (once; ~2 min, dggal under Rosetta)
just survey      # aspect-ratio survey -> dggs_cache/out/*.png
just calibrate   # area-match resolutions across systems
just dnc-check   # assert the DNC invariants (pass/fail)
just web         # interactive explorer (histograms + synced globes) at :8000
```

See `dggs_cache/README.md` and `dggs_old/README.md` for the full tour.

## Relationship to skar_py

skar comes from the released tag pinned in `pyproject.toml`
(`[tool.uv.sources]`); bump it when skar_py cuts a release. This repo also
serves as **skar's pre-release regression gate**: before tagging a skar
release, point the pin at the candidate, `uv sync`, and run `just dnc-check`
(and ideally `just survey` for a numeric eyeball) — it exercises ~5.7M real
DGGS cells against invariants that have caught real solver differences.
