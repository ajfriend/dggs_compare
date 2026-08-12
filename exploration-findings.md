# Exploration findings (archived)

Conclusions from one-shot investigations whose scripts have been
removed from the repo (they lived in `scripts/explorations/` and
`notebooks/`; retrieve with `git log --diff-filter=D --name-only` if
ever needed). Recorded here because the findings existed only in those
files' docstrings. Numbers were measured 2026-08 against the data-v10
tables ("working resolution" = that release's TARGET_RES calibration);
this file is an archive and is not re-verified against later releases.

## csar AR vs vertex-PCA AR (metric validity)

Per-cell comparison on the ISEA7H working-resolution tables: bulk
corr(csar enclosing-cone AR, PCA second-moment AR) = **1.0000**,
median |diff| = 0.0000 — the published AR fields are the real, smooth
Tissot distortion field, not an artifact of choosing the enclosing
metric. The near-1.0 tail peels off the diagonal for a knowable
reason: isolated icosahedral-seam cells are markedly *irregular*
hexagons (edge ratio ~1.25, vertex-radius ratio ~1.32) whose bounding
ellipse happens to be near-circular — csar reports ~1.0 faithfully
(the enclosing shape *is* round) while PCA stays ~1.1–1.16. These are
grid accidents, not round cells and not metric error. Genuinely
isotropic cells exist only at the ~20 face centroids, measure ~zero
under uniform sampling. (Relevant to the methodology discussion in
issue #42.)

## The ISEA7H "dark spots" (and csar faithfulness)

Sharp low-AR cells in ISEA7H are *razor-thin* dips embedded in the
high-distortion seams — AR ~1.0 at a point, ~1.34 just 0.05° away —
not smooth face-centroid lows. Known specimen: the r10 cell containing
lat/lng **(-71.90959, 140.97260)**. At that adversarial cell, csar's
AR matches an independent Khachiyan minimum-area-enclosing-ellipse
computation on the gnomonic-projected vertices to ~5 decimals
(1.00191 vs 1.00191): not a solver bug. The MVEE cross-check code and
this cell were proposed to csar_py as standing validation material
(filed upstream).

## Projection distortion signatures (ISEA vs IVEA)

Cell AR is resolution-independent for small cells, so per-cell AR maps
picture each projection's distortion field. The global view is on the
site's AR-colored globes; the archived part is the finer face zoom
(gnomonic, r10): ISEA is a low core with a sharp 6-pointed seam star;
IVEA a smooth, bounded radial pinwheel. Generator: `ar_heatmap.py` in
git history (needs a live dggal engine env).

## Finest-resolution AR spot checks (from the retired notebook)

Sampled at each system's finest resolution (H3 r15, S2 L30, A5 r30),
AR distributions match the working-resolution survey — AR is
~resolution-stable — with A5 tightly clustered around 2:1: its cells
are never circular. (Per-resolution detail is better served by the
published tables and the site's by-resolution plots.)
