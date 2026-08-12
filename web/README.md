# DGGS shape & area site

A static page comparing the grids by per-cell shape (aspect ratio) and
area: the survey's matplotlib plots (cross-system + by-resolution
distributions, the shape-vs-area tradeoff, best/worst cells) plus a grid
of [ajglobe](https://github.com/ajfriend/ajglobe) globes, one per system,
with cells colored by a selectable metric (AR, or area relative to the
resolution's mean) on a shared scale. No build step, no framework, no
server — just `index.html` + `style.css` + `globe.js` (ES module) and the
vendored ajglobe bundle under `vendor/` (refresh with `just web-vendor`).

```sh
just site   # build everything into web/out/ from the tables (survey PNGs + globes + manifest)
just web    # build, then serve at http://localhost:8000
```

## How it works

Everything the page shows is generated from the Parquet tables — nothing is
hand-authored or pre-baked into git. `just site` runs two steps into `web/out/`
(gitignored):

- **`scripts/survey.py`** → `histograms.png`, `area_histograms.png`,
  `area_ratio_by_res.png`, `tradeoff.png`, `extremes.png`,
  `by_res_<sys>.png` — the distribution, tradeoff, and best/worst-cell
  plots (matplotlib).
- **`scripts/web_data.py`** (via `dggs_compare.webdata`) → `manifest.json` and
  `globe/{sys}_r{res}_{pos.f32,idx.u32,ar.f32,area.f32,ids.json}` — ajglobe's
  flat-binary polygons plus one f32 per cell per metric, one globe per system.

**Globe resolution is area-matched.** Each system's globe uses the resolution
whose cell count is closest (in log-ratio) to H3's at `config.GLOBE_H3_RES`
(default r3). These are ~equal-area grids, so `avg cell area = 4πR²/N(res)` —
matching cell counts matches cell size, computed from the closed-form counts in
`config.CELLS_PER_RES` (no table reads). Raise `GLOBE_H3_RES` for finer/heavier
globes, lower for coarser/lighter ones; changing it needs a full site rebuild.

## Viewer features (`globe.js`)

- **Tabs** — the page is four hash-routed tabs (`#globes` default,
  `#shape`, `#area`, `#tradeoff`): linkable, reload-stable, and plain
  in-page anchors switch tabs. The globes build lazily and resumably on
  showings of their tab (canvases size themselves from the DOM at
  construction, so building only happens while the panel is visible);
  a plots-tab visit costs no WebGL and no globe binaries.
- **Globes** — one per system, cells colored by the selected metric on a
  shared domain so systems compare directly. Drag to rotate, scroll to zoom;
  **all globes are synced** (moving one moves all). Hover a cell for its
  id + value.
- **Metric dropdown** — aspect ratio, or cell area relative to the
  resolution's exact mean; switching re-domains the legend, histograms, and
  per-globe stat lines together.
- **Color-scale dropdowns** — pick a `(colormap, value-transform)`: viridis
  (linear / γ0.4 / p99 / log), cividis, magma, turbo, grayscale. The globes and
  the legend recolor together. The value axis stays linear; only the color
  mapping changes (so, e.g., linear-viridis vs. the stretched default is
  directly comparable — the stretched one is the historical default that keeps
  the low-AR bulk legible).
- **Hovered-globe histogram** — above the color bar, the hovered globe's
  distribution of the active metric (log-count bars, shared value axis), with
  a line marking the hovered cell — aligned to the color bar below.
- **Full-screen plots** — click any survey plot to view it full-screen; there's
  an "open full-size" link to the original PNG. Esc / backdrop / × to close.

## Publishing

Deployed to GitHub Pages from a data release, via `pages.yml`:

```sh
gh workflow run pages.yml -f tag=data-vN                    # full: rebuild from the release tables
gh workflow run pages.yml -f tag=data-vN -f rebuild=false   # fast: front-end-only, reuse cached build
```

The **full** run fetches the release's tables (`just fetch-data`), runs
`just site` so every plot and globe comes from those tables, publishes, and
caches the built `web/out/` (keyed per tag). The manifest records the release
tag (shown in the footer). A **`rebuild=false`** run restores that cache and
only re-copies the static files — no fetch, no rebuild — so an HTML/CSS/JS
change deploys in ~30s. (Run a full deploy first after a new data release to
seed the cache, and whenever you change `GLOBE_H3_RES` or the data.)
