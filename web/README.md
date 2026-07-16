# DGGS aspect-ratio site

A static page comparing the grids by per-cell aspect ratio: the survey's
matplotlib plots (cross-system + by-resolution distributions, best/worst
cells) plus a grid of [ajglobe](https://github.com/ajfriend/ajglobe) globes,
one per system, with cells colored by AR on a shared scale. No build step, no
framework, no server — just `index.html` + `style.css` + `globe.js` (ES
module) and the vendored ajglobe bundle under `vendor/` (refresh with
`just web-vendor`).

```sh
just site   # build everything into web/out/ from the tables (survey PNGs + globes + manifest)
just web    # build, then serve at http://localhost:8000
```

## How it works

Everything the page shows is generated from the Parquet tables — nothing is
hand-authored or pre-baked into git. `just site` runs two steps into `web/out/`
(gitignored):

- **`scripts/survey.py`** → `histograms.png`, `extremes.png`,
  `by_res_<sys>.png` — the distribution and best/worst-cell plots (matplotlib).
- **`scripts/web_data.py`** (via `dggs_compare.webdata`) → `manifest.json` and
  `globe/{sys}_r{res}_{pos.f32,idx.u32,ar.f32,ids.json}` — ajglobe's flat-binary
  polygons for each system's coarsest full-coverage resolution.

`globe.js` reads the manifest, draws one globe per system (rotate/zoom, hover
for a cell's id + AR), and colors every cell by AR on the shared `[1, max]`
scale in the manifest so the systems compare directly.

## Publishing

The hosted site is deployed to GitHub Pages from a data release, built from
that release's tables at deploy time:

```sh
gh workflow run pages.yml -f tag=data-vN
```

`pages.yml` fetches the release's tables (`just fetch-data`), runs `just site`
so every plot and globe comes from those tables, and publishes. The manifest
records the release tag (shown in the page footer).
