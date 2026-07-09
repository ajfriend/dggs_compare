# DGGS aspect-ratio web explorer

An interactive view of the comparison: dynamic, overlaid histograms for any
system/resolution, and two synced orthographic globes with cells colored by
aspect ratio. A static page — no build step, no server framework. The globes
are rendered by [ajglobe](https://github.com/ajfriend/ajglobe) (vendored under
`vendor/`; refresh with `just web-vendor`).

```sh
just web        # build the data (if needed) + serve at http://localhost:8000
just web-data   # just (re)build the data
```

## How it works

`scripts/web_data.py` (via `dggs_compare.webdata`) reshapes the Parquet
tables' columns into browser-friendly data under `out/` (gitignored) —
nothing is solved here; the stats were computed at table-generation time:

- `histograms.json` — for every `(system, resolution)`, a **fixed** fine-bin
  histogram plus summary stats; the fixed grid lets the page re-aggregate to
  any coarser bin width in the browser.
- `globe/{sys}_r{res}_{pos.f32,idx.u32,ar.f32,ids.json}` — ajglobe's native
  flat-binary polygon format for the coarse resolutions (open rings, any
  winding — ajglobe triangulates by ring topology, so there's no orientation
  or antimeridian preprocessing).
- `manifest.json` — what exists, per-system colors/labels, the solve
  tolerance, and the shared globe AR max.

The page (`index.html` + `app.js`, shared helpers in `_shared.js`; d3 +
Observable Plot from a CDN for the histograms, vendored ajglobe for the
globes):

- **Histograms** (Observable Plot) — pick any `system · resolution` series;
  one shared AR axis, a bins slider, density/count, linear/log.
- **Globes** (ajglobe `Orb`, WebGL2) — two GPU-rendered globes that rotate
  and zoom together, GPU hover picking (cell id + AR), country outlines,
  per-globe/shared color domains; toggling restyles in place.

## Full-globe experiment

`globe_full.html` renders **every** ivea7h cell (up to 1.18M at r6) rather
than the sampled tables — same binary format, produced by `just web-full`
(geometry enumeration + native skar solve; both passes skip levels already
on disk) into `out/full/`.
