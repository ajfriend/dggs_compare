// Helpers shared by the two viewer pages (app.js and globe_full.html) so the
// color pipeline can't drift between them. App-level code by design — the
// ajglobe library owns no color scales or data loading.

// Alternate data source: pass ?data=<base-url> to load the viewer's data
// from a flat-named mirror instead of the local out/ directory. NOTE:
// GitHub release assets send no CORS headers, so a release URL only works
// same-origin (the hosted page ships its data alongside via Pages); the
// param is for local/same-origin mirrors. Flat names encode the directory:
// out/globe/<f> -> globe--<f>, out/full/<f> -> full--<f>.
export const DATA_BASE =
  new URLSearchParams(location.search).get('data') ?? '';

export function dataURL(path) {
  if (!DATA_BASE) return path;                      // local out/ layout
  const flat = path.replace(/^out\//, '').replaceAll('/', '--');
  return `${DATA_BASE}/${flat}`;
}

export const DNC_GREY = [68, 68, 68, 255];

// viridis via control stops — no per-cell color-string parsing.
const STOPS = [[68, 1, 84], [71, 44, 122], [59, 81, 139], [44, 113, 142],
               [33, 144, 141], [39, 173, 129], [92, 200, 99], [253, 231, 37]];
export function viridis(t) {
  t = Math.max(0, Math.min(1, t)) * (STOPS.length - 1);
  const i = Math.min(STOPS.length - 2, t | 0), f = t - i, a = STOPS[i], b = STOPS[i + 1];
  return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f, 255];
}

// 256-entry LUT of preallocated RGBA arrays: fill callbacks run once per cell
// (up to 1.18M), so color lookups should be indexed, not allocated.
export const VIRIDIS = Array.from({ length: 256 }, (_, i) => viridis(i / 255));
export const lut = (t) => VIRIDIS[Math.min(255, (t * 255) | 0)];

export async function fetchBin(path, Ctor) {
  const r = await fetch(dataURL(path));
  if (!r.ok) throw new Error(`${path} ${r.status}`);
  return new Ctor(await r.arrayBuffer());
}

// Min/max/span of the finite values, for linear value->color scales
// (value-based, never rank: equal AR ⇒ equal color). A span that is
// near-zero relative to the values is float noise (e.g. a level of congruent
// cells, AR all equal) — return span Infinity so the caller's
// (a - min) / span collapses to one color instead of amplifying the noise.
export function finiteRange(values) {
  let min = Infinity, max = -Infinity;
  for (const v of values) {
    if (Number.isFinite(v)) {
      if (v < min) min = v;
      if (v > max) max = v;
    }
  }
  const span = max - min;
  const noise = 1e-4 * Math.max(Math.abs(min), Math.abs(max), 1);
  return { min, max, span: span > noise ? span : Infinity };
}
