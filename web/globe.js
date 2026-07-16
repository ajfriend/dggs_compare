// Static DGGS aspect-ratio site: a grid of ajglobe globes, one per system,
// cells colored by aspect ratio on a SHARED scale (so systems compare directly).
// Data comes from the tables via `just site` (out/manifest.json + out/globe/*).
// No dynamic histograms, no synced panels — the distribution plots are the
// pre-rendered matplotlib PNGs on the page; this file only draws the globes.
import { Orb } from './vendor/ajglobe.min.js';

// ---- value-based color: viridis over [1, shared max] with a gamma stretch ----
// (equal AR ⇒ equal color; the sub-linear stretch keeps the low-AR bulk legible)
const STOPS = [[68, 1, 84], [71, 44, 122], [59, 81, 139], [44, 113, 142],
               [33, 144, 141], [39, 173, 129], [92, 200, 99], [253, 231, 37]];
function viridis(t) {
  t = Math.max(0, Math.min(1, t)) * (STOPS.length - 1);
  const i = Math.min(STOPS.length - 2, t | 0), f = t - i, a = STOPS[i], b = STOPS[i + 1];
  return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f, 255];
}
const LUT = Array.from({ length: 256 }, (_, i) => viridis(i / 255));
const lut = (t) => LUT[Math.min(255, (t * 255) | 0)];
const DNC_GREY = [68, 68, 68, 255];
const GAMMA = 0.4;

const fmt = (x, d = 4) => (x == null || !Number.isFinite(x) ? '—' : x.toFixed(d));
const $ = (s) => document.querySelector(s);

async function fetchBin(path, Ctor) {
  const r = await fetch(path);
  if (!r.ok) throw new Error(`${path} ${r.status}`);
  return new Ctor(await r.arrayBuffer());
}

function legend(hi) {
  $('#legLo').textContent = '1.0';
  $('#legHi').textContent = fmt(hi, 2);
  const n = 48;
  const stops = Array.from({ length: n }, (_, i) => {
    const f = i / (n - 1), c = viridis(Math.pow(f, GAMMA));
    return `rgb(${c[0] | 0},${c[1] | 0},${c[2] | 0}) ${100 * f}%`;
  }).join(',');
  $('#legendGrad').style.background = `linear-gradient(90deg, ${stops})`;
}

function byResGrid(M) {
  const grid = $('#byresGrid');
  for (const sys of M.systems) {
    const fig = document.createElement('figure');
    const img = document.createElement('img');
    img.src = `out/by_res_${sys}.png`;
    img.alt = `${M.labels[sys]} aspect-ratio distribution by resolution`;
    img.loading = 'lazy';
    const cap = document.createElement('figcaption');
    cap.innerHTML = `<span class="swatch" style="background:${M.colors[sys]}"></span>${M.labels[sys]}`;
    fig.append(img, cap);
    grid.appendChild(fig);
  }
}

async function drawGlobe(M, sys) {
  const res = M.globe_res[sys].at(-1);   // most-detailed full-coverage resolution
  const hi = M.globe_ar_max;

  const card = document.createElement('div');
  card.className = 'globe-card';
  card.innerHTML =
    `<div class="globe-label"><span class="swatch" style="background:${M.colors[sys]}"></span>`
    + `${M.labels[sys]} <span class="res">${M.res_prefix[sys]}${res}</span></div>`
    + `<div class="globe-holder"><canvas></canvas></div><div class="globe-meta"></div>`;
  $('#globeGrid').appendChild(card);

  const canvas = card.querySelector('canvas');
  const orb = new Orb(canvas, { background: '#0b0e13', sphere: '#11151c' });
  orb.lookAt(0, 20);
  orb.borders({ color: '#2a3340', width: 1 });   // Natural Earth outlines (CDN)

  const key = `${sys}_r${res}`;
  const [pos, starts, ar, ids] = await Promise.all([
    fetchBin(`out/globe/${key}_pos.f32`, Float32Array),
    fetchBin(`out/globe/${key}_idx.u32`, Uint32Array),
    fetchBin(`out/globe/${key}_ar.f32`, Float32Array),
    fetch(`out/globe/${key}_ids.json`).then((r) => r.json()),
  ]);

  let amax = 1;
  for (const a of ar) if (Number.isFinite(a) && a > amax) amax = a;
  card.querySelector('.globe-meta').textContent =
    `${ids.length.toLocaleString()} cells · max AR ${fmt(amax, 3)}`;

  const fill = (i) => {
    const a = ar[i];
    return Number.isFinite(a) ? lut(Math.pow((Math.min(a, hi) - 1) / (hi - 1), GAMMA)) : DNC_GREY;
  };
  const layer = orb.polygons({ lnglat: pos, starts, fill });

  // Light hover: highlight the picked cell + a tooltip with its id and AR.
  const tip = $('#tooltip');
  orb.on('hover', (e) => {
    orb.highlight(e.index ?? -1, layer);
    if (e.index == null) { tip.style.opacity = 0; return; }
    const a = ar[e.index];
    tip.innerHTML = `${ids[e.index]}<br>AR ${Number.isFinite(a) ? fmt(a, 4) : 'DNC'}`;
    const r = canvas.getBoundingClientRect();
    tip.style.left = `${r.left + e.x + 14}px`;
    tip.style.top = `${r.top + e.y + 14}px`;
    tip.style.opacity = 1;
  });
  canvas.addEventListener('pointerleave', () => { tip.style.opacity = 0; });
}

async function main() {
  const M = await fetch('out/manifest.json').then((r) => r.json());
  $('#subtitle').textContent =
    `${M.systems.length} systems · aspect ratio via csar's enclosing-cone solver`
    + ` · gap_tol ${M.gap_tol.toExponential()}`;
  if (M.tag) $('#tag').textContent = M.tag;

  byResGrid(M);
  legend(M.globe_ar_max);
  for (const sys of M.systems) await drawGlobe(M, sys);   // one at a time; 6 WebGL panels
}
main();
