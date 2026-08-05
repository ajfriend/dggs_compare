// Static DGGS aspect-ratio site: a grid of ajglobe globes, one per system,
// cells colored by aspect ratio on a SHARED scale so systems compare directly.
// - a Color-scale dropdown picks (colormap, value-transform); globes + legend
//   recolor together.
// - hovering a globe shows THAT globe's AR distribution as a histogram above
//   the legend, with a line marking the hovered cell's AR (aligned to the
//   color bar, so you see where that cell sits in both distribution and color).
// - all globes rotate/zoom together.
import { Orb } from './vendor/ajglobe.min.js';

// ---- colormaps: RGB control points (0–255), linearly interpolated ----
const CMAPS = {
  viridis: [[68,1,84],[72,40,120],[62,74,137],[49,104,142],[38,130,142],
            [31,158,137],[53,183,121],[110,206,88],[181,222,43],[253,231,37]],
  cividis: [[0,34,78],[0,55,110],[62,73,106],[112,94,100],[150,116,94],
            [190,140,86],[230,169,73],[255,201,58],[255,233,69]],
  magma:   [[0,0,4],[28,16,68],[79,18,123],[129,37,129],[181,54,122],
            [229,80,100],[251,135,97],[254,194,135],[252,253,191]],
  turbo:   [[48,18,59],[54,88,196],[36,133,237],[30,183,208],[45,224,152],
            [122,247,79],[189,238,52],[245,197,45],[251,131,42],[226,72,28],
            [175,31,17],[122,4,3]],
  gray:    [[20,20,20],[245,245,245]],
};
const clamp01 = (x) => (x < 0 ? 0 : x > 1 ? 1 : x);
function interp(stops, t) {
  t = clamp01(t) * (stops.length - 1);
  const i = Math.min(stops.length - 2, t | 0), f = t - i, a = stops[i], b = stops[i + 1];
  return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f, 255];
}
const LUTS = {};
const lutFor = (name) => (LUTS[name] ??= Array.from({ length: 256 }, (_, i) => interp(CMAPS[name], i / 255)));

// ---- value transforms: aspect ratio (>= 1) -> t in [0,1] for the COLOR.
// The AR axis (histogram + legend position) stays linear over [1, max]; only
// the color placement changes, so a stretched scale shows saturation visibly.
const TF = {
  linear: (ar, d) => (ar - 1) / (d.max - 1),
  power:  (ar, d) => Math.pow((ar - 1) / (d.max - 1), 0.4),   // the "modified" stretch
  p99:    (ar, d) => clamp01((ar - 1) / (d.p99 - 1)),         // linear, saturates at p99
  log:    (ar, d) => (d.max > 1 ? Math.log(ar) / Math.log(d.max) : 0),
};
const SCALES = [
  { label: 'Viridis · γ0.4 (modified — the default)', cmap: 'viridis', tf: 'power' },
  { label: 'Viridis · linear (perceptually uniform in AR)', cmap: 'viridis', tf: 'linear' },
  { label: 'Viridis · linear, saturates at p99', cmap: 'viridis', tf: 'p99' },
  { label: 'Viridis · log', cmap: 'viridis', tf: 'log' },
  { label: 'Cividis · linear (uniform, colorblind-optimized)', cmap: 'cividis', tf: 'linear' },
  { label: 'Magma · γ0.4', cmap: 'magma', tf: 'power' },
  { label: 'Turbo · linear (rainbow — NOT perceptually uniform)', cmap: 'turbo', tf: 'linear' },
  { label: 'Grayscale · linear', cmap: 'gray', tf: 'linear' },
];

const DNC_GREY = [68, 68, 68, 255];
const HBINS = 64;
const fmt = (x, d = 4) => (x == null || !Number.isFinite(x) ? '—' : x.toFixed(d));
const $ = (s) => document.querySelector(s);

async function fetchBin(path, Ctor) {
  const r = await fetch(path);
  if (!r.ok) throw new Error(`${path} ${r.status}`);
  return new Ctor(await r.arrayBuffer());
}

const DOMAIN = { max: 1, p99: 1 };   // shared AR domain over all globe cells
const PANELS = [];                   // { sys, label, color, orb, layer, ar, hist }
let scale = SCALES[0];

const makeFill = (ar, sc) => {
  const lut = lutFor(sc.cmap), tf = TF[sc.tf];
  return (i) => {
    const a = ar[i];
    if (!Number.isFinite(a)) return DNC_GREY;
    return lut[Math.min(255, Math.max(0, (tf(a, DOMAIN) * 255) | 0))];
  };
};

function drawLegend(sc) {
  $('#legLo').textContent = '1.0';
  $('#legHi').textContent = fmt(DOMAIN.max, 2);
  const lut = lutFor(sc.cmap), tf = TF[sc.tf], n = 96;
  const stops = Array.from({ length: n }, (_, i) => {
    const p = i / (n - 1), ar = 1 + p * (DOMAIN.max - 1);
    const c = lut[Math.min(255, (tf(ar, DOMAIN) * 255) | 0)];
    return `rgb(${c[0] | 0},${c[1] | 0},${c[2] | 0}) ${100 * p}%`;
  }).join(',');
  $('#legendGrad').style.background = `linear-gradient(90deg, ${stops})`;
}

function applyScale(sc) {
  scale = sc;
  drawLegend(sc);
  for (const p of PANELS) p.layer.update({ fill: makeFill(p.ar, sc) });
}

// ---- shared AR domain + per-globe histograms (bars are scale-independent) ----
function computeDomainAndHists() {
  const all = [];
  for (const p of PANELS) for (const a of p.ar) if (Number.isFinite(a)) all.push(a);
  all.sort((x, y) => x - y);
  DOMAIN.max = all.length ? all[all.length - 1] : 1;
  DOMAIN.p99 = all.length ? all[Math.floor(0.99 * (all.length - 1))] : DOMAIN.max;
  const span = DOMAIN.max - 1 || 1;
  for (const p of PANELS) {
    const counts = new Array(HBINS).fill(0);
    for (const a of p.ar) {
      if (!Number.isFinite(a)) continue;
      counts[Math.min(HBINS - 1, Math.max(0, Math.floor((a - 1) / span * HBINS)))]++;
    }
    p.hist = { counts, cmax: Math.max(1, ...counts) };
  }
}

// ---- the dynamic histogram canvas (shows the hovered globe) ----
let hctx, hcw = 460, hch = 84;
function initHist() {
  const c = $('#histCanvas');
  hcw = c.clientWidth || 460;
  const dpr = Math.min(2, window.devicePixelRatio || 1);
  c.width = hcw * dpr; c.height = hch * dpr;
  c.style.height = hch + 'px';
  hctx = c.getContext('2d');
  hctx.scale(dpr, dpr);
}
function drawHist(panel, hlAR) {
  $('#histLabel').innerHTML = panel
    ? `<span class="swatch" style="background:${panel.color}"></span>${panel.label} — ${panel.ar.length.toLocaleString()} cells`
    : 'hover a globe to see its distribution';
  hctx.clearRect(0, 0, hcw, hch);
  if (!panel) return;
  const { counts, cmax } = panel.hist, bw = hcw / HBINS, logmax = Math.log(cmax + 1);
  hctx.fillStyle = panel.color;
  for (let i = 0; i < HBINS; i++) {
    if (!counts[i]) continue;
    const h = (Math.log(counts[i] + 1) / logmax) * (hch - 3);
    hctx.fillRect(i * bw, hch - h, Math.max(1, bw - 0.6), h);
  }
  if (hlAR != null && Number.isFinite(hlAR)) {
    const x = clamp01((hlAR - 1) / (DOMAIN.max - 1)) * hcw;
    hctx.strokeStyle = '#fff'; hctx.lineWidth = 1.5;
    hctx.beginPath(); hctx.moveTo(x, 0); hctx.lineTo(x, hch); hctx.stroke();
    hctx.fillStyle = '#fff'; hctx.font = '11px ui-monospace, Menlo, monospace';
    const t = `AR ${hlAR.toFixed(3)}`, tw = hctx.measureText(t).width;
    hctx.fillText(t, x + tw + 6 > hcw ? x - 4 - tw : x + 4, 11);
  }
}

// ---- rotate/zoom sync across all globes ----
let syncing = false;
function syncFrom(orb) {
  if (syncing) return;
  syncing = true;
  const v = orb.getView();
  for (const p of PANELS) if (p.orb !== orb) p.orb.setView(v);
  syncing = false;
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
    fig.append(cap, img);   // caption above the plot
    grid.appendChild(fig);
  }
}

async function buildGlobe(M, sys) {
  // globe_res is now the single area-matched resolution per system; tolerate
  // the older list form (a cached manifest from a previous build).
  const gr = M.globe_res[sys];
  const res = Array.isArray(gr) ? gr.at(-1) : gr;
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
  orb.borders({ color: '#2a3340', width: 1 });

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

  const layer = orb.polygons({ lnglat: pos, starts, fill: makeFill(ar, scale) });
  const panel = { sys, label: `${M.labels[sys]} ${M.res_prefix[sys]}${res}`, color: M.colors[sys], orb, layer, ar };
  PANELS.push(panel);

  orb.on('viewchange', () => syncFrom(orb));

  const tip = $('#tooltip');
  orb.on('hover', (e) => {
    orb.highlight(e.index ?? -1, layer);
    const a = e.index == null ? null : ar[e.index];
    drawHist(panel, Number.isFinite(a) ? a : null);   // this globe's distribution + line
    if (e.index == null) { tip.style.opacity = 0; return; }
    tip.innerHTML = `${ids[e.index]}<br>AR ${Number.isFinite(a) ? fmt(a, 4) : 'DNC'}`;
    const r = canvas.getBoundingClientRect();
    tip.style.left = `${r.left + e.x + 14}px`;
    tip.style.top = `${r.top + e.y + 14}px`;
    tip.style.opacity = 1;
  });
  canvas.addEventListener('pointerleave', () => { tip.style.opacity = 0; drawHist(panel, null); });
}

// Click any plot (static histogram/extremes or a dynamically-added by-res
// image) to view it full-screen; click/Esc to close, or open the original.
function initLightbox() {
  const box = $('#lightbox'), img = $('#lightboxImg'), raw = $('#lightboxRaw');
  const open = (src, alt) => {
    img.src = src; img.alt = alt || ''; raw.href = src;
    box.classList.add('open'); document.body.style.overflow = 'hidden';
  };
  const close = () => {
    box.classList.remove('open'); img.removeAttribute('src'); document.body.style.overflow = '';
  };
  document.addEventListener('click', (e) => {
    const t = e.target;
    if (t.matches('.plot, .byres-grid img')) open(t.currentSrc || t.src, t.alt);
    else if (t === box || t.id === 'lightboxImg' || t.id === 'lightboxClose') close();
  });
  document.addEventListener('keydown', (e) => {
    if (e.key === 'Escape' && box.classList.contains('open')) close();
  });
}

async function main() {
  const M = await fetch('out/manifest.json').then((r) => r.json());
  $('#subtitle').textContent =
    `${M.systems.length} systems · aspect ratio via csar's enclosing-cone solver`
    + ` · gap_tol ${M.gap_tol.toExponential()}`;
  if (M.tag) $('#tag').textContent = M.tag;

  byResGrid(M);
  initLightbox();
  initHist();

  const sel = $('#scaleSelect');
  SCALES.forEach((s, i) => sel.add(new Option(s.label, i)));
  sel.value = '0';

  DOMAIN.max = M.globe_ar_max || 1;   // provisional, so the first fill is sane
  DOMAIN.p99 = DOMAIN.max;
  for (const sys of M.systems) await buildGlobe(M, sys);   // one at a time; 6 WebGL panels
  computeDomainAndHists();
  applyScale(SCALES[0]);
  drawHist(PANELS[0], null);          // seed the histogram before any hover

  sel.addEventListener('change', () => applyScale(SCALES[+sel.value]));
}
main();
