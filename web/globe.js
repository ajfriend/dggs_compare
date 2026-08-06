// Static DGGS survey site: a grid of ajglobe globes, one per system, cells
// colored by a per-cell metric on a SHARED scale so systems compare directly.
// - three dropdowns pick the metric (aspect ratio / area ratio), the colormap,
//   and the value transform independently (any combination is valid); globes +
//   legend recolor together.
// - hovering a globe shows THAT globe's distribution of the current metric as
//   a histogram above the legend, with a line marking the hovered cell (aligned
//   to the color bar, so you see where that cell sits in both distribution and
//   color).
// - all globes rotate/zoom together.
import { Orb } from './vendor/ajglobe.min.js';

// ---- the three dropdown axes. One entry per option, in menu order; the FIRST
// entry of each list is the default (labeled so in the UI at build time).
// Adding a metric, colormap or transform is one new object here — the
// dropdowns and lookups all derive from these lists.
// Colormap stops: RGB control points (0–255), linearly interpolated.
const CMAPS = [
  { key: 'viridis', label: 'Viridis (perceptually uniform)',
    stops: [[68,1,84],[72,40,120],[62,74,137],[49,104,142],[38,130,142],
            [31,158,137],[53,183,121],[110,206,88],[181,222,43],[253,231,37]] },
  { key: 'cividis', label: 'Cividis (uniform, colorblind-optimized)',
    stops: [[0,34,78],[0,55,110],[62,73,106],[112,94,100],[150,116,94],
            [190,140,86],[230,169,73],[255,201,58],[255,233,69]] },
  { key: 'magma', label: 'Magma (uniform)',
    stops: [[0,0,4],[28,16,68],[79,18,123],[129,37,129],[181,54,122],
            [229,80,100],[251,135,97],[254,194,135],[252,253,191]] },
  { key: 'cetr2', label: 'CET-R2 (rainbow, perceptually uniform)',
    // Kovesi's uniform rainbow; stops from colorcet.
    stops: [[0,51,245],[0,78,214],[0,99,183],[0,115,154],[40,127,127],
            [57,138,101],[63,151,69],[64,162,33],[87,171,14],[116,177,19],
            [142,184,24],[167,190,29],[191,195,34],[215,201,39],[238,205,43],
            [251,198,42],[253,184,36],[254,168,30],[255,151,23],[255,135,17],
            [255,118,10],[255,99,4],[254,77,0],[252,48,0]] },
  { key: 'turbo', label: 'Turbo (rainbow — NOT perceptually uniform)',
    stops: [[48,18,59],[54,88,196],[36,133,237],[30,183,208],[45,224,152],
            [122,247,79],[189,238,52],[245,197,45],[251,131,42],[226,72,28],
            [175,31,17],[122,4,3]] },
  { key: 'gray', label: 'Grayscale', stops: [[20,20,20],[245,245,245]] },
];
const clamp01 = (x) => (x < 0 ? 0 : x > 1 ? 1 : x);
function interp(stops, t) {
  t = clamp01(t) * (stops.length - 1);
  const i = Math.min(stops.length - 2, t | 0), f = t - i, a = stops[i], b = stops[i + 1];
  return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f, 255];
}
const LUTS = {};
const lutFor = (cmap) => (LUTS[cmap.key] ??= Array.from({ length: 256 }, (_, i) => interp(cmap.stops, i / 255)));

// Per-cell metrics: which value colors the globes. `file` names the flat
// binary ({sys}_r{res}_{file}.f32); `fmtVal` is the tooltip/hist-line text.
// The metric's shared domain {lo, hi, p99} lives in DOMAINS[key], computed
// exactly over all loaded cells (AR anchors lo at its floor of 1; the area
// ratio spans its observed min..max around the ideal 1).
const METRICS = [
  { key: 'ar', label: 'aspect ratio', file: 'ar', floor: 1,
    fmtVal: (v) => `AR ${v.toFixed(3)}` },
  { key: 'area', label: 'area ratio (cell / ideal)', file: 'area', floor: null,
    fmtVal: (v) => `area ×${v.toFixed(3)}` },
];

// Value transforms: metric value -> t in [0,1] for the COLOR. The value
// axis (histogram + legend position) stays linear over [lo, hi]; only the
// color placement changes, so a stretched scale shows saturation visibly.
const TFS = [
  { key: 'power', label: 'γ0.4 stretch',
    fn: (v, d) => Math.pow(clamp01((v - d.lo) / (d.hi - d.lo)), 0.4) },
  { key: 'linear', label: 'linear in value',
    fn: (v, d) => (v - d.lo) / (d.hi - d.lo) },
  { key: 'p99', label: 'linear, saturates at p99',
    fn: (v, d) => clamp01((v - d.lo) / (d.p99 - d.lo)) },
  { key: 'log', label: 'log',
    fn: (v, d) => (d.hi > d.lo ? Math.log(v / d.lo) / Math.log(d.hi / d.lo) : 0) },
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

const DOMAINS = {};   // per metric key: shared {lo, hi, p99} over all globe cells
const PANELS = [];    // { sys, label, color, orb, layer, vals, hists }
let scale = { metric: METRICS[0], cmap: CMAPS[0], tf: TFS[0] };   // the current pick

const domain = () => DOMAINS[scale.metric.key];

const makeFill = (vals, sc) => {
  const lut = lutFor(sc.cmap), tf = sc.tf.fn, d = DOMAINS[sc.metric.key];
  return (i) => {
    const v = vals[i];
    if (!Number.isFinite(v)) return DNC_GREY;
    return lut[Math.min(255, Math.max(0, (tf(v, d) * 255) | 0))];
  };
};

function drawLegend(sc) {
  const d = DOMAINS[sc.metric.key];
  $('#legLo').textContent = fmt(d.lo, 2);
  $('#legHi').textContent = fmt(d.hi, 2);
  const lut = lutFor(sc.cmap), tf = sc.tf.fn, n = 96;
  const stops = Array.from({ length: n }, (_, i) => {
    const p = i / (n - 1), v = d.lo + p * (d.hi - d.lo);
    const c = lut[Math.min(255, (tf(v, d) * 255) | 0)];
    return `rgb(${c[0] | 0},${c[1] | 0},${c[2] | 0}) ${100 * p}%`;
  }).join(',');
  $('#legendGrad').style.background = `linear-gradient(90deg, ${stops})`;
}

function applyScale() {
  drawLegend(scale);
  for (const p of PANELS)
    p.layer.update({ fill: makeFill(p.vals[scale.metric.key], scale) });
}

// ---- shared domains + per-globe histograms (bars are scale-independent) ----
function computeDomainAndHists() {
  for (const m of METRICS) {
    const all = [];
    for (const p of PANELS)
      for (const v of p.vals[m.key]) if (Number.isFinite(v)) all.push(v);
    all.sort((x, y) => x - y);
    const lo = m.floor ?? (all.length ? all[0] : 0);
    const hi = all.length ? all[all.length - 1] : lo + 1;
    const p99 = all.length ? all[Math.floor(0.99 * (all.length - 1))] : hi;
    DOMAINS[m.key] = { lo, hi, p99 };
    const span = hi - lo || 1;
    for (const p of PANELS) {
      const counts = new Array(HBINS).fill(0);
      for (const v of p.vals[m.key]) {
        if (!Number.isFinite(v)) continue;
        counts[Math.min(HBINS - 1, Math.max(0, Math.floor((v - lo) / span * HBINS)))]++;
      }
      (p.hists ??= {})[m.key] = { counts, cmax: Math.max(1, ...counts) };
    }
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
let histPanel = null;   // last-drawn panel, so a metric switch can redraw it
function drawHist(panel, hlVal) {
  histPanel = panel;
  $('#histLabel').innerHTML = panel
    ? `<span class="swatch" style="background:${panel.color}"></span>${panel.label} — ${panel.vals.ar.length.toLocaleString()} cells`
    : 'hover a globe to see its distribution';
  hctx.clearRect(0, 0, hcw, hch);
  if (!panel) return;
  const d = domain();
  const { counts, cmax } = panel.hists[scale.metric.key];
  const bw = hcw / HBINS, logmax = Math.log(cmax + 1);
  hctx.fillStyle = panel.color;
  for (let i = 0; i < HBINS; i++) {
    if (!counts[i]) continue;
    const h = (Math.log(counts[i] + 1) / logmax) * (hch - 3);
    hctx.fillRect(i * bw, hch - h, Math.max(1, bw - 0.6), h);
  }
  if (hlVal != null && Number.isFinite(hlVal)) {
    const x = clamp01((hlVal - d.lo) / (d.hi - d.lo)) * hcw;
    hctx.strokeStyle = '#fff'; hctx.lineWidth = 1.5;
    hctx.beginPath(); hctx.moveTo(x, 0); hctx.lineTo(x, hch); hctx.stroke();
    hctx.fillStyle = '#fff'; hctx.font = '11px ui-monospace, Menlo, monospace';
    const t = scale.metric.fmtVal(hlVal), tw = hctx.measureText(t).width;
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
  const [pos, starts, ar, area, ids] = await Promise.all([
    fetchBin(`out/globe/${key}_pos.f32`, Float32Array),
    fetchBin(`out/globe/${key}_idx.u32`, Uint32Array),
    fetchBin(`out/globe/${key}_ar.f32`, Float32Array),
    fetchBin(`out/globe/${key}_area.f32`, Float32Array),
    fetch(`out/globe/${key}_ids.json`).then((r) => r.json()),
  ]);
  const vals = { ar, area };

  let amax = 1, rlo = Infinity, rhi = -Infinity;
  for (const a of ar) if (Number.isFinite(a) && a > amax) amax = a;
  for (const v of area) { if (v < rlo) rlo = v; if (v > rhi) rhi = v; }
  card.querySelector('.globe-meta').textContent =
    `${ids.length.toLocaleString()} cells · max AR ${fmt(amax, 3)}`
    + ` · area ×${fmt(rlo, 2)}–${fmt(rhi, 2)}`;

  const layer = orb.polygons({ lnglat: pos, starts,
                               fill: makeFill(vals[scale.metric.key], scale) });
  const panel = { sys, label: `${M.labels[sys]} ${M.res_prefix[sys]}${res}`, color: M.colors[sys], orb, layer, vals };
  PANELS.push(panel);

  orb.on('viewchange', () => syncFrom(orb));

  const tip = $('#tooltip');
  orb.on('hover', (e) => {
    orb.highlight(e.index ?? -1, layer);
    const v = e.index == null ? null : vals[scale.metric.key][e.index];
    drawHist(panel, Number.isFinite(v) ? v : null);   // this globe's distribution + line
    if (e.index == null) { tip.style.opacity = 0; return; }
    const a = ar[e.index];
    tip.innerHTML = `${ids[e.index]}<br>AR ${Number.isFinite(a) ? fmt(a, 4) : 'DNC'}`
      + ` · area ×${fmt(area[e.index], 4)}`;
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

  const metricSel = $('#metricSelect'), cmapSel = $('#cmapSelect'), tfSel = $('#tfSelect');
  const onPick = () => {
    scale = { metric: METRICS.find((m) => m.key === metricSel.value),
              cmap: CMAPS.find((c) => c.key === cmapSel.value),
              tf: TFS.find((t) => t.key === tfSel.value) };
    applyScale();
    drawHist(histPanel ?? PANELS[0], null);   // the hist shows the current metric
  };
  for (const [sel, opts] of [[metricSel, METRICS], [cmapSel, CMAPS], [tfSel, TFS]]) {
    opts.forEach((o, i) =>
      sel.add(new Option(i ? o.label : `${o.label} — the default`, o.key)));
    sel.addEventListener('change', onPick);
  }
  // a fresh <select> starts on its first option, which IS the default

  // provisional domains from the manifest, so the first fill is sane
  DOMAINS.ar = { lo: 1, hi: M.globe_ar_max || 1, p99: M.globe_ar_max || 1 };
  DOMAINS.area = { lo: M.globe_area_min || 1, hi: M.globe_area_max || 1,
                   p99: M.globe_area_max || 1 };
  for (const sys of M.systems) await buildGlobe(M, sys);   // one at a time; 6 WebGL panels
  computeDomainAndHists();
  onPick();                           // apply the selects' current state
}
main();
