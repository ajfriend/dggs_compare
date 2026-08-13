// Static DGGS comparison site: a grid of ajglobe globes, one per system,
// cells colored by the selected METRIC (aspect ratio, or area relative to
// the resolution's exact mean) on a SHARED scale so systems compare directly.
// - four dropdowns pick the shared cell size (the RESOLUTION anchor, each
//   system count-matched to it), the metric, the colormap, and the value
//   transform; globes + legend + histogram update together.
// - hovering a globe shows THAT globe's distribution of the active metric as
//   a histogram above the legend, with a line marking the hovered cell
//   (aligned to the color bar, so you see where the cell sits in both the
//   distribution and the color scale).
// - all globes rotate/zoom together.
import { Orb } from './vendor/ajglobe.min.js';

// ---- the dropdown axes (colormap, transform, metric — METRICS sits
// below the transforms it feeds). One entry per option, in menu order; the
// FIRST entry of each list is the default. Adding an option is one new
// object in its list — the dropdowns and lookups all derive from these.
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

// Value transforms: metric value in [lo, max] -> t in [0,1] for the COLOR.
// The value axis (histogram + legend position) stays linear over [lo, max];
// only the color placement changes, so a stretched scale shows saturation
// visibly. For AR lo == 1 (its floor); for relative area lo is the smallest
// cell across the globes.
const TFS = [
  { key: 'power', label: 'γ0.4 stretch',
    fn: (v, d) => Math.pow((v - d.lo) / (d.max - d.lo), 0.4) },
  { key: 'linear', label: 'linear in the value',
    fn: (v, d) => (v - d.lo) / (d.max - d.lo) },
  { key: 'p99', label: 'linear, saturates at p99',
    fn: (v, d) => clamp01((v - d.lo) / (d.p99 - d.lo)) },
  { key: 'log', label: 'log',
    fn: (v, d) => (d.max > d.lo ? Math.log(v / d.lo) / Math.log(d.max / d.lo) : 0) },
];

// The metric axis: which per-cell value the globes color by. `key` doubles
// as the binary suffix (out/globe/*_{key}.f32) and the panel channel name;
// `name` prefixes tooltips/histograms; `stat` renders the per-globe line;
// `fixedLo` pins the domain floor (AR starts at 1 by definition; area's
// floor is data-driven).
const METRICS = [
  { key: 'ar', label: 'Aspect ratio', name: 'AR', fixedLo: 1,
    stat: (lo, hi) => `max AR ${fmt(hi, 3)}` },
  { key: 'area', label: 'Relative area (cell / mean)', name: 'area',
    fixedLo: null, stat: (lo, hi) => `area ${fmt(lo, 3)}–${fmt(hi, 3)}` },
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

const DOMAIN = { lo: 1, max: 1, p99: 1 };  // shared domain of the ACTIVE metric
const PANELS = [];   // { sys, label, color, orb, layer, vals, ids, hist, meta, resEl, byAnchor }
let scale = { cmap: CMAPS[0], tf: TFS[0] };   // the current pick (defaults)
let metric = METRICS[0];
let anchorIdx = 0;   // active Resolution-menu index into globe_res lists
let nSystems = 0;    // panels expected once fully built (set from the manifest)

const vals = (p) => p.vals[metric.key];

const makeFill = (v, sc) => {
  const lut = lutFor(sc.cmap), tf = sc.tf.fn;
  return (i) => {
    const a = v[i];
    if (!Number.isFinite(a)) return DNC_GREY;
    return lut[Math.min(255, Math.max(0, (tf(a, DOMAIN) * 255) | 0))];
  };
};

function drawLegend(sc) {
  $('#legLo').textContent = fmt(DOMAIN.lo, 2);
  $('#legHi').textContent = fmt(DOMAIN.max, 2);
  const lut = lutFor(sc.cmap), tf = sc.tf.fn, n = 96;
  const stops = Array.from({ length: n }, (_, i) => {
    const p = i / (n - 1), v = DOMAIN.lo + p * (DOMAIN.max - DOMAIN.lo);
    const c = lut[Math.min(255, (tf(v, DOMAIN) * 255) | 0)];
    return `rgb(${c[0] | 0},${c[1] | 0},${c[2] | 0}) ${100 * p}%`;
  }).join(',');
  $('#legendGrad').style.background = `linear-gradient(90deg, ${stops})`;
}

function applyScale() {
  drawLegend(scale);
  for (const p of PANELS)
    if (p.layer) p.layer.update({ fill: makeFill(vals(p), scale) });
}

// ---- shared metric domain + per-globe histograms (scale-independent).
// Memoized per (metric, anchor): the data is static, so revisiting a
// combination is a lookup, not a ~400k-value re-sort. The domain
// recomputes per anchor deliberately — each resolution view is
// self-consistent rather than color-comparable across resolutions.
const DOMAINS = {};   // `${metric}@${anchor}` -> { lo, max, p99, hists: Map }
function computeDomainAndHists() {
  const key = `${metric.key}@${anchorIdx}`;
  let c = DOMAINS[key];
  if (!c) {
    const ready = PANELS.filter((p) => p.vals);   // mid-build: some pending
    const all = [];
    for (const p of ready) for (const a of vals(p)) if (Number.isFinite(a)) all.push(a);
    const s = Float32Array.from(all).sort();   // typed: numeric, no comparator
    const lo = metric.fixedLo ?? (s.length ? s[0] : 0);
    const max = s.length ? s[s.length - 1] : lo + 1;
    const p99 = s.length ? s[Math.floor(0.99 * (s.length - 1))] : max;
    const span = max - lo || 1;
    const hists = new Map();
    for (const p of ready) {
      const counts = new Array(HBINS).fill(0);
      for (const a of vals(p)) {
        if (!Number.isFinite(a)) continue;
        counts[Math.min(HBINS - 1, Math.max(0, Math.floor((a - lo) / span * HBINS)))]++;
      }
      hists.set(p, { counts, cmax: Math.max(1, ...counts) });
    }
    c = { lo, max, p99, hists };
    // memoize only complete views — a partial domain must not masquerade
    // as this (metric, anchor) combination's answer on later visits
    if (ready.length === nSystems) DOMAINS[key] = c;
  }
  Object.assign(DOMAIN, { lo: c.lo, max: c.max, p99: c.p99 });
  for (const p of PANELS) p.hist = c.hists.get(p);
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
function drawHist(panel, hl) {
  $('#histLabel').innerHTML = panel
    ? `<span class="swatch" style="background:${panel.color}"></span>${panel.label} — ${vals(panel).length.toLocaleString()} cells`
    : 'hover a globe to see its distribution';
  hctx.clearRect(0, 0, hcw, hch);
  if (!panel?.hist) return;
  const { counts, cmax } = panel.hist, bw = hcw / HBINS, logmax = Math.log(cmax + 1);
  hctx.fillStyle = panel.color;
  for (let i = 0; i < HBINS; i++) {
    if (!counts[i]) continue;
    const h = (Math.log(counts[i] + 1) / logmax) * (hch - 3);
    hctx.fillRect(i * bw, hch - h, Math.max(1, bw - 0.6), h);
  }
  if (hl != null && Number.isFinite(hl)) {
    const x = clamp01((hl - DOMAIN.lo) / (DOMAIN.max - DOMAIN.lo)) * hcw;
    hctx.strokeStyle = '#fff'; hctx.lineWidth = 1.5;
    hctx.beginPath(); hctx.moveTo(x, 0); hctx.lineTo(x, hch); hctx.stroke();
    hctx.fillStyle = '#fff'; hctx.font = '11px ui-monospace, Menlo, monospace';
    const t = `${metric.name} ${hl.toFixed(3)}`, tw = hctx.measureText(t).width;
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

function summaryTables() {
  // The two renderings of the working-resolution summary stats (fragments
  // written by survey.py); a missing file (an older cached build) just
  // leaves its block empty.
  for (const [id, file] of [['#statsTable', 'summary_table.html'],
                            ['#statsTableGT', 'summary_table_gt.html']]) {
    fetch(`out/${file}`)
      .then((r) => (r.ok ? r.text() : ''))
      .then((html) => { $(id).innerHTML = html; });
  }
}

function buildGlobe(M, sys) {   // card + Orb + hover wiring; data comes via loadAnchor
  const card = document.createElement('div');
  card.className = 'globe-card';
  card.innerHTML =
    `<div class="globe-label"><span class="swatch" style="background:${M.colors[sys]}"></span>`
    + `${M.labels[sys]} <span class="res"></span></div>`
    + `<div class="globe-holder"><canvas></canvas></div><div class="globe-meta"></div>`;
  $('#globeGrid').appendChild(card);

  const canvas = card.querySelector('canvas');
  const orb = new Orb(canvas, { background: '#0b0e13', sphere: '#11151c' });
  orb.lookAt(0, 20);
  orb.borders({ color: '#2a3340', width: 1 });

  // Data (and the polygon layer) attach per anchor via loadAnchor; the
  // hover handlers read the panel's ACTIVE fields, so they bind once.
  const panel = { sys, label: '', color: M.colors[sys], orb,
                  layer: null, vals: null, ids: null, byAnchor: {},
                  meta: card.querySelector('.globe-meta'),
                  resEl: card.querySelector('.res') };
  PANELS.push(panel);

  orb.on('viewchange', () => syncFrom(orb));

  const tip = $('#tooltip');
  orb.on('hover', (e) => {
    orb.highlight(e.index ?? -1, panel.layer);
    const a = e.index == null ? null : vals(panel)[e.index];
    drawHist(panel, Number.isFinite(a) ? a : null);   // this globe's distribution + line
    if (e.index == null) { tip.style.opacity = 0; return; }
    tip.innerHTML = `${panel.ids[e.index]}<br>${metric.name} ${Number.isFinite(a) ? fmt(a, 4) : 'DNC'}`;
    const r = canvas.getBoundingClientRect();
    tip.style.left = `${r.left + e.x + 14}px`;
    tip.style.top = `${r.top + e.y + 14}px`;
    tip.style.opacity = 1;
  });
  canvas.addEventListener('pointerleave', () => { tip.style.opacity = 0; drawHist(panel, null); });
  return panel;
}

// Fetch (once, cached per panel) and activate the CURRENT anchor's data
// on `p`: swap the polygon layer and point the panel's active fields at it.
async function loadAnchor(M, p) {
  let d = p.byAnchor[anchorIdx];
  if (!d) {
    const res = M.globe_res[p.sys][anchorIdx];
    const key = `${p.sys}_r${res}`;
    const [pos, starts, ar, area, ids] = await Promise.all([
      fetchBin(`out/globe/${key}_pos.f32`, Float32Array),
      fetchBin(`out/globe/${key}_idx.u32`, Uint32Array),
      fetchBin(`out/globe/${key}_ar.f32`, Float32Array),
      // Tolerate a cached pre-area build (see the guard in main()).
      fetchBin(`out/globe/${key}_area.f32`, Float32Array).catch(() => null),
      fetch(`out/globe/${key}_ids.json`).then((r) => r.json()),
    ]);
    d = p.byAnchor[anchorIdx] = { res, pos, starts, vals: { ar, area }, ids };
  }
  p.layer?.remove();
  Object.assign(p, { vals: d.vals, ids: d.ids });
  p.layer = p.orb.polygons({ lnglat: d.pos, starts: d.starts,
                             fill: makeFill(vals(p) ?? d.vals.ar, scale) });
  const rlabel = `${M.res_prefix[p.sys]}${d.res}`;
  p.resEl.textContent = rlabel;
  p.label = `${M.labels[p.sys]} ${rlabel}`;
}

// Per-globe stat line under each canvas, restated for the active metric.
function setMeta(p) {
  if (!p.vals) return;   // mid-build: not loaded yet
  let lo = Infinity, hi = -Infinity;
  for (const a of vals(p)) if (Number.isFinite(a)) { if (a < lo) lo = a; if (a > hi) hi = a; }
  p.meta.textContent = `${p.ids.length.toLocaleString()} cells · ${metric.stat(lo, hi)}`;
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

// Plot-group tabs. The active tab rides in the URL hash (#globes / #shape /
// #area / #tradeoff) so tabs are linkable and survive reload; unknown or
// absent hashes fall back to the first tab. `onShow` fires with the tab
// name after every switch (including the initial one).
function initTabs(onShow) {
  const tabs = [...document.querySelectorAll('.tab')];
  const show = (name) => {
    if (!tabs.some((t) => t.dataset.tab === name)) name = tabs[0].dataset.tab;
    for (const t of tabs) t.classList.toggle('active', t.dataset.tab === name);
    for (const p of document.querySelectorAll('.tab-panel'))
      p.classList.toggle('active', p.id === `panel-${name}`);
    onShow(name);
  };
  for (const t of tabs)
    t.addEventListener('click', () => {
      history.replaceState(null, '', `#${t.dataset.tab}`);
      show(t.dataset.tab);
    });
  window.addEventListener('hashchange', () => show(location.hash.slice(1)));
  show(location.hash.slice(1));
}

async function main() {
  summaryTables();   // independent of the manifest — don't wait on it
  const M = await fetch('out/manifest.json').then((r) => r.json());
  $('#subtitle').textContent =
    `${M.systems.length} systems · shape via csar's enclosing-cone solver, area via sparea`;
  if (M.tag) $('#tag').textContent = M.tag;
  if (M.gap_tol) $('#gap').textContent = M.gap_tol.toExponential();

  byResGrid(M);
  initLightbox();

  // Everything globe-panel runs lazily and RESUMABLY, driven by showings
  // of the globes tab: canvases (the histogram's and each Orb's) size
  // themselves from the DOM at construction, which only works while
  // their panel is displayed. So the one-time setup happens on the
  // first showing, and the build loop checks visibility before each
  // globe — switching away mid-build pauses it, switching back resumes.
  // Normalize an older cached manifest (one resolution per system, no
  // anchor labels) into the per-anchor list form.
  if (!Array.isArray(M.globe_res[M.systems[0]])) {
    for (const s of M.systems) M.globe_res[s] = [M.globe_res[s]];
  }
  M.globe_anchor_labels ??= [''];
  anchorIdx = M.globe_res[M.systems[0]].length - 1;   // default: the finest
  nSystems = M.systems.length;

  const panelShown = () => $('#panel-globes').classList.contains('active');
  let onMetric, resSel, wired = false, building = false, finalized = false, next = 0;
  const ensureGlobes = async () => {
    if (!wired) {
      wired = true;
      const metricSel = $('#metricSelect'), cmapSel = $('#cmapSelect'), tfSel = $('#tfSelect');
      const onPick = () => {
        scale = { cmap: CMAPS.find((c) => c.key === cmapSel.value),
                  tf: TFS.find((t) => t.key === tfSel.value) };
        applyScale();
      };
      // The metric changes the DOMAIN and the histograms, not just colors.
      onMetric = () => {
        metric = METRICS.find((m) => m.key === metricSel.value);
        computeDomainAndHists();
        for (const p of PANELS) setMeta(p);
        onPick();
        drawHist(PANELS[0], null);
      };
      for (const [sel, opts, on] of [[metricSel, METRICS, onMetric],
                                     [cmapSel, CMAPS, onPick], [tfSel, TFS, onPick]]) {
        opts.forEach((o) => sel.add(new Option(o.label, o.key)));
        sel.addEventListener('change', on);
      }
      // a fresh <select> starts on its first option, which IS the default
      // The anchor changes every panel's DATA: all four menus disable for
      // the swap's duration (and until the initial build finishes) so no
      // recolor or domain computation can run over mixed-anchor panels.
      resSel = $('#resSelect');
      M.globe_anchor_labels.forEach((l) => resSel.add(new Option(l)));
      resSel.selectedIndex = anchorIdx;
      resSel.disabled = true;
      if (resSel.options.length < 2) resSel.closest('label').hidden = true;
      const allSels = [resSel, metricSel, cmapSel, tfSel];
      resSel.addEventListener('change', async () => {
        for (const x of allSels) x.disabled = true;
        anchorIdx = resSel.selectedIndex;
        try {
          await Promise.all(PANELS.map((p) => loadAnchor(M, p)));
          onMetric();   // per-anchor domain/hists/legend/meta
        } finally {
          for (const x of allSels) x.disabled = false;
        }
      });
      initHist();
      DOMAIN.max = M.globe_ar_max || 1;   // provisional, so the first fill is sane
      DOMAIN.p99 = DOMAIN.max;
    }
    if (building) return;
    building = true;
    while (next < M.systems.length && panelShown()) {   // serial builds; pause point
      const p = buildGlobe(M, M.systems[next]);
      await loadAnchor(M, p);
      setMeta(p);
      next++;
    }
    building = false;
    if (next === M.systems.length && !finalized) {
      finalized = true;
      // A cached pre-area build has no area binaries: keep the option
      // visible but unusable, rather than silently showing wrong data.
      if (PANELS.some((p) => !p.vals.area)) {
        $('#metricSelect').querySelector('option[value="area"]').disabled = true;
      }
      onMetric();   // seed domain/hists/legend/meta for the default metric
      resSel.disabled = false;
    }
  };

  initTabs((name) => { if (name === 'globes') ensureGlobes(); });
}
main();
