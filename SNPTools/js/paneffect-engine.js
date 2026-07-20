/* =====================================================================
 *  paneffect-engine.js — native PanEffect engine for SNPTools.
 *
 *  PanEffect's scripts (support.js, genome.js, pan.js, dom.js, main.js) are
 *  classic scripts, so they all share ONE global scope. This engine (also a
 *  classic script) therefore reads and writes their globals directly
 *  (gene_model, main_option, gene_model_length, GN_array, …) and reassigns
 *  renderHeatmap / renderHeatmapPan to canvas versions.
 *
 *  Load order (see index.html):
 *    d3  ->  support.js -> genome.js -> pan.js -> dom.js -> main.js
 *        ->  paneffect-heatmap.js -> paneffect-engine.js -> snppaneffect.js
 *
 *  Exposes window.PanEffectEngine.render(container, opts) / .teardown().
 *    opts = { gene, esm, option, variant:{pos,wt,sub}, dataBase }
 * ===================================================================== */
(function () {
  'use strict';
  if (typeof window === 'undefined') return;

  var DATA_DIRS = ['csv','heatmap','target','query','pfam','uniprot','synonym','traits','dssp','structures'];
  var ESM_OK = ['ESM1','ESM2','ESM3'];

  var state = {
    installedFetch: false, origFetch: null,
    base: './paneffect/', variant: null,
    grids: { b73: null, pan: null },
  };

  /* ---------------- fetch rebaser -------------------------------------
     Rewrites only relative paths whose first segment is a known PanEffect
     data directory, so SNPTools' own requests are never touched. */
  function normalizeBase(b) {
    b = String(b || './paneffect/');
    if (b.slice(-1) !== '/') b += '/';
    return b;
  }
  function rebase(url) {
    if (typeof url !== 'string') return url;
    if (/^[a-z]+:\/\//i.test(url) || url.indexOf('//') === 0) return url; // absolute
    var u = url.replace(/^\.\//, '').replace(/^\//, '');
    var seg = u.split('/')[0].split('?')[0];
    if (DATA_DIRS.indexOf(seg) === -1) return url;
    return state.base + u;
  }
  function installFetch() {
    if (state.installedFetch) return;
    var orig = window.fetch.bind(window);
    state.origFetch = orig;
    window.fetch = function (input, init) {
      try {
        var url = (typeof input === 'string') ? input : (input && input.url);
        if (typeof url === 'string') {
          var r = rebase(url);
          if (r !== url) input = (typeof input === 'string') ? r : new Request(r, input);
        }
      } catch (e) { /* fall through to original */ }
      return orig(input, init);
    };
    state.installedFetch = true;
  }

  /* ---------------- id parsing (writes shared globals) --------------- */
  function parseId(id) {
    id = String(id || '').trim();
    if (id.indexOf('_T') >= 0)      { gene_model = id.split('_')[0]; transcript = id;                   protein = id.replace('_T', '_P'); can_flag = false; }
    else if (id.indexOf('_P') >= 0) { gene_model = id.split('_')[0]; transcript = id.replace('_P','_T'); protein = id;                     can_flag = false; }
    else                            { gene_model = id; can_flag = true; transcript = id + '_T001';        protein = id + '_P001'; }
  }

  /* ---------------- canvas heatmap renderers (override) -------------- */
  var Y_ROWS = 20;
  function subToRow(sub) {
    if (!sub) return null;
    sub = String(sub).toUpperCase();
    for (var y = 1; y <= Y_ROWS; y++) { if (numberToAminoAcid(y) === sub) return y; }
    return null;
  }

  function canvasRenderGene(data) {
    var host = document.getElementById('full-heatmap');
    if (!host || !window.PanEffectHeatmap) return;
    if (state.grids.b73) { state.grids.b73.destroy(); state.grids.b73 = null; }
    Array.prototype.forEach.call(host.querySelectorAll('canvas.pe-hm'), function (c) { c.remove(); });

    var cells = data.map(function (d) {
      return { x: +d.X, y: +d.Y, score: +d.Score, wt: d.WT, wgs2024: +d.WGS2024, wgs2026: +d.WGS2026 };
    });
    var vis = function (c) {
      return wgs_status || (wgs2024_status && c.wgs2024 === 1) || (wgs2026_status && c.wgs2026 === 1);
    };
    var hl = state.variant && state.variant.pos ? { x: state.variant.pos, y: subToRow(state.variant.sub) } : null;

    var grid = PanEffectHeatmap.renderGrid(host, cells, {
      cols: gene_model_length, rows: Y_ROWS,
      cellW: window_length / gene_model_length, cellH: 14,
      color: function (score, c) { return vis(c) ? colorScale(score) : '#FFFFFF'; },
      tooltip: function (c) {
        return 'Position: ' + c.x + '<br>Substitution: ' + c.wt + ' &rarr; ' +
               numberToAminoAcid(c.y) + '<br>Score: ' + c.score;
      },
      highlight: (hl && hl.y) ? hl : null,
    });
    grid.canvas.classList.add('pe-hm');
    state.grids.b73 = grid;
    legend(host, '.heatmap-container', 'pe-legend-gene');
  }

  function canvasRenderPan(data) {
    var host = document.getElementById('full-heatmap-pan');
    if (!host || !window.PanEffectHeatmap) return;
    if (state.grids.pan) { state.grids.pan.destroy(); state.grids.pan = null; }
    Array.prototype.forEach.call(host.querySelectorAll('canvas.pe-hm'), function (c) { c.remove(); });

    var rows = GN_size || 0;
    var cells = data.map(function (d) {
      var wt = String(d.WT == null ? '' : d.WT).trim();
      var sub = String(d.Sub == null ? '' : d.Sub).trim();
      var raw = d.Score;
      var sc = (raw === '' || raw == null) ? NaN : +raw;
      var gap = wt.charAt(0) === '-' || sub.charAt(0) === '-' || raw === '' || raw == null || isNaN(sc);
      if (!rows && +d.Y > rows) rows = +d.Y;
      return { x: +d.X, y: +d.Y, score: sc, gap: gap, wt: d.WT, sub: d.Sub, x2: d.X2, x3: d.X3 };
    });

    var grid = PanEffectHeatmap.renderGrid(host, cells, {
      cols: alignment_length, rows: rows || 1,
      cellW: window_length / alignment_length, cellH: 5,
      color: function (score, c) { return c.gap ? '#e6e6e6' : colorScalePan(score); },
      tooltip: function (c) {
        return 'B73 Position: ' + (c.x2 == null ? '' : c.x2) +
               '<br>Target Position: ' + (c.x3 == null ? '' : c.x3) +
               '<br>Genome: ' + (GN_array[c.y] == null ? '' : GN_array[c.y]) +
               '<br>G.M.: ' + (GM_array[c.y] == null ? '' : GM_array[c.y]) +
               '<br>Substitution: ' + c.wt + ' to ' + c.sub +
               '<br>Score: ' + (isNaN(c.score) ? '' : c.score);
      },
    });
    grid.canvas.classList.add('pe-hm');
    state.grids.pan = grid;
    var box = document.getElementById('colorBox-pan');
    if (box) { box.innerHTML = ''; box.appendChild(buildLegendBar()); }
  }

  /* compact horizontal score legend */
  var LEGEND_COLORS = ['#00429d','#3860aa','#587fb3','#78a0b7','#9ac0b3','#c1e19e',
    '#ffff00','#ffd337','#fea447','#f1784d','#db4c4d','#bd2147','#93003a'];
  function buildLegendBar() {
    var wrap = document.createElement('div');
    wrap.className = 'pe-legend';
    var bar = document.createElement('div'); bar.className = 'pe-legend-bar';
    LEGEND_COLORS.forEach(function (col) {
      var s = document.createElement('span'); s.style.background = col; bar.appendChild(s);
    });
    var labels = document.createElement('div'); labels.className = 'pe-legend-labels';
    labels.innerHTML = '<span>&gt; 0 (benign)</span><span>&minus;11</span><span>&lt; &minus;22 (strong)</span>';
    wrap.appendChild(bar); wrap.appendChild(labels);
    return wrap;
  }
  function legend(hostEl, containerSel, cls) {
    var container = hostEl.closest ? hostEl.closest(containerSel) : null;
    if (!container || !container.parentNode) return;
    var prev = container.parentNode.querySelector('.' + cls);
    if (prev) prev.remove();
    var el = buildLegendBar(); el.classList.add(cls);
    container.parentNode.insertBefore(el, container.nextSibling);
  }

  /* ---------------- view visibility --------------------------------- */
  function applyVisibility(option) {
    var sm = document.getElementById('summary');
    var b = document.getElementById('b73');
    var p = document.getElementById('pan-genome');
    if (sm) sm.classList.add('active');
    if (b) b.classList.toggle('active', option !== 'pan');
    if (p) p.classList.toggle('active', option !== 'b73');
  }

  function centerOnVariant(pos) {
    var s = document.getElementById('zoom-slider');
    if (!s) return;
    var start = Math.max(1, Math.round(pos - 25));
    var max = +s.max || start;
    s.value = Math.min(max, start);
    s.dispatchEvent(new Event('input'));
  }

  /* ---------------- WGS view default -------------------------------- */
  function applyWgs(mode) {
    var is2026 = (mode === 'maize2026');
    wgs_status = !is2026;         /* 'all variants' unless a 2026 handoff */
    wgs2024_status = false;
    wgs2026_status = is2026;
    var rAll = document.getElementById('allVariants');
    var r26 = document.getElementById('maizeWGS2026');
    if (r26) r26.checked = is2026;
    if (rAll) rAll.checked = !is2026;
  }

  /* ---------------- clean summary (override) ------------------------ */
  function cleanSummary() {
    var gnPrint = (typeof gn === 'string' && gn.trim() && gn.trim() !== 'N/A') ? ' (' + gn.trim() + ')' : '';
    var file = (typeof gene_model_file !== 'undefined' && gene_model_file)
      ? gene_model_file : ('./csv/' + currentESM + '/' + protein + '.csv');
    var csv = state.base + String(file).replace(/^\.\//, '');
    var sm = document.getElementById('summary');
    if (sm) {
      sm.innerHTML =
        '<div class="pe-sum">' +
          '<div class="pe-sum-kicker">Missense variant effects · ' + currentESM + ' · B73 v5</div>' +
          '<h2 class="pe-sum-title">' + (gene_model || '') + gnPrint + '</h2>' +
          '<a class="pe-sum-dl" href="' + csv + '" download>' +
            '<svg viewBox="0 0 24 24" width="16" height="16" fill="none" aria-hidden="true">' +
            '<path d="M12 4v10m0 0l-4-4m4 4l4-4M5 20h14" stroke="currentColor" stroke-width="2" ' +
            'stroke-linecap="round" stroke-linejoin="round"/></svg>' +
            'Download variant effects file</a>' +
        '</div>';
    }
    var r = document.getElementById('reference_gm'); if (r) r.innerHTML = 'B73 Reference View';
    var p = document.getElementById('pan_gm'); if (p) p.innerHTML = 'Pan-genome View';
  }

  /* ---------------- teardown ---------------------------------------- */
  function teardown() {
    if (state.grids.b73) { try { state.grids.b73.destroy(); } catch (e) {} state.grids.b73 = null; }
    if (state.grids.pan) { try { state.grids.pan.destroy(); } catch (e) {} state.grids.pan = null; }
    /* remove any leaked tooltips from the DOM zoomed views + stray canvas tips */
    Array.prototype.forEach.call(
      document.querySelectorAll('body > .tooltip, body > .tooltipdomain, .pe-canvas-tip'),
      function (n) { n.remove(); }
    );
  }

  /* ---------------- render ------------------------------------------ */
  function render(container, opts) {
    opts = opts || {};
    return (async function () {
      teardown();

      /* Fail loudly & clearly if the ported PanEffect scripts aren't loaded.
         typeof is ReferenceError-safe for undeclared identifiers, unlike a
         bare read/write of main_id. */
      if (typeof runPanEffect !== 'function' ||
          typeof renderHeatmap !== 'function' ||
          typeof renderHeatmapPan !== 'function' ||
          typeof colorScale !== 'function' ||
          typeof numberToAminoAcid !== 'function') {
        var err = new Error(
          'PanEffect scripts not loaded. Expected ./js/pe/{support,genome,pan,dom,main}.js ' +
          '(main.js must be the edited version that defines runPanEffect). ' +
          'Check the <script src> paths in index.html and the Network tab for 404s.');
        err.code = 'PE_SCRIPTS_MISSING';
        console.error('[PanEffectEngine] ' + err.message);
        throw err;
      }

      installFetch();
      state.base = normalizeBase(opts.dataBase || state.base);
      state.variant = opts.variant || null;

      /* drive PanEffect's shared globals */
      main_id = opts.gene || '';
      main_option = opts.option || 'both';
      currentESM = ESM_OK.indexOf(opts.esm) >= 0 ? opts.esm : 'ESM2';
      parseId(opts.gene || '');

      /* swap the two full heatmaps to canvas (perf) */
      renderHeatmap = canvasRenderGene;
      renderHeatmapPan = canvasRenderPan;
      /* clean, minimal summary + skip the GWAS traits fetch (that lives in SNPFunction) */
      populateSummary = cleanSummary;
      loadAndDisplayTraits = function () {};

      /* inject the view skeleton */
      if (container.classList) container.classList.add('pe-native-root');
      container.innerHTML = skeleton();
      applyVisibility(main_option);

      /* external entries (SNPVersity / SNPFold) open in the MaizeGDB 2026 view */
      applyWgs(opts.wgs);

      /* run the (now callable) PanEffect pipeline */
      var ok = false;
      try { ok = await runPanEffect(); }
      catch (e) { console.error('[PanEffectEngine] pipeline error:', e); }

      if (state.variant && state.variant.pos) {
        /* let the sliders finish wiring, then centre on the variant */
        setTimeout(function () { centerOnVariant(state.variant.pos); }, 0);
      }
      return ok;
    })();
  }

  /* ---------------- skeleton ---------------------------------------- */
  function skeleton() {
    return '' +
'<div id="loading-icon" style="display:none">Loading…</div>' +
'<div id="reference_gm_top" class="pe-gm-head"></div>' +

/* summary content is rendered by the clean populateSummary override */
'<div id="summary" class="content"></div>' +

'<div id="b73" class="content">' +
  '<span class="gene" id="reference_gm"></span>' +
  '<span id="wgs_span"><br><br>' +
    '<input type="radio" id="allVariants" name="variantEffect" value="all" checked>' +
    '<label for="allVariants">Show all variant effects</label>' +
    /* MaizeGDB 2024 option disabled for now
    '<input type="radio" id="maizeWGS" name="variantEffect" value="maize2024">' +
    '<label for="maizeWGS">MaizeGDB 2024 High Coverage variant effects</label>' +
    */
    '<input type="radio" id="maizeWGS2026" name="variantEffect" value="maize2026">' +
    '<label for="maizeWGS2026">MaizeGDB 2026 High Coverage variant effects</label><br>' +
  '</span><br>' +
  '<div class="sectionHeader">PFAM Domains</div>' +
  '<div id="pfam-wrap"><div id="pfamNumberLine" class="numberLine"></div><div id="pfamGeneModel" class="geneModel"></div></div>' +
  '<div class="sectionHeader">Secondary Structure</div>' +
  '<div id="dssp"><canvas id="proteinStructure" width="1400" height="80"></canvas></div>' +
  '<div class="sectionHeader">Variant Effects of full protein (heatmap)</div>' +
  '<div id="heatNumberLine" class="numberLine"></div>' +
  '<div class="heatmap-container"><div id="full-heatmap"><div id="highlight-box"></div></div></div>' +
  '<div class="slider-container" id="slider-container"><span id="slider"></span><span id="slider-value">1</span></div>' +
  '<div class="sectionHeader">Variant Effects of zoomed in region (heatmap)</div>' +
  '<div id="zoomNumberLine" class="numberLine"></div>' +
  '<div id="heatmap-container-zoom" class="heatmap-container-zoom"><div id="zoomed-heatmap"></div></div>' +
  '<div id="zoomWTLine" class="numberLine"></div>' +
'</div>' +

'<div id="pan-genome" class="content">' +
  '<span class="gene" id="pan_gm"></span>' +
  '<div class="sectionHeader">PFAM Domains</div>' +
  '<div id="pfam-wrap-pan"><div id="pfamNumberLine-pan" class="numberLine"></div><div id="pfamGeneModel-pan" class="geneModel"></div></div>' +
  '<div class="sectionHeader">Secondary Structure</div>' +
  '<div id="dssp-pan"><canvas id="proteinStructure-pan" width="1400" height="80"></canvas></div>' +
  '<div class="sectionHeader" id="heatheader-pan">Heatmap of full protein</div>' +
  '<div id="heatNumberLine-pan" class="numberLine"></div>' +
  '<div id="colorBox-pan"></div>' +
  '<div class="heatmap-container-pan"><div id="full-heatmap-pan"><div id="highlight-box-pan"></div></div></div>' +
  '<div class="slider-container" id="slider-container-pan"><span id="slider-pan"></span><span id="slider-value-pan">1</span></div>' +
  '<div class="sectionHeader">Heatmap of zoomed in region</div>' +
  '<div id="zoomNumberLine-pan" class="numberLine"></div>' +
  '<div id="heatmap-container-zoom-pan" class="heatmap-container-zoom-pan"><div id="zoomed-heatmap-pan"></div></div>' +
  '<div id="zoomWTLine-pan" class="numberLine"></div>' +
'</div>';
  }

  window.PanEffectEngine = { render: render, teardown: teardown };
})();
