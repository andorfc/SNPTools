/* =====================================================================
 *  snpfunction-ontology.js — SNPFunction's Gene Ontology and Pathways
 *  views, read live from the MaizeGDB record API.
 *
 *  Ported from MaizeGDB's gene record ("Function at a glance",
 *  mgdb_design src/js/mgdb-gene-function.js as of 2026-09-30) so the two
 *  sites draw a gene's function the same way:
 *    - Gene Ontology: one row per aspect with its ancestry graph, the
 *      plant GO-slim fingerprint (a fixed strip of squares, names at 45
 *      degrees), and the terms with their evidence beneath. Squares,
 *      names, term rows and graph nodes highlight each other on hover and
 *      pin on click.
 *    - Pathways: the pathway explorer's E2P2 pathways as step strips
 *      (this gene's step in gold, NAM-founder coverage under each step,
 *      founder presence dots), then EnTAP's KEGG and EggNOG block.
 *  Left out on purpose: MaizeGDB's protein family / domain-atlas block
 *  (SNPFunction draws the domain architecture itself).
 *
 *  Differences from the MaizeGDB file: the two views sit behind a switch
 *  instead of one long column, the four tiles count GO and pathways only,
 *  and every MaizeGDB link the API sends relative is made absolute and
 *  opens in a new tab, so following one never drops the SNPTools session.
 *
 *  Data: GET <base>/api/v1/records/gene/<id>?fields=function, the same
 *  payload the MaizeGDB page draws. The API sends
 *  Access-Control-Allow-Origin: *, so no proxy is needed.
 *
 *  window.SNPFunctionOntology = {
 *    base,                       // MaizeGDB origin the API is read from
 *    fetch(gene) -> Promise<{status:'ok'|'missing', fn, id, url}>
 *                                // rejects on network / server failure
 *    render(container, spec) -> {setView(v), view(), destroy()} | null
 *      spec = { gene, fn, base, view:'go'|'pathways', onView(v) }
 *    recordUrl(gene), apiUrl(gene)
 *  }
 * ===================================================================== */
(function (window, document) {
  'use strict';

  /* The MaizeGDB origin to read from. claude.maizegdb.org is the
     development instance of the new site; alpha/beta.maizegdb.org serve
     the same API (beta has no KEGG block yet) and www.maizegdb.org has no
     /api/v1 until the new site goes live there. Override before this
     script loads: window.SNPTOOLS_MAIZEGDB_BASE = 'https://…' ('' = the
     page's own origin, once SNPTools is served from a MaizeGDB host). */
  var BASE = (typeof window.SNPTOOLS_MAIZEGDB_BASE === 'string')
    ? window.SNPTOOLS_MAIZEGDB_BASE.replace(/\/+$/, '')
    : 'https://claude.maizegdb.org';
  var TIMEOUT_MS = 20000;

  var SVG_NS = 'http://www.w3.org/2000/svg';

  var ASPECTS = [
    { key: 'biological_process', label: 'Biological process', short: 'BP', cls: 'bp' },
    { key: 'molecular_function', label: 'Molecular function', short: 'MF', cls: 'mf' },
    { key: 'cellular_component', label: 'Cellular component', short: 'CC', cls: 'cc' }
  ];
  var ASPECT_BY_KEY = {};
  ASPECTS.forEach(function (a) { ASPECT_BY_KEY[a.key] = a; });

  var EVIDENCE_KIND = {};
  [['experimental', ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP', 'HGI', 'HEP']],
   ['similarity', ['ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'IBA', 'IBD', 'IKR', 'IRD', 'RCA']],
   ['author', ['TAS', 'NAS', 'IC']],
   ['computational', ['IEA', 'COMP']],
   ['none', ['ND']]].forEach(function (pair) {
    pair[1].forEach(function (code) { EVIDENCE_KIND[code] = pair[0]; });
  });
  var EVIDENCE_LABEL = { experimental: 'experimental', similarity: 'by similarity', author: 'author statement',
                         computational: 'computational', none: 'no data', unknown: 'unstated' };

  /* ---------------- the API ---------------- */
  var cache = {};
  function apiUrl(gene) {
    return BASE + '/api/v1/records/gene/' + encodeURIComponent(gene) + '?fields=function';
  }
  function recordUrl(gene) {
    return BASE + '/gene_center/gene/' + encodeURIComponent(gene);
  }
  /* One request per gene, kept for the session. A 404 is an answer
     ("MaizeGDB has no such gene"), not a failure; anything else that is
     not a 200 with a function section rejects, and the caller falls back
     to SNPTools' own annotation file. */
  function fetchFunction(gene) {
    var key = String(gene || '').trim();
    if (!key) { return Promise.reject(new Error('no gene')); }
    if (cache[key]) { return cache[key]; }
    var ctl = (typeof AbortController === 'function') ? new AbortController() : null;
    var timer = ctl ? setTimeout(function () { ctl.abort(); }, TIMEOUT_MS) : null;
    var url = apiUrl(key);
    var p = fetch(url, { headers: { Accept: 'application/json' }, signal: ctl ? ctl.signal : undefined })
      .then(function (r) {
        if (r.status === 404) { return { status: 'missing', fn: null, id: key, url: url }; }
        if (!r.ok) { throw new Error('MaizeGDB answered HTTP ' + r.status); }
        return r.json().then(function (j) {
          var fn = j && j.data && j.data.sections && j.data.sections['function'];
          if (!fn) { throw new Error('MaizeGDB sent no function section'); }
          return { status: 'ok', fn: fn, id: (j.data.id || key), url: url };
        });
      })
      .catch(function (e) {
        delete cache[key];                       // let a later visit retry
        if (e && e.name === 'AbortError') { throw new Error('MaizeGDB did not answer within ' + (TIMEOUT_MS / 1000) + ' s'); }
        throw e;
      })
      .finally(function () { if (timer) { clearTimeout(timer); } });
    cache[key] = p;
    return p;
  }

  /* ---------------- helpers ---------------- */
  function esc(value) {
    return String(value == null ? '' : value).replace(/[&<>"']/g, function (c) {
      return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c];
    });
  }
  function num(v) { return v == null ? '—' : Number(v).toLocaleString(); }
  function html(tag, cls, inner) {
    var n = document.createElement(tag);
    if (cls) { n.className = cls; }
    if (inner != null) { n.innerHTML = inner; }
    return n;
  }
  function el(name, attrs, text) {
    var n = document.createElementNS(SVG_NS, name);
    Object.keys(attrs || {}).forEach(function (k) { n.setAttribute(k, attrs[k]); });
    if (text != null) { n.textContent = text; }
    return n;
  }
  /* The explorer's equations carry HTML entities (&rarr;, &beta;); decode
     them as text, never as markup. */
  function decodeEntities(text) {
    if (!text) { return ''; }
    var t = document.createElement('textarea');
    t.innerHTML = String(text).replace(/<[^>]*>/g, '');
    return t.value;
  }
  function shortEc(ec) { return ec ? String(ec).replace(/^EC-/, '') : ''; }
  function shortReaction(r) { return r ? String(r).replace(/-RXN$/, '').replace(/^RXN-?/, 'RXN ') : ''; }
  function evidenceKind(code) { return code ? (EVIDENCE_KIND[String(code).toUpperCase()] || 'unknown') : 'unknown'; }
  function truncate(text, max) {
    text = String(text || '');
    return text.length > max ? text.slice(0, max - 1).replace(/\s+\S*$/, '') + '…' : text;
  }
  /* Zm-B73-REFERENCE-NAM-5.0 -> B73 v5; Zm-B97-REFERENCE-NAM-1.0 -> B97. */
  function shortGenome(g) {
    if (/^Zm-B73-REFERENCE-NAM-5\.0$/.test(g || '')) { return 'B73 v5'; }
    var m = /^Zm-([^-]+)-REFERENCE-NAM/.exec(g || '');
    return m ? m[1] : String(g || '').replace(/-(DRAFT|REFERENCE)-PanAnd.*$/, '');
  }
  function debounce(fn, ms) {
    var t = null;
    return function () { clearTimeout(t); t = setTimeout(fn, ms); };
  }
  /* MaizeGDB links that leave SNPTools: absolute, new tab. */
  function outLink(href, label, title, cls) {
    return '<a' + (cls ? ' class="' + cls + '"' : '') + ' href="' + esc(href) + '" target="_blank" rel="noopener"' +
      (title ? ' title="' + esc(title) + '"' : '') + '>' + label + '</a>';
  }

  /* ---- the label fan: measured, not guessed ----
     A name rotated 45 degrees about the bottom-right corner of its box,
     that corner pinned above the centre of its square, reaches
     width / sqrt 2 to the left of the square and stands
     (width + font size) / sqrt 2 tall. The strip's top padding takes the
     tallest name; its left inset takes the furthest reach of the first
     names past the left edge, which depends on the pitch the width
     allows -- solved directly for the pitch the width gives, or for the
     smallest pitch the strip accepts when it will scroll instead. One
     inset for every strip keeps the three rows of squares on one left
     edge. (Constants sit above the functions that read them: a hoisted
     var read before its line is undefined and the inset silently NaN.) */
  var MIN_PITCH = 22;
  var measureCtx;
  function textWidth(font, text) {
    if (measureCtx === undefined) {
      var canvas = document.createElement('canvas');
      measureCtx = (canvas.getContext && canvas.getContext('2d')) || null;
    }
    if (!measureCtx) { return String(text).length * 7; }
    measureCtx.font = font;
    return measureCtx.measureText(text).width;
  }
  function fitStrips(strips) {
    if (!strips.length) { return; }
    var inset = 0;
    var heights = strips.map(function (sp) {
      var lbl = sp.grid.querySelector('.gf-lbl');
      var cs = lbl ? window.getComputedStyle(lbl) : null;
      var fontPx = (cs && parseFloat(cs.fontSize)) || 12;
      var font = '700 ' + fontPx + 'px ' + ((cs && cs.fontFamily) || 'sans-serif');
      var n = sp.names.length || 1;
      var W = sp.scroll.clientWidth;
      var maxW = 0, needAtWidth = 0, needAtMin = 0;
      sp.names.forEach(function (name, i) {
        var w = textWidth(font, name);
        if (w > maxW) { maxW = w; }
        var reach = Math.SQRT1_2 * w;
        var frac = (i + 0.5) / n;
        if (W > 0) { needAtWidth = Math.max(needAtWidth, (reach - W * frac) / (1 - frac)); }
        needAtMin = Math.max(needAtMin, reach - MIN_PITCH * (i + 0.5));
      });
      var need = (W > 0 && (W - needAtWidth) / n >= MIN_PITCH) ? needAtWidth : needAtMin;
      if (need > inset) { inset = need; }
      return Math.ceil(Math.SQRT1_2 * (maxW + fontPx)) + 8;
    });
    inset = Math.ceil(inset);
    strips.forEach(function (sp, i) {
      sp.grid.style.setProperty('--gf-lbl-h', heights[i] + 'px');
      sp.grid.style.setProperty('--gf-inset', inset + 'px');
    });
  }

  /* ======================================================================
     render
     ====================================================================== */
  function render(container, spec) {
    spec = spec || {};
    var fn = spec.fn || {};
    var gene = spec.gene || '';
    var base = spec.base != null ? spec.base : BASE;
    var abs = function (u) { return (u && u.charAt(0) === '/') ? base + u : u; };
    var go = fn.go && fn.go.available ? fn.go : null;
    var pathways = fn.pathways && fn.pathways.available ? fn.pathways : null;
    /* EnTAP's functional annotation (api: function.kegg): EggNOG description,
       COG categories, KEGG orthologs and the KEGG maps kept for maize. */
    var kegg = fn.kegg && fn.kegg.available ? fn.kegg : null;

    var hasGo = !!(go && ((go.terms && go.terms.length) || (go.implied && go.implied.length)));
    var hasPathways = !!(pathways && pathways.pathways && pathways.pathways.length);
    var hasKegg = !!(kegg && (kegg.description || (kegg.kegg_orthologs || []).length || (kegg.kegg_pathways || []).length));
    var nMaps = kegg ? (kegg.kegg_pathways || []).length : 0;
    var nKo = kegg ? (kegg.kegg_orthologs || []).length : 0;

    container.innerHTML = '';
    container.classList.add('gf');
    var tip = html('div', 'gf-tip');
    tip.hidden = true;
    container.appendChild(tip);
    var observers = [];

    /* ---- tooltip ---- */
    function attachTip(node, build) {
      function show(e) { tip.innerHTML = build(); tip.hidden = false; place(e); }
      function place(e) {
        var rect = container.getBoundingClientRect();
        var x = (e && e.clientX != null ? e.clientX : rect.left + 40) - rect.left + 14;
        var y = (e && e.clientY != null ? e.clientY : rect.top + 40) - rect.top + 14;
        if (x + 340 > rect.width) { x = Math.max(0, x - 360); }
        tip.style.left = x + 'px';
        tip.style.top = y + 'px';
      }
      node.addEventListener('mouseenter', show);
      node.addEventListener('mousemove', place);
      node.addEventListener('mouseleave', function () { tip.hidden = true; });
      node.addEventListener('focus', function () { show(null); });
      node.addEventListener('blur', function () { tip.hidden = true; });
    }

    /* ---- cross-highlighting between squares, term rows and graph nodes ---- */
    var pinned = null;
    function clearHot() {
      Array.prototype.forEach.call(container.querySelectorAll('.is-hot'), function (n) { n.classList.remove('is-hot'); });
    }
    function setHot(termIds, slimIds) {
      clearHot();
      (termIds || []).forEach(function (id) {
        Array.prototype.forEach.call(container.querySelectorAll('[data-term="' + id + '"]'), function (n) { n.classList.add('is-hot'); });
      });
      (slimIds || []).forEach(function (id) {
        Array.prototype.forEach.call(container.querySelectorAll('[data-slim="' + id + '"]'), function (n) { n.classList.add('is-hot'); });
      });
    }
    function hotOn(node, key, getTerms, getSlims) {
      node.addEventListener('mouseenter', function () { if (!pinned) { setHot(getTerms(), getSlims()); } });
      node.addEventListener('mouseleave', function () { if (!pinned) { clearHot(); } });
      node.addEventListener('click', function () {
        if (pinned === key) { pinned = null; clearHot(); }
        else { pinned = key; setHot(getTerms(), getSlims()); }
        Array.prototype.forEach.call(container.querySelectorAll('[aria-pressed]'), function (b) {
          if (b.hasAttribute('data-hot-key')) { b.setAttribute('aria-pressed', b.getAttribute('data-hot-key') === pinned ? 'true' : 'false'); }
        });
      });
    }

    /* ---- the numbers first ---- */
    var tiles = html('div', 'gf-tiles');
    var aspectCounts = go ? (go.aspects || {}) : {};
    var litSlim = go ? (go.slim || []).filter(function (s) { return s.terms.length; }).length : 0;
    var totalSlim = go ? (go.slim || []).length : 0;
    var evidenceKinds = {};
    if (go) {
      go.terms.forEach(function (t) {
        var kinds = {};
        (t.evidence || []).forEach(function (c) { kinds[evidenceKind(c)] = true; });
        if (!(t.evidence || []).length) { kinds.unknown = true; }
        Object.keys(kinds).forEach(function (k) { evidenceKinds[k] = (evidenceKinds[k] || 0) + 1; });
      });
    }
    var evNote = Object.keys(evidenceKinds).sort().map(function (k) { return evidenceKinds[k] + ' ' + EVIDENCE_LABEL[k]; }).join(' · ');

    tiles.appendChild(html('div', 'gf-tile',
      '<span class="gf-tile-label">GO terms</span>' +
      '<span class="gf-tile-value">' + (go ? go.terms.length : '—') + '</span>' +
      '<span class="gf-tile-note">' + (go
        ? ASPECTS.map(function (a) { return '<strong>' + (aspectCounts[a.key] || 0) + '</strong> ' + a.short; }).join(' · ') +
          (go.implied.length ? ' · <strong>' + go.implied.length + '</strong> suggested by domains' : '')
        : 'GO reference index not on file') + '</span>'));
    tiles.appendChild(html('div', 'gf-tile',
      '<span class="gf-tile-label">Ontology footprint</span>' +
      '<span class="gf-tile-value">' + (go ? litSlim + ' <small>of ' + totalSlim + '</small>' : '—') + '</span>' +
      '<span class="gf-tile-note">plant GO-slim categories the terms fall under' + (evNote ? ' · ' + esc(evNote) : '') + '</span>'));
    tiles.appendChild(html('div', 'gf-tile',
      '<span class="gf-tile-label">Metabolic pathways</span>' +
      '<span class="gf-tile-value">' + (pathways ? pathways.pathways.length : '—') + '</span>' +
      '<span class="gf-tile-note">' + (pathways
        ? (pathways.pathways.length
            ? '<strong>' + pathways.counts.core + '</strong> core across the NAM founders · ' + pathways.counts.reactions + ' reaction' + (pathways.counts.reactions === 1 ? '' : 's')
            : 'no pathway assignment in the explorer')
        : 'pathway explorer not on file') + '</span>'));
    tiles.appendChild(html('div', 'gf-tile',
      '<span class="gf-tile-label">KEGG maps</span>' +
      '<span class="gf-tile-value">' + (kegg ? nMaps : '—') + '</span>' +
      '<span class="gf-tile-note">' + (kegg
        ? '<strong>' + nKo + '</strong> KEGG ortholog' + (nKo === 1 ? '' : 's') + ' · EnTAP on ' + esc(shortGenome(kegg.genome))
        : 'EnTAP annotation not on file at this MaizeGDB host') + '</span>'));
    container.appendChild(tiles);

    /* ---- the view switch: Gene Ontology | Pathways ---- */
    var uid = 'gf' + Math.random().toString(36).slice(2, 8);
    var bar = html('div', 'gf-switchbar');
    var seg = html('div', 'seg gf-switch');
    seg.setAttribute('role', 'tablist');
    seg.setAttribute('aria-label', 'Function views');
    bar.appendChild(seg);
    container.appendChild(bar);
    var panes = {}, tabs = {};
    var VIEWS = [
      { key: 'go', label: 'Gene Ontology', count: go ? go.terms.length + (go.implied.length ? ' + ' + go.implied.length : '') : '—' },
      { key: 'pathways', label: 'Pathways', count: (pathways ? pathways.pathways.length : 0) + (kegg ? ' + ' + nMaps + ' KEGG' : '') }
    ];
    VIEWS.forEach(function (v) {
      var b = html('button', 'seg-b gf-tab', esc(v.label) + ' <span class="gf-tab-n">' + esc(v.count) + '</span>');
      b.type = 'button';
      b.id = uid + '-tab-' + v.key;
      b.setAttribute('role', 'tab');
      b.setAttribute('aria-controls', uid + '-pane-' + v.key);
      b.addEventListener('click', function () { setView(v.key, true); });
      b.addEventListener('keydown', function (e) {
        if (e.key !== 'ArrowRight' && e.key !== 'ArrowLeft') { return; }
        e.preventDefault();
        var next = VIEWS[(VIEWS.indexOf(v) + (e.key === 'ArrowRight' ? 1 : VIEWS.length - 1)) % VIEWS.length].key;
        setView(next, true);
        tabs[next].focus();
      });
      seg.appendChild(b);
      tabs[v.key] = b;
      var pane = html('div', 'gf-pane gf-pane-' + v.key);
      pane.id = uid + '-pane-' + v.key;
      pane.setAttribute('role', 'tabpanel');
      pane.setAttribute('aria-labelledby', b.id);
      panes[v.key] = pane;
      container.appendChild(pane);
    });
    var current = null;
    function setView(key, user) {
      if (!panes[key]) { key = 'go'; }
      current = key;
      VIEWS.forEach(function (v) {
        var on = v.key === key;
        tabs[v.key].classList.toggle('on', on);
        tabs[v.key].setAttribute('aria-selected', on ? 'true' : 'false');
        tabs[v.key].tabIndex = on ? 0 : -1;
        panes[v.key].hidden = !on;
      });
      tip.hidden = true;
      if (user && typeof spec.onView === 'function') { spec.onView(key); }
    }

    /* ==================================================================
       Gene Ontology: fingerprint, terms, ancestry
       ================================================================== */
    var slimById = {};
    var termById = {};
    if (hasGo) {
      go.slim.forEach(function (s) { slimById[s.id] = s; });
      go.terms.forEach(function (t) { termById[t.term] = t; });
      go.implied.forEach(function (t) { if (!termById[t.term]) { termById[t.term] = t; } });

      var block = html('div', 'gf-block gf-go');
      block.appendChild(html('div', 'gf-block-head',
        '<h4>Gene Ontology <small>' + go.terms.length + ' term' + (go.terms.length === 1 ? '' : 's') + ' · GO release ' + esc((go.release || '').replace(/^releases\//, '')) + ' · from MaizeGDB</small></h4>'));

      /* one chip shows or hides the three ancestry graphs together */
      var graphs = [];
      var hasGraph = !!(go.graph && go.graph.nodes && go.graph.nodes.length);
      if (hasGraph) {
        var toggle = html('button', 'gf-chip', 'Hide ancestry');
        toggle.type = 'button';
        toggle.setAttribute('aria-pressed', 'true');
        block.firstChild.appendChild(toggle);
        toggle.addEventListener('click', function () {
          var open = !!(graphs.length && graphs[0].stage.hidden);
          graphs.forEach(function (g) { g.stage.hidden = !open; if (open) { drawGraph(g.stage, go, g.aspect.key); } });
          toggle.textContent = open ? 'Hide ancestry' : 'Show ancestry';
          toggle.setAttribute('aria-pressed', open ? 'true' : 'false');
        });
      }

      /* one full-width row per aspect: its ancestry graph, then the
         plant-slim categories named above their squares at 45 degrees,
         the squares in one unbroken strip across the width (a fixed order,
         so the pattern is comparable between genes), and the terms
         beneath */
      var rows = html('div', 'gf-aspect-rows');
      rows.appendChild(html('div', 'gf-legend',
        '<span><i class="gf-sq is-lit gf-sq-bp"></i> category with a term</span>' +
        '<span><i class="gf-sq is-implied gf-sq-bp"></i> suggested by domains only</span>' +
        '<span><i class="gf-sq gf-sq-bp"></i> not touched</span>' +
        '<span>ancestry: each term traced to its root through the plant-slim categories named beneath it</span>' +
        '<span>hover a square, a name, a term or a node to see what belongs together; click to pin</span>'));
      var strips = [];
      ASPECTS.forEach(function (a) {
        var terms = go.terms.filter(function (t) { return t.aspect === a.key; });
        var implied = go.implied.filter(function (t) { return t.aspect === a.key; });
        var slims = go.slim.filter(function (s) { return s.aspect === a.key; });
        var lit = slims.filter(function (s) { return s.terms.length; });
        var row = html('section', 'gf-aspect gf-aspect-' + a.cls);
        row.setAttribute('aria-label', a.label);
        row.appendChild(html('h5', null, esc(a.label) + ' <small>' + terms.length + ' term' + (terms.length === 1 ? '' : 's') +
          (implied.length ? ' · ' + implied.length + ' suggested by domains' : '') +
          ' · ' + lit.length + ' of ' + slims.length + ' plant-slim categories</small>'));

        /* this aspect's ancestry, drawn once the row is in the document */
        if (hasGraph && go.graph.nodes.some(function (n) { return n.namespace === a.key; })) {
          var stage = html('div', 'gf-graph gf-graph-' + a.cls);
          stage.setAttribute('role', 'img');
          stage.setAttribute('aria-label', a.label + ' ancestry of the terms of this gene');
          row.appendChild(stage);
          graphs.push({ stage: stage, aspect: a });
        }

        var scroll = html('div', 'gf-strip-scroll');
        var grid = html('div', 'gf-fpx');
        grid.setAttribute('role', 'group');
        grid.setAttribute('aria-label', a.label + ' plant GO-slim categories');
        grid.style.setProperty('--gf-n', String(slims.length || 1));
        slims.forEach(function (s, i) {
          var n = s.terms.length;
          var cell = html('div', 'gf-cell');
          var title = s.name + (n ? ': ' + n + ' term' + (n === 1 ? '' : 's') : (s.implied.length ? ': suggested by domains' : ''));
          var tipFor = function () {
            return '<strong>' + esc(s.name) + '</strong><span class="gf-tip-muted">' + esc(s.id) + ' · plant slim, ' + esc(a.label.toLowerCase()) + '</span>' +
              (n ? '<ul>' + s.terms.map(function (id) { return '<li>' + esc(termById[id] ? termById[id].name : id) + '</li>'; }).join('') + '</ul>' : '') +
              (s.implied.length ? '<div class="gf-tip-muted">suggested by domains: ' + s.implied.map(function (id) { return esc(termById[id] ? termById[id].name : id); }).join(', ') + '</div>' : '') +
              (!n && !s.implied.length ? '<div class="gf-tip-muted">no term of this gene falls here</div>' : '');
          };
          /* the name, rotated above the square; the button is the control
             and carries the same name for assistive technology */
          var lbl = html('span', 'gf-lbl' + (n ? ' is-lit' : (s.implied.length ? ' is-implied' : '')), esc(s.name));
          lbl.setAttribute('data-slim', s.id);
          lbl.setAttribute('aria-hidden', 'true');
          var b = html('button', 'gf-sq gf-sq-' + a.cls + (n ? ' is-lit is-lit-' + Math.min(n, 3) : (s.implied.length ? ' is-implied' : '')));
          b.type = 'button';
          b.setAttribute('data-slim', s.id);
          b.setAttribute('data-hot-key', 'slim:' + s.id);
          b.setAttribute('aria-pressed', 'false');
          b.style.transitionDelay = (i * 18) + 'ms';
          b.title = title;
          b.setAttribute('aria-label', title);
          [lbl, b].forEach(function (node) {
            hotOn(node, 'slim:' + s.id, function () { return s.terms.concat(s.implied); }, function () { return [s.id]; });
            attachTip(node, tipFor);
          });
          cell.appendChild(lbl);
          cell.appendChild(b);
          grid.appendChild(cell);
        });
        scroll.appendChild(grid);
        row.appendChild(scroll);
        strips.push({ grid: grid, scroll: scroll, names: slims.map(function (s) { return s.name; }) });

        /* the terms themselves, beneath the strip */
        if (!terms.length && !implied.length) {
          row.appendChild(html('p', 'gf-note', 'No term of this gene in this aspect.'));
        } else {
          if (terms.length && !lit.length) {
            row.appendChild(html('p', 'gf-note', 'These terms sit outside the plant slim, so no category lights.'));
          }
          var list = html('ul', 'gf-terms');
          terms.forEach(function (t) { list.appendChild(termRow(t, a, false)); });
          implied.forEach(function (t) { list.appendChild(termRow(t, a, true)); });
          row.appendChild(list);
        }
        rows.appendChild(row);
      });
      block.appendChild(rows);

      var unplaced = go.terms.filter(function (t) { return !t.aspect; });
      if (unplaced.length) {
        block.appendChild(html('p', 'gf-note', unplaced.length + ' term' + (unplaced.length === 1 ? '' : 's') +
          ' the GO release no longer carries (' + unplaced.map(function (t) { return esc(t.term); }).join(', ') + '): not placed here.'));
      }
      panes.go.appendChild(block);
    } else {
      panes.go.appendChild(html('p', 'gf-note gf-empty', go
        ? 'MaizeGDB records no Gene Ontology term for ' + esc(gene) + ', and no InterPro entry on its protein implies one.'
        : 'MaizeGDB could not place GO terms here: its GO reference index is not on file at this host.'));
    }

    function termRow(t, a, implied) {
      var li = html('li', 'gf-term' + (implied ? ' gf-term-implied' : '') + (t.obsolete ? ' gf-term-obsolete' : ''));
      li.setAttribute('data-term', t.term);
      var badges = '';
      if (implied) {
        badges += '<span class="gf-badge gf-badge-implied" title="Not an annotation of this gene: InterPro2GO maps a domain on the protein to this term">domain-implied</span>';
      } else {
        var kinds = {};
        (t.evidence || []).forEach(function (c) { kinds[evidenceKind(c)] = c; });
        var keys = Object.keys(kinds);
        if (!keys.length) { keys = ['unknown']; }
        keys.forEach(function (k) {
          badges += '<span class="gf-badge gf-badge-' + k + '" title="evidence ' + esc(kinds[k] || 'not stated') + '">' + esc(EVIDENCE_LABEL[k]) + '</span>';
        });
        if (t.implied_by && t.implied_by.length) {
          badges += '<span class="gf-badge gf-badge-implied" title="also implied by ' + esc(t.implied_by.map(function (x) { return x.accession; }).join(', ')) + '">◆ domain</span>';
        }
      }
      if (t.obsolete) { badges += '<span class="gf-badge gf-badge-obsolete">retired' + (t.replaced_by ? ' → ' + esc(t.replaced_by) : '') + '</span>'; }
      if (t.merged_into) { badges += '<span class="gf-badge gf-badge-merged" title="The annotation carries the old id; the GO release merged it into this term">merged into ' + esc(t.merged_into) + '</span>'; }
      li.innerHTML = '<span class="gf-term-name">' + esc(t.name || t.db_name || t.term) + '</span>' +
        '<span class="gf-term-meta">' + outLink(t.url, esc(t.term)) + badges +
        (t.pmids || []).map(function (id) { return outLink('https://pubmed.ncbi.nlm.nih.gov/' + encodeURIComponent(id) + '/', 'PMID ' + esc(id)); }).join('') + '</span>';
      hotOn(li, 'term:' + t.term, function () { return [t.term]; }, function () { return (t.slim_ancestors || []).map(function (s) { return s.id; }); });
      attachTip(li, function () {
        var path = (t.slim_ancestors || []).map(function (s) { return esc(s.name); });
        return '<strong>' + esc(t.name || t.db_name || t.term) + '</strong>' +
          '<span class="gf-tip-muted">' + esc(t.term) + (a ? ' · ' + esc(a.label.toLowerCase()) : '') + (t.depth != null ? ' · depth ' + t.depth : '') + '</span>' +
          (t.definition ? '<p>' + esc(t.definition) + '</p>' : '') +
          (path.length ? '<div class="gf-tip-path">' + path.join(' › ') + '</div>' : '') +
          (implied
            ? '<div class="gf-tip-muted">suggested by ' + (t.implied_by || []).map(function (x) { return esc(x.accession) + ' ' + esc(x.name); }).join('; ') + ' (InterPro2GO)</div>'
            : '<div class="gf-tip-muted">' +
                ((t.evidence || []).length ? 'evidence ' + esc(t.evidence.join(', ')) + ' · ' : '') +
                ((t.sources || []).length ? 'source ' + esc(t.sources.join(', ')) : 'source not stated') +
                ((t.proteins || []).length ? ' · on ' + esc(t.proteins.join(', ')) : '') +
                ((t.comments || []).length ? ' · ' + esc(t.comments.join('; ')) : '') +
              '</div>');
      });
      return li;
    }

    /* ---- the ancestry graph: one aspect per card, layered by depth,
       ordered by barycentre, edges child -> parent. A sparse graph is drawn
       compact and centred rather than flung to the corners of a wide card. ---- */
    function drawGraph(stage, go, only) {
      stage.innerHTML = '';
      var nodes = go.graph.nodes;
      var edges = go.graph.edges;
      var byNs = {};
      nodes.forEach(function (n) { if (n.namespace && (!only || n.namespace === only)) { (byNs[n.namespace] = byNs[n.namespace] || []).push(n); } });
      var colsPresent = ASPECTS.filter(function (a) { return byNs[a.key] && byNs[a.key].length; });
      if (!colsPresent.length) { return; }
      var width = Math.max(320, (stage.clientWidth || 900) - 16);
      var single = colsPresent.length === 1;
      var narrow = !single && width < 700;
      var colW = narrow ? width : Math.floor(width / colsPresent.length);
      var titleH = single ? 0 : 18;
      var nodeMax = single ? 300 : 210;
      var rowH = 44, nodeH = 22, padX = 10, padY = 14;

      var parentsOf = {}, childrenOf = {};
      edges.forEach(function (e) {
        (parentsOf[e[0]] = parentsOf[e[0]] || []).push(e[1]);
        (childrenOf[e[1]] = childrenOf[e[1]] || []).push(e[0]);
      });

      var layouts = [];
      var maxRows = 0;
      colsPresent.forEach(function (a) {
        var list = byNs[a.key];
        var maxDepth = 0;
        list.forEach(function (n) { if (n.depth != null && n.depth > maxDepth) { maxDepth = n.depth; } });
        var layers = [];
        for (var d = 0; d <= maxDepth + 1; d++) { layers.push([]); }
        list.forEach(function (n) { layers[n.depth == null ? maxDepth + 1 : n.depth].push(n); });
        if (!layers[maxDepth + 1].length) { layers.pop(); }
        layers.forEach(function (layer) { layer.sort(function (x, y) { return x.name < y.name ? -1 : 1; }); });
        var pos = {};
        function assign() { layers.forEach(function (layer) { layer.forEach(function (n, i) { pos[n.id] = i; }); }); }
        assign();
        function bary(n, rel) {
          var ids = (rel[n.id] || []).filter(function (id) { return pos[id] != null; });
          if (!ids.length) { return pos[n.id]; }
          return ids.reduce(function (s, id) { return s + pos[id]; }, 0) / ids.length;
        }
        for (var pass = 0; pass < 4; pass++) {
          for (var d2 = 1; d2 < layers.length; d2++) {
            layers[d2].sort(function (x, y) { return bary(x, parentsOf) - bary(y, parentsOf) || (x.name < y.name ? -1 : 1); });
            assign();
          }
          for (var d3 = layers.length - 2; d3 >= 0; d3--) {
            layers[d3].sort(function (x, y) { return bary(x, childrenOf) - bary(y, childrenOf) || (x.name < y.name ? -1 : 1); });
            assign();
          }
        }
        layouts.push({ aspect: a, layers: layers });
        if (layers.length > maxRows) { maxRows = layers.length; }
      });

      var height = narrow ? layouts.reduce(function (h, L) { return h + L.layers.length * rowH + padY * 2 + titleH; }, 0)
                          : maxRows * rowH + padY * 2 + titleH;
      var svg = el('svg', { viewBox: '0 0 ' + width + ' ' + height, width: '100%', height: height, 'class': 'gf-graph-svg', role: 'img',
                            'aria-label': 'Ancestry of the GO terms of this gene' });
      var coords = {};
      var yOffset = 0;
      layouts.forEach(function (L, ci) {
        var x0 = narrow ? 0 : ci * colW;
        var y0 = narrow ? yOffset : 0;
        if (!single) {
          svg.appendChild(el('text', { x: x0 + padX, y: y0 + 12, 'class': 'gf-graph-title gf-graph-title-' + L.aspect.cls }, L.aspect.label));
        }
        if (!narrow && ci > 0) {
          svg.appendChild(el('line', { x1: x0, x2: x0, y1: 0, y2: height, 'class': 'gf-graph-divider' }));
        }
        var widest = L.layers.reduce(function (m, layer) { return layer.length > m ? layer.length : m; }, 1);
        var used = single ? Math.min(colW, Math.max(440, widest * 250)) : colW;
        var xBase = x0 + Math.floor((colW - used) / 2);
        L.layers.forEach(function (layer, d) {
          var slot = (used - padX * 2) / layer.length;
          var w = Math.min(nodeMax, Math.max(56, slot - 8));
          layer.forEach(function (n, i) {
            var cx = xBase + padX + slot * (i + 0.5);
            var cy = y0 + titleH + padY + d * rowH + rowH / 2;
            coords[n.id] = { x: cx, y: cy, w: w, h: nodeH };
          });
        });
        yOffset += L.layers.length * rowH + padY * 2 + titleH;
      });

      var gEdges = el('g', { 'class': 'gf-graph-edges' });
      edges.forEach(function (e) {
        var c = coords[e[0]], p = coords[e[1]];
        if (!c || !p) { return; }
        var y1 = c.y - c.h / 2, y2 = p.y + p.h / 2;
        var my = (y1 + y2) / 2;
        gEdges.appendChild(el('path', { d: 'M' + c.x + ',' + y1 + ' C' + c.x + ',' + my + ' ' + p.x + ',' + my + ' ' + p.x + ',' + y2,
                                        'class': 'gf-graph-edge' }));
      });
      svg.appendChild(gEdges);

      nodes.forEach(function (n) {
        var c = coords[n.id];
        if (!c) { return; }
        var a = ASPECT_BY_KEY[n.namespace] || ASPECTS[0];
        var g = el('g', { 'class': 'gf-node gf-node-' + n.kind + ' gf-node-' + a.cls, tabindex: 0, role: 'img' });
        if (n.annotated) { g.setAttribute('data-term', n.id); }
        if (n.slim || n.kind === 'root') { g.setAttribute('data-slim', n.id); }
        g.appendChild(el('rect', { x: c.x - c.w / 2, y: c.y - c.h / 2, width: c.w, height: c.h, rx: 11 }));
        var maxChars = Math.max(6, Math.floor((c.w - 12) / 5.9));
        g.appendChild(el('text', { x: c.x, y: c.y + 4, 'text-anchor': 'middle' }, truncate(n.name, maxChars)));
        var title = n.name + ' (' + n.id + ')' + (n.kind === 'root' ? ', root' : n.kind === 'slim' ? ', plant-slim category' : ', annotated term');
        g.appendChild(el('title', {}, title));
        g.setAttribute('aria-label', title);
        var t = termById[n.id];
        hotOn(g, 'node:' + n.id,
          function () { return n.annotated ? [n.id] : (slimById[n.id] ? slimById[n.id].terms.concat(slimById[n.id].implied) : []); },
          function () { return n.annotated && t ? [n.id].concat((t.slim_ancestors || []).map(function (s) { return s.id; })) : [n.id]; });
        attachTip(g, function () {
          return '<strong>' + esc(n.name) + '</strong><span class="gf-tip-muted">' + esc(n.id) + ' · ' +
            (n.kind === 'root' ? 'root of ' + esc(a.label.toLowerCase()) : n.kind === 'slim' ? 'plant-slim category' : 'annotated term') +
            (n.depth != null ? ' · depth ' + n.depth : '') + '</span>' +
            (t && t.definition ? '<p>' + esc(t.definition) + '</p>' : '');
        });
        svg.appendChild(g);
      });
      stage.appendChild(svg);
    }

    /* ------------------------------------------------------------------
       Genome groups (js/mgdb-genome-groups.js, vendored from MaizeGDB):
       the pathway presence dots color each NAM founder by its heterotic
       group or population and carry a key naming the groups drawn.
       ------------------------------------------------------------------ */
    var GG = (window.MGDB && window.MGDB.genomeGroups) || null;
    var groupsSeen = [];
    function groupStyle(g) {
      return '--mgdb-genome-fill:' + g.color + ';--mgdb-genome-ink:' + g.ink + ';--mgdb-genome-edge:' + g.edge;
    }
    /* Lineage order, then label; a copy, so the API's arrays are untouched. */
    function lineageOrder(list, idOf) {
      if (!GG) { return list; }
      return list.slice().sort(function (a, b) {
        return (GG.rank(idOf(a)) - GG.rank(idOf(b))) ||
          String(idOf(a)).localeCompare(String(idOf(b)), undefined, { numeric: true });
      });
    }
    function groupKey(label) {
      var keyHtml = GG && groupsSeen.length ? GG.legendHtml(groupsSeen, { label: label }) : '';
      groupsSeen = [];
      return keyHtml ? html('div', 'gf-groups', keyHtml) : null;
    }

    /* ==================================================================
       Metabolic pathways (pathway explorer, E2P2 on the NAM founders)
       ================================================================== */
    if (hasPathways) {
      var pb = html('div', 'gf-block gf-pathways');
      var founders = pathways.pathways[0].founders || 26;
      pb.appendChild(html('div', 'gf-block-head',
        '<h4>Metabolic pathways <small>' + pathways.pathways.length + ' pathway' + (pathways.pathways.length === 1 ? '' : 's') +
        ' · E2P2 assignment on ' + esc(pathways.genome_label || pathways.genome) + ', compared across ' + founders + ' NAM founders</small></h4>' +
        outLink(abs(pathways.explorer), 'Pathway explorer ↗', 'Open the MaizeGDB pan-genome pathway explorer', 'gf-out') +
        '<div class="gf-legend">' +
          '<span><i class="gf-step-key is-this"></i> this gene’s step</span>' +
          '<span><i class="gf-step-key is-filled"></i> another ' + esc(pathways.genome_label || 'B73') + ' gene fills it</span>' +
          '<span><i class="gf-step-key is-empty"></i> no gene in ' + esc(pathways.genome_label || 'B73') + '</span>' +
          '<span><i class="gf-step-key is-bar"></i> founders with a gene</span>' +
        '</div>'));
      var pwList = html('div', 'gf-pw-list');
      groupsSeen = [];
      pathways.pathways.forEach(function (p) { pwList.appendChild(pathwayCard(p, pathways)); });
      /* The presence dots' group key, under the marks key. */
      var pwKey = groupKey('Presence dots: heterotic group or population');
      if (pwKey) { pb.appendChild(pwKey); }
      pb.appendChild(pwList);
      if (pathways.counts.truncated) {
        pb.appendChild(html('p', 'gf-note', 'Steps are drawn for the first ' + pathways.counts.files_read + ' pathways; the rest are listed by name. The explorer shows every one.'));
      }
      panes.pathways.appendChild(pb);
    } else if (pathways) {
      /* the explorer answers genome: null for a gene it never assigns, so
         say so whenever it answered at all, not only when it named a genome */
      panes.pathways.appendChild(html('div', 'gf-block gf-pathways',
        '<div class="gf-block-head"><h4>Metabolic pathways <small>E2P2 on the NAM founders</small></h4>' +
        outLink(abs(pathways.explorer), 'Pathway explorer ↗', null, 'gf-out') + '</div>' +
        '<p class="gf-note gf-empty">The pathway explorer assigns no reaction step to this gene.</p>'));
    }

    /* ==================================================================
       KEGG and EggNOG (EnTAP)
       ================================================================== */
    if (hasKegg) {
      var kb = html('div', 'gf-block gf-kegg');
      kb.appendChild(html('div', 'gf-block-head',
        '<h4>KEGG and EggNOG <small>EnTAP on ' + esc(shortGenome(kegg.genome)) + ' · ' + nMaps + ' KEGG map' + (nMaps === 1 ? '' : 's') +
        ' · ' + nKo + ' ortholog' + (nKo === 1 ? '' : 's') + '</small></h4>' +
        outLink(abs(kegg.links.hub), 'Metabolic pathways hub ↗', 'This gene in the MaizeGDB Metabolic Pathways hub', 'gf-out')));

      if (kegg.description || (kegg.cog || []).length) {
        kb.appendChild(html('p', 'gf-kegg-desc',
          (kegg.description ? '<span class="gf-muted">EggNOG</span> ' + esc(kegg.description) : '') +
          (kegg.cog || []).map(function (c) {
            return ' <span class="gf-cog" title="COG functional category">' + esc(c.letter) + (c.name ? ' · ' + esc(c.name) : '') + '</span>';
          }).join('')));
      }
      if (nKo) {
        var kol = html('ul', 'gf-ko-list');
        kegg.kegg_orthologs.forEach(function (k) {
          kol.appendChild(html('li', null,
            outLink(k.links.kegg, esc(k.id), null, 'gf-ko-id') +
            (k.symbol ? ' <strong>' + esc(k.symbol) + '</strong>' : '') +
            (k.name ? ' <span class="gf-muted">' + esc(String(k.name).replace(/\s*\[EC:[^\]]*\]/, '')) + '</span>' : '') +
            ((k.ec || []).length ? ' <span class="gf-ko-ec">EC ' + esc(k.ec.join(', ')) + '</span>' : '')));
        });
        kb.appendChild(kol);
      }
      if (nMaps) {
        var ml = html('div', 'gf-pw-list gf-kegg-maps');
        kegg.kegg_pathways.forEach(function (m) {
          var card = html('div', 'gf-pw');
          card.appendChild(html('div', 'gf-pw-head',
            outLink(abs(m.links.hub), esc(m.name || m.id), 'This map in the MaizeGDB Metabolic Pathways hub', 'gf-pw-name') +
            '<span class="gf-muted">' + esc(m.id) + '</span>' +
            (m.subclass ? '<span class="gf-pw-class">' + esc(m['class'] ? m['class'] + ' › ' + m.subclass : m.subclass) + '</span>' : '')));
          var links = [];
          if (m.links.kegg_highlighted) { links.push(outLink(m.links.kegg_highlighted, 'KEGG map, this gene’s orthologs in red', 'KEGG draws this gene’s orthologs in red')); }
          if (m.links.kegg_maize) { links.push(outLink(m.links.kegg_maize, 'maize map')); }
          card.appendChild(html('div', 'gf-pw-foot',
            (m.genes_in_genome != null ? '<span><strong>' + num(m.genes_in_genome) + '</strong> ' + esc(shortGenome(kegg.genome)) + ' genes on this map</span>' : '') +
            (links.length ? '<span class="gf-pw-links">' + links.join(' · ') + '</span>' : '')));
          ml.appendChild(card);
        });
        kb.appendChild(ml);
      } else {
        kb.appendChild(html('p', 'gf-note', 'EnTAP places this gene on no KEGG map that KEGG lists for maize.'));
      }
      var notes2 = [];
      if ((kegg.excluded_maps || []).length) {
        notes2.push('EggNOG also names ' + kegg.excluded_maps.length + ' map' + (kegg.excluded_maps.length === 1 ? '' : 's') +
          ' KEGG does not list for maize, or its global overview maps (' + esc(kegg.excluded_maps.slice(0, 4).join(', ')) +
          (kegg.excluded_maps.length > 4 ? ', …' : '') + '); they are left out.');
      }
      if (kegg.best_hit && kegg.best_hit.subject) {
        var h = kegg.best_hit;
        notes2.push('Closest protein ' + esc(h.description || h.subject) +
          (h.identity != null ? ', ' + h.identity + '% identity' : '') + (h.coverage != null ? ' over ' + h.coverage + '% coverage' : '') +
          (h.database ? ' (' + esc(h.database.replace(/_/g, ' ')) + ')' : '') + '.');
      }
      if (notes2.length) { kb.appendChild(html('p', 'gf-note', notes2.join(' '))); }
      panes.pathways.appendChild(kb);
    } else if (kegg) {
      panes.pathways.appendChild(html('p', 'gf-note gf-empty', 'EnTAP annotates no KEGG ortholog or map for this gene.'));
    }
    if (!panes.pathways.childNodes.length) {
      panes.pathways.appendChild(html('p', 'gf-note gf-empty', 'MaizeGDB has no pathway data for ' + esc(gene) + ' at this host.'));
    }

    function pathwayCard(p, ctx) {
      var card = html('div', 'gf-pw');
      var mine = p.this_gene || [];
      var head = html('div', 'gf-pw-head',
        outLink(abs(p.url), esc(p.name), 'Open this pathway in the MaizeGDB pathway explorer', 'gf-pw-name') +
        '<span class="gf-pan gf-pan-' + esc(p.pan) + '" title="pan-genome status across the NAM founders">' + esc(p.pan) + '</span>' +
        (p.variability ? '<span class="gf-var" title="completeness variability across founders">' + esc(p.variability) + ' variability</span>' : '') +
        (p.class_tail ? '<span class="gf-pw-class">' + esc(p.class_tail) + '</span>' : ''));
      card.appendChild(head);

      if (p.steps && p.steps.length) {
        var strip = html('div', 'gf-steps');
        strip.setAttribute('role', 'list');
        p.steps.forEach(function (s, i) {
          if (i > 0) { strip.appendChild(html('span', 'gf-step-arrow', '<span aria-hidden="true">&rarr;</span>')); }
          var state = s.this_gene ? 'is-this' : (s.filled ? 'is-filled' : 'is-empty');
          var b = html('button', 'gf-step ' + state);
          b.type = 'button';
          b.setAttribute('role', 'listitem');
          var label = shortEc(s.ec) || shortReaction(s.reaction);
          var share = ctx.pathways[0].founders ? s.founders_filled / ctx.pathways[0].founders : 0;
          b.innerHTML = '<span class="gf-step-ec">' + esc(label) + '</span>' +
            '<span class="gf-step-enz">' + esc(truncate(s.enzyme || s.common_name || s.reaction, 28)) + '</span>' +
            '<span class="gf-step-bar"><span style="width:' + Math.round(share * 100) + '%"></span></span>';
          var geneNames = (s.genes || []).map(function (g) { return g.symbol || g.gene; });
          var title = (s.enzyme || s.reaction) + (s.ec ? ' (' + s.ec + ')' : '') + ': ' +
            (s.this_gene ? 'this gene’s step; ' : '') + (s.filled ? s.genes_total_here + ' gene' + (s.genes_total_here === 1 ? '' : 's') + ' in ' + ctx.genome_label : 'no gene in ' + ctx.genome_label) +
            '; a gene in ' + s.founders_filled + ' of ' + ctx.pathways[0].founders + ' founders';
          b.title = title;
          b.setAttribute('aria-label', title);
          attachTip(b, function () {
            return '<strong>' + esc(s.enzyme || s.common_name || s.reaction) + '</strong>' +
              '<span class="gf-tip-muted">' + esc(s.reaction) + (s.ec ? ' · ' + esc(s.ec) : '') + (s.occurrence ? ' · ' + esc(s.occurrence) + ' step' : '') + '</span>' +
              (s.equation ? '<p class="gf-tip-eq">' + esc(decodeEntities(s.equation)) + '</p>' : '') +
              '<div>' + (s.this_gene ? '<strong>This gene fills this step.</strong> ' : '') +
                (s.filled ? esc(ctx.genome_label) + ' genes: ' + esc(geneNames.join(', ')) + (s.genes_total_here > geneNames.length ? ' +' + (s.genes_total_here - geneNames.length) + ' more' : '')
                          : 'No ' + esc(ctx.genome_label) + ' gene is assigned here' + (s.gap_class ? ' (' + esc(s.gap_class) + ')' : '') + '.') + '</div>' +
              '<div class="gf-tip-muted">a gene in ' + s.founders_filled + ' of ' + ctx.pathways[0].founders + ' NAM founders</div>';
          });
          strip.appendChild(b);
        });
        card.appendChild(strip);
      } else {
        card.appendChild(html('p', 'gf-note', 'Steps not drawn: ' + mine.map(function (m) { return esc(shortReaction(m.reaction)); }).join(', ') + '.'));
      }

      var foot = html('div', 'gf-pw-foot');
      if (p.presence) {
        var dots = html('span', 'gf-presence');
        dots.title = 'pathway present in ' + p.present_in + ' of ' + p.founders + ' NAM founders';
        /* A present founder's dot is filled in its group color, the dots run
           in lineage order, and each names its founder on hover. */
        lineageOrder(p.presence, function (g) { return g.genome; }).forEach(function (g) {
          var grp = GG ? GG.of(g.genome) : null;
          var dot = html('i', 'gf-presence-dot' + (g.present ? ' is-on' : ''), '');
          /* Absent: grey, edged in the group color (css). */
          if (grp) { dot.setAttribute('style', groupStyle(grp)); groupsSeen.push(g.genome); }
          if (GG) { dot.title = g.genome + (grp ? ' (' + grp.label + ')' : '') + ': ' + (g.present ? 'pathway present' : 'absent'); }
          dots.appendChild(dot);
        });
        foot.appendChild(dots);
        foot.appendChild(html('span', null, '<strong>' + p.present_in + '</strong> of ' + p.founders + ' founders have the pathway'));
      }
      if (p.genome && p.genome.steps_filled != null) {
        foot.appendChild(html('span', null, esc(ctx.genome_label) + ': <strong>' + p.genome.steps_filled + '</strong> of ' + (p.reactions != null ? p.reactions : (p.steps ? p.steps.length : '?')) +
          ' steps have a gene' + (p.genome.genes != null ? ' (' + p.genome.genes + ' genes)' : '')));
      }
      if (mine.length) {
        foot.appendChild(html('span', 'gf-muted', 'evidence ' + esc(mine.map(function (m) { return m.evidence || 'unstated'; }).filter(function (v, i, arr) { return arr.indexOf(v) === i; }).join(', '))));
      }
      var links = [];
      if (p.links && p.links.metacyc) { links.push(outLink(p.links.metacyc, 'MetaCyc')); }
      if (p.links && p.links.plantcyc) { links.push(outLink(p.links.plantcyc, 'PlantCyc')); }
      if (links.length) { foot.appendChild(html('span', 'gf-pw-links', links.join(' · '))); }
      card.appendChild(foot);
      return card;
    }

    /* ---- footer: where the numbers come from ---- */
    var footer = html('div', 'gf-footer');
    var notes = [];
    if (go) { notes = notes.concat(go.notes || []); }
    if (pathways) { notes.push(pathways.source + ' A step is "filled" when the explorer assigns any gene of the genome to its reaction.'); }
    if (kegg) { notes.push('KEGG and EggNOG: ' + kegg.source + '; KEGG maps are assigned through orthology and kept to the maps KEGG lists for maize.'); }
    footer.appendChild(html('p', null, esc(notes.join(' '))));
    var flinks = html('div', 'gf-footer-links');
    flinks.innerHTML = outLink(recordUrlAt(base, gene), 'MaizeGDB gene record ↗') + outLink(base + '/api/v1/records/gene/' + encodeURIComponent(gene) + '?fields=function', 'JSON ↗', 'The API answer these views are drawn from');
    footer.appendChild(flinks);
    container.appendChild(footer);

    /* Pick the opening view: the one asked for, else Gene Ontology unless
       there is no GO and there are pathways. */
    var startView = spec.view || ((!hasGo && (hasPathways || hasKegg)) ? 'pathways' : 'go');
    setView(startView, false);

    /* The strips need two measurements the stylesheet cannot make: how
       tall the longest name stands at 45 degrees, and how far the first
       names of a strip reach past its left edge. Both come from the text
       itself, so they are taken now and again whenever the width changes --
       including the change from nothing to something when the GO view (or
       the whole tool) is shown after being hidden. */
    if (hasGo) {
      var refit = function () {
        fitStrips(strips);
        graphs.forEach(function (g) { if (!g.stage.hidden) { drawGraph(g.stage, go, g.aspect.key); } });
      };
      refit();
      var refitLater = debounce(refit, 150);
      if (typeof window.ResizeObserver === 'function') {
        var lastWidth = rows.clientWidth;
        var ro = new window.ResizeObserver(function () {
          if (rows.clientWidth !== lastWidth) { lastWidth = rows.clientWidth; refitLater(); }
        });
        ro.observe(rows);
        observers.push(function () { ro.disconnect(); });
      } else {
        window.addEventListener('resize', refitLater);
        observers.push(function () { window.removeEventListener('resize', refitLater); });
      }
    }

    /* light the fingerprint up once it is on screen (a timer as well as a
       frame: a background tab paints no frames until it is shown) */
    function ready() { container.classList.add('is-ready'); }
    window.requestAnimationFrame(ready);
    window.setTimeout(ready, 60);

    return {
      setView: function (v) { setView(v, false); },
      view: function () { return current; },
      destroy: function () { observers.forEach(function (off) { off(); }); observers = []; }
    };
  }

  function recordUrlAt(base, gene) { return base + '/gene_center/gene/' + encodeURIComponent(gene); }

  window.SNPFunctionOntology = {
    base: BASE,
    fetch: fetchFunction,
    render: render,
    apiUrl: apiUrl,
    recordUrl: recordUrl
  };
}(window, document));
