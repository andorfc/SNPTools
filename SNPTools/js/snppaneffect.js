/* =====================================================================
 *  snppaneffect.js — PanEffect as a native SNPTools panel.
 *
 *  Loads AFTER core.js (needs SNPTools, S, ICONS, go()).  Registers the
 *  'paneffect' tool.  One left-menu item hosts both views (Variant effects
 *  in B73 / Pan-genome) as a Both / B73 / Pan-genome toggle, plus a gene +
 *  ESM picker.  Other tools deep-link straight to a highlighted missense
 *  variant via the global goPanEffect().
 *
 *  Three panel behaviours make the single-item model work well:
 *    1. Breadcrumb reflects the active view   (SNPTools › PanEffect › B73)
 *    2. Smart default                          (canonical → Both, else B73;
 *                                               arriving from a variant → B73)
 *    3. Canonical-aware Pan option             (Pan/Both disabled off-canonical,
 *                                               with a link to the canonical id)
 *
 *  Rendering engine is pluggable:
 *    - 'embed'  : iframes the existing PanEffect site (default, works today).
 *    - 'native' : window.PanEffectEngine.render(...) once a port exists.
 * ===================================================================== */
(function () {
  'use strict';

  /* ------------------------------------------------------------------ *
   *  CONFIG                                                             *
   * ------------------------------------------------------------------ */
  const PE = {
    cfg: {
      /* 'native' renders PanEffect directly in the page via PanEffectEngine
         (paneffect-engine.js). Falls back to 'embed' automatically if that
         engine isn't loaded. Data is served same-origin from dataBase. */
      engine:     'native',
      /* local data folder (same-origin) used by BOTH the native engine and
         the panel's own canonical-transcript check */
      dataBase:   './paneffect/',
      /* only used by the 'embed' fallback */
      baseUrl:    'https://maizegdb.org/effect/maize_v2/index.html',
      esmModels:  ['ESM1', 'ESM2', 'ESM3'],
      defaultESM: 'ESM2',
      examples:   ['Zm00001eb260000', 'Zm00001eb268770_T001', 'lg1', 'wx1'],
    },
    /* current panel state; hydrated from S.pe on each render() */
    state: { gene: '', esm: 'ESM2', option: 'both', variant: null, wgs: null, resolved: null },
  };
  PE.cfg.dataBase = PE.cfg.dataBase ||
    PE.cfg.baseUrl.replace(/[^\/]*$/, '');   /* …/index.html -> …/ */

  const VIEW_LABEL = { both: 'Both views', b73: 'B73', pan: 'Pan-genome' };

  /* ------------------------------------------------------------------ *
   *  PUBLIC HANDOFF API                                                 *
   *  goPanEffect('Zm00001eb378140', {                                   *
   *      variant:{ pos:123, wt:'M', sub:'K' } | 'M123K',   // optional  *
   *      esm:'ESM2', option:'both'                          // optional  *
   *  })                                                                 *
   * ------------------------------------------------------------------ */
  window.goPanEffect = function (gene, opts) {
    opts = opts || {};
    const variant = normalizeVariant(opts.variant);
    S.pe = {
      gene:    gene || '',
      esm:     (opts.esm || PE.cfg.defaultESM),
      /* explicit option wins; else a variant jump lands on B73 (where the ESM
         substitution ring lives), a plain gene jump defaults to Both */
      option:  (opts.option || (variant ? 'b73' : 'both')),
      variant: variant,
      /* external modules open in the MaizeGDB 2026 view, not "all variants" */
      wgs:     (opts.wgs || 'maize2026'),
    };
    go('paneffect');
  };

  function variantLabel(v) {
    if (!v) return '';
    return `${v.wt || ''}${v.pos}${v.sub || ''}`;
  }

  /* variant may be a {pos,wt,sub} object or a raw SUB token (M123K /
     p.Met123Lys / Met123Lys). Returns null if no position can be recovered. */
  const AA3 = { Ala:'A',Arg:'R',Asn:'N',Asp:'D',Cys:'C',Gln:'Q',Glu:'E',Gly:'G',
    His:'H',Ile:'I',Leu:'L',Lys:'K',Met:'M',Phe:'F',Pro:'P',Ser:'S',Thr:'T',
    Trp:'W',Tyr:'Y',Val:'V',Ter:'*' };
  function normalizeVariant(v) {
    if (!v) return null;
    if (typeof v === 'object') return (v.pos != null) ? v : null;
    const s = String(v).trim().replace(/^p\./i, '');
    let m = s.match(/^([A-Za-z]{3})(\d+)([A-Za-z]{3})$/);
    if (m) return { wt: AA3[cap(m[1])] || '', pos: +m[2], sub: AA3[cap(m[3])] || '' };
    m = s.match(/^([A-Za-z])(\d+)([A-Za-z*])$/);
    if (m) return { wt: m[1].toUpperCase(), pos: +m[2], sub: m[3].toUpperCase() };
    return null;
  }
  function cap(w){ return w.charAt(0).toUpperCase() + w.slice(1).toLowerCase(); }

  /* ------------------------------------------------------------------ *
   *  RENDER — router entry point                                       *
   * ------------------------------------------------------------------ */
  function render(page) {
    const h = S.pe || {};
    PE.state.gene    = h.gene    != null ? h.gene    : PE.state.gene;
    PE.state.esm     = h.esm     || PE.state.esm || PE.cfg.defaultESM;
    PE.state.option  = h.option  || PE.state.option || 'both';
    PE.state.variant = ('variant' in h) ? h.variant : PE.state.variant;
    PE.state.wgs     = ('wgs' in h) ? h.wgs : PE.state.wgs;
    PE.state.resolved = null;
    S.pe = null; /* consume handoff */

    page.className = 'page fade';
    page.innerHTML = shellHTML();
    wireControls();
    updateCrumb();

    if (PE.state.gene && PE.state.gene.trim()) loadGene();
  }

  /* ------------------------------------------------------------------ *
   *  PANEL MARKUP                                                       *
   * ------------------------------------------------------------------ */
  function shellHTML() {
    const st = PE.state;
    const seg = (val, lbl) =>
      `<button class="pe-seg ${st.option === val ? 'on' : ''}" data-opt="${val}">${lbl}</button>`;
    const esmOpts = PE.cfg.esmModels
      .map(m => `<option value="${m}" ${m === st.esm ? 'selected' : ''}>${m}</option>`).join('');
    const examples = PE.cfg.examples
      .map(g => `<a href="#" class="pe-ex" data-gene="${g}">${g}</a>`)
      .join('<span class="pe-ex-sep">·</span>');
    const vLbl = variantLabel(st.variant);

    return `
    <style>${panelCSS()}</style>

    <div class="sec"><div class="bar"></div><div style="width:100%">
      <div class="n">MISSENSE VARIANT EFFECTS · ESM · B73 v5</div>
      <h2>Predicted effects of amino-acid substitutions</h2>
      <p>Pick a gene model and a protein language-model, then read every possible
         substitution as a heatmap — in B73 and across the pan-genome. Jump here
         from SNPVersity or SNPFold on a specific missense call and it lands
         pre-highlighted.</p>
    </div></div>

    <div class="card pe-controls">
      <div class="pe-row">
        <div class="field pe-gene">
          <label>Gene model, transcript, or protein</label>
          <input type="text" id="peGene" class="mono-in" placeholder="Zm00001eb…"
                 value="${escAttr(st.gene)}" autocomplete="off">
        </div>
        <div class="field pe-esm">
          <label>Variant score model</label>
          <select id="peEsm">${esmOpts}</select>
        </div>
        <div class="field pe-view">
          <label>Views</label>
          <div class="pe-segwrap" id="peSeg">
            ${seg('both', 'Both')}${seg('b73', 'B73 only')}${seg('pan', 'Pan-genome only')}
          </div>
        </div>
        <button class="btn pe-load" id="peLoad">Load views</button>
      </div>
      <div class="pe-ex-row">Examples: ${examples}</div>
      ${vLbl ? `<div class="pe-vchip" id="peVchip">Highlighting <b>${vLbl}</b>
        <button class="pe-vclear" id="peVclear" title="Clear highlight">×</button></div>` : ''}
      <div class="pe-canon-hint" id="peCanonHint" hidden></div>
      <div class="status" id="peStatus"></div>
    </div>

    <div id="peView" class="pe-view-host"></div>
    `;
  }

  /* ------------------------------------------------------------------ *
   *  CONTROL WIRING                                                     *
   * ------------------------------------------------------------------ */
  function wireControls() {
    const gene = document.getElementById('peGene');
    const esm  = document.getElementById('peEsm');
    const load = document.getElementById('peLoad');

    gene.addEventListener('keydown', e => { if (e.key === 'Enter') { e.preventDefault(); readAndLoad(); } });
    esm.addEventListener('change', () => {
      PE.state.esm = esm.value;
      if (PE.state.gene) loadViews();  /* ESM change doesn't affect canonical status */
    });
    load.addEventListener('click', readAndLoad);

    document.querySelectorAll('.pe-ex').forEach(a => a.addEventListener('click', e => {
      e.preventDefault();
      gene.value = a.dataset.gene;
      PE.state.variant = null; /* a fresh gene, not a variant jump */
      readAndLoad();
    }));

    document.querySelectorAll('.pe-seg').forEach(b => b.addEventListener('click', () => {
      if (b.classList.contains('disabled')) { flashCanonHint(); return; }
      PE.state.option = b.dataset.opt;
      setSegStates(currentDisabled());
      updateCrumb();
      if (PE.state.gene) loadViews();
    }));

    const vclear = document.getElementById('peVclear');
    if (vclear) vclear.addEventListener('click', () => {
      PE.state.variant = null;
      const chip = document.getElementById('peVchip');
      if (chip) chip.remove();
      if (PE.state.gene) loadViews();
    });
  }

  function readAndLoad() {
    const gene = document.getElementById('peGene').value.trim();
    PE.state.gene = gene;
    PE.state.esm  = document.getElementById('peEsm').value;
    PE.state.wgs  = null;   /* a manual in-panel search starts at "all variants" */
    if (!gene) { status('Enter a gene model, transcript, or protein id to load.'); return; }
    loadGene();
  }

  function status(msg) {
    const el = document.getElementById('peStatus');
    if (el) el.textContent = msg || '';
  }

  /* ------------------------------------------------------------------ *
   *  LOAD FLOW — resolve canonical, apply policy, then render          *
   * ------------------------------------------------------------------ */
  async function loadGene() {
    status('Checking transcript…');
    let info = null;
    try { info = await resolveGene(PE.state.gene); } catch (e) { info = null; }
    PE.state.resolved = info;
    applyCanonicalPolicy(info);   /* may coerce option + disable Pan */
    updateCrumb();
    loadViews();
  }

  /* Smart default + canonical-aware Pan. If we couldn't resolve (file missing,
     offline), we leave everything enabled so the user is never blocked. */
  function applyCanonicalPolicy(info) {
    const hint = document.getElementById('peCanonHint');

    if (info && info.canonical === false) {
      /* Pan and Both both surface the pan view, which needs the canonical
         transcript — so restrict to B73 and point the user at the canonical id */
      if (PE.state.option !== 'b73') PE.state.option = 'b73';
      setSegStates(new Set(['both', 'pan']));
      if (hint) {
        hint.hidden = false;
        hint.innerHTML =
          `Pan-genome view is only available for the canonical transcript. ` +
          `Showing <b>B73</b> for ${escAttr(info.rawId || PE.state.gene)}. ` +
          (info.canonicalTranscript
            ? `<a href="#" id="peCanonLink">Load canonical ${escAttr(info.canonicalTranscript)} →</a>`
            : '');
        const link = document.getElementById('peCanonLink');
        if (link) link.addEventListener('click', e => { e.preventDefault(); loadCanonical(); });
      }
    } else {
      setSegStates(new Set());               /* all enabled */
      if (hint) { hint.hidden = true; hint.innerHTML = ''; }
    }
  }

  function currentDisabled() {
    const r = PE.state.resolved;
    return (r && r.canonical === false) ? new Set(['both', 'pan']) : new Set();
  }

  function setSegStates(disabledSet) {
    document.querySelectorAll('.pe-seg').forEach(b => {
      const opt = b.dataset.opt;
      b.classList.toggle('disabled', disabledSet.has(opt));
      b.classList.toggle('on', opt === PE.state.option);
      if (disabledSet.has(opt)) b.title = 'Only available for the canonical transcript';
      else b.removeAttribute('title');
    });
  }

  function flashCanonHint() {
    const hint = document.getElementById('peCanonHint');
    if (!hint || hint.hidden) return;
    hint.classList.remove('flash'); void hint.offsetWidth; hint.classList.add('flash');
  }

  function loadCanonical() {
    const t = PE.state.resolved && PE.state.resolved.canonicalTranscript;
    if (!t) return;
    const input = document.getElementById('peGene');
    if (input) input.value = t;
    PE.state.gene = t;
    PE.state.variant = null;                 /* canonical, not the original variant call */
    const chip = document.getElementById('peVchip'); if (chip) chip.remove();
    loadGene();
  }

  function updateCrumb() {
    const el = document.getElementById('crumbTool');
    if (!el) return;
    const label = VIEW_LABEL[PE.state.option] || '';
    const caret = (typeof ICONS !== 'undefined' && ICONS.caret) ? ICONS.caret : '›';
    el.innerHTML = '<b>PanEffect</b>' +
      (label ? ` <span class="pe-crumb-sep">${caret}</span> <span class="pe-crumb-view">${label}</span>` : '');
  }

  /* ------------------------------------------------------------------ *
   *  CANONICAL RESOLUTION — small fetch against PanEffect's data dir    *
   *  (same host). Best-effort: any failure returns canonical:true.      *
   * ------------------------------------------------------------------ */
  async function resolveGene(rawId) {
    const out = { rawId: rawId, geneModel: '', transcript: '', protein: '',
                  canonicalTranscript: null, canonical: true };
    let id = String(rawId || '').trim();
    if (!id) return out;

    let canFlag = false;
    if (id.includes('_T'))      { out.geneModel = id.split('_')[0]; out.transcript = id;                       out.protein = id.replace('_T', '_P'); }
    else if (id.includes('_P')) { out.geneModel = id.split('_')[0]; out.transcript = id.replace('_P', '_T');   out.protein = id; }
    else                        { out.geneModel = id; canFlag = true; out.transcript = id + '_T001';           out.protein = id + '_P001'; }

    /* names like lg1 / wx1 → gene model via the synonym table */
    if (!/^zm\d/i.test(out.geneModel)) {
      const syn = await lookupSynonym(out.geneModel);
      if (syn) {
        out.geneModel = syn;
        if (canFlag) { out.transcript = syn + '_T001'; out.protein = syn + '_P001'; }
      }
    }

    const gmCan = await fetchCanonicalTranscript(out.geneModel);
    if (gmCan) {
      out.canonicalTranscript = gmCan;
      if (canFlag) {
        /* a bare gene adopts the canonical transcript → always canonical */
        out.transcript = gmCan; out.protein = gmCan.replace('_T', '_P'); out.canonical = true;
      } else {
        out.canonical = (out.protein === gmCan.replace('_T', '_P'));
      }
    }
    return out;
  }

  async function fetchCanonicalTranscript(geneModel) {
    try {
      const r = await fetch(PE.cfg.dataBase + 'uniprot/' + geneModel + '.tsv');
      if (!r.ok) return null;
      const txt = await r.text();
      if (/<html|<script/i.test(txt)) return null;
      let gmCan = null;
      txt.split('\n').forEach(row => {
        const c = row.split('\t');
        if (c.length >= 7 && c[1]) gmCan = c[1];   /* col 1 = canonical transcript */
      });
      return gmCan;
    } catch (e) { return null; }
  }

  async function lookupSynonym(name) {
    try {
      const r = await fetch(PE.cfg.dataBase + 'synonym/maize_synonym.tsv');
      if (!r.ok) return null;
      const txt = await r.text();
      if (/<html|<script/i.test(txt)) return null;
      const lines = txt.split('\n');
      for (const line of lines) {
        const c = line.split(/\s+|\t+/);
        if (c.length >= 2 && c[1] === name) return c[0];
      }
      return null;
    } catch (e) { return null; }
  }

  /* ------------------------------------------------------------------ *
   *  ENGINE DISPATCH                                                    *
   * ------------------------------------------------------------------ */
  function loadViews() {
    const host = document.getElementById('peView');
    if (!host) return;

    const nativeReady =
      window.PanEffectEngine && typeof window.PanEffectEngine.render === 'function' &&
      typeof window.runPanEffect === 'function';   /* ported pipeline present */

    if (PE.cfg.engine === 'native' && nativeReady) {
      renderNative(host);
    } else {
      if (PE.cfg.engine === 'native') {
        status(window.PanEffectEngine
          ? 'PanEffect scripts (./js/pe/*.js) aren’t loaded — showing embedded view. See console.'
          : 'Native engine not loaded — showing embedded view.');
      }
      renderEmbed(host);
    }
  }

  function renderEmbed(host) {
    const st = PE.state;
    const p = new URLSearchParams();
    p.set('id', st.gene);
    p.set('option', st.option);
    p.set('esm', st.esm);
    p.set('embed', '1');
    if (st.variant) {
      p.set('highlight', variantLabel(st.variant));
      if (st.variant.pos != null) p.set('pos', String(st.variant.pos));
      if (st.variant.sub)         p.set('sub', st.variant.sub);
    }
    const url = PE.cfg.baseUrl + '?' + p.toString();

    host.innerHTML =
      `<div class="pe-frame-wrap">
         <iframe id="peFrame" class="pe-frame" title="PanEffect — ${escAttr(st.gene)}"
                 src="${escAttr(url)}" loading="lazy"></iframe>
       </div>
       <div class="pe-embed-note">
         Rendered by the PanEffect site in an isolated frame.
         <a href="${escAttr(url)}" target="_blank" rel="noopener">Open full page ↗</a>
       </div>`;
    status('Loading ' + st.gene + ' …');
    const frame = document.getElementById('peFrame');
    if (frame) frame.addEventListener('load', () => status(''));
  }

  function renderNative(host) {
    host.innerHTML = '<div id="peNativeRoot" class="pe-native-root"></div>';
    const root = document.getElementById('peNativeRoot');
    onLoadingIcon_safe();
    Promise.resolve(
      window.PanEffectEngine.render(root, {
        gene:     PE.state.gene,
        transcript: PE.state.resolved && PE.state.resolved.transcript,
        protein:  PE.state.resolved && PE.state.resolved.protein,
        esm:      PE.state.esm,
        option:   PE.state.option,
        variant:  PE.state.variant,
        wgs:      PE.state.wgs,
        dataBase: PE.cfg.dataBase,
      })
    ).then(() => { status(''); offLoadingIcon_safe(); })
     .catch(err => {
        console.error('[paneffect] native render error:', err);
        status('That id could not be loaded. Check the gene model and try again.');
        offLoadingIcon_safe();
     });
  }

  function onLoadingIcon_safe()  { try { if (typeof onLoadingIcon  === 'function') onLoadingIcon();  } catch (e) {} }
  function offLoadingIcon_safe() { try { if (typeof offLoadingIcon === 'function') offLoadingIcon(); } catch (e) {} }

  /* ------------------------------------------------------------------ *
   *  utils                                                              *
   * ------------------------------------------------------------------ */
  function escAttr(s) {
    return String(s == null ? '' : s)
      .replace(/&/g, '&amp;').replace(/"/g, '&quot;')
      .replace(/</g, '&lt;').replace(/>/g, '&gt;');
  }

  function panelCSS() {
    return `
    .pe-controls{padding:16px;margin-top:6px}
    .pe-row{display:flex;gap:14px;align-items:flex-end;flex-wrap:wrap}
    .pe-gene{flex:2 1 320px;margin:0}
    .pe-esm{flex:0 0 160px;margin:0}
    .pe-view{flex:1 1 260px;margin:0}
    .pe-load{align-self:flex-end;white-space:nowrap}
    .pe-segwrap{display:inline-flex;border:1px solid var(--line,#e6e9ef);border-radius:9px;overflow:hidden}
    .pe-seg{appearance:none;border:0;background:#fff;padding:9px 12px;font:inherit;font-size:13px;
      cursor:pointer;color:var(--muted,#5b6472);border-right:1px solid var(--line,#e6e9ef)}
    .pe-seg:last-child{border-right:0}
    .pe-seg.on{background:var(--blue-600,#2563eb);color:#fff}
    .pe-seg.disabled{color:var(--faint,#9aa3b2);background:#f6f7f9;cursor:not-allowed}
    .pe-seg.disabled.on{background:#f6f7f9;color:var(--faint,#9aa3b2)}
    .pe-ex-row{margin-top:12px;font-size:12.5px;color:var(--muted,#5b6472)}
    .pe-ex{color:var(--blue-600,#2563eb);text-decoration:none;margin:0 6px}
    .pe-ex:hover{text-decoration:underline}
    .pe-ex-sep{color:var(--faint,#9aa3b2)}
    .pe-vchip{display:inline-flex;align-items:center;gap:8px;margin-top:12px;padding:6px 10px;
      background:#eef4ff;border:1px solid #d5e2ff;border-radius:999px;font-size:12.5px;color:#1b3b7a}
    .pe-vclear{border:0;background:transparent;cursor:pointer;font-size:15px;line-height:1;color:#1b3b7a}
    .pe-canon-hint{margin-top:12px;padding:9px 12px;border-radius:9px;font-size:12.5px;
      background:#fff8ec;border:1px solid #f3e3c2;color:#7a5a17}
    .pe-canon-hint a{color:var(--blue-600,#2563eb);margin-left:4px}
    .pe-canon-hint.flash{animation:pe-flash .6s ease}
    @keyframes pe-flash{0%{background:#ffe9bf}100%{background:#fff8ec}}
    .pe-crumb-sep{display:inline-flex;vertical-align:middle;opacity:.5}
    .pe-crumb-view{color:var(--muted,#5b6472)}
    .pe-view-host{margin-top:18px}
    .pe-frame-wrap{border:1px solid var(--line,#e6e9ef);border-radius:12px;overflow:hidden;background:#fff}
    .pe-frame{width:100%;height:1400px;border:0;display:block}
    .pe-embed-note{margin-top:8px;font-size:12px;color:var(--faint,#9aa3b2)}
    .pe-embed-note a{color:var(--blue-600,#2563eb)}
    @media (max-width:640px){ .pe-frame{height:1100px} .pe-load{width:100%} }
    `;
  }

  /* ------------------------------------------------------------------ *
   *  register with the suite shell                                     *
   * ------------------------------------------------------------------ */
  SNPTools.register('paneffect', { render: render });
})();
