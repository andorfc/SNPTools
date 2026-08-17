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
      engine:     'native',
      /* Fusarium PanEffect data is per-genome under ./paneffect/<refDir>/.
         The reference is chosen from the gene-id prefix (FGSG_/FVEG_/FVERT4_),
         the same convention SNPFold uses for structures/domains. `dataBase`
         (state) is set per load; `dataRoot` is the shared parent. */
      dataRoot:   './paneffect/',
      refDirs:    { graminearum:'Fgram_ph1', vert7600:'Fvert_7600', vertMRC826:'Fvert_mrc' },
      refLabels:  { graminearum:'F. graminearum PH-1', vert7600:'F. verticillioides 7600', vertMRC826:'F. verticillioides MRC826' },
      /* four protein language models; ESM4 = ESM-C (shown as "ESM-C") */
      esmModels:  ['ESM1', 'ESM2', 'ESM3', 'ESM4'],
      esmLabel:   { ESM4:'ESM-C' },
      defaultESM: 'ESM2',
      examples:   ['FGSG_00134', 'FVEG_000452', 'FVERT4_000426'],
    },
    /* current panel state; hydrated from S.pe on each render() */
    state: { gene: '', esm: 'ESM2', option: 'both', variant: null, wgs: null, resolved: null,
             dataBase: './paneffect/Fgram_ph1/', refLabel: '', coverage: 'hq' },
  };

  const COVERAGE_LABELS = { hq: 'High Quality', hc: 'High Coverage' };
  /* observed-variant cache keyed "gene|coverage" — the missense set for a gene
     depends only on the gene + HQ/HC dataset, not on ESM or the active view */
  const observedCache = {};

  const VIEW_LABEL = { both: 'Both views', b73: 'Reference', pan: 'Across species' };

  /* ---- reference (species) resolution from the gene-id prefix ---- */
  function familyFromGene(gene) {
    const g = String(gene || '').toUpperCase();
    if (g.indexOf('FGSG')   === 0) return 'graminearum';
    if (g.indexOf('FVERT4') === 0) return 'vertMRC826';   // test before FVEG
    if (g.indexOf('FVEG')   === 0) return 'vert7600';
    try { if (typeof Data !== 'undefined' && Data.familyOf && typeof S !== 'undefined' && S) return Data.familyOf(S.dataset); } catch (e) {}
    return 'graminearum';
  }
  /* strip a transcript/protein suffix (_T001 / _P001 / -T1) to the gene model */
  function geneModelOf(id) { return String(id == null ? '' : id).trim().replace(/[_-][TP]\d+$/, ''); }
  /* canonical (padding-correct) gene id, via the app's shared normalizer */
  function canonGene(id) {
    try { if (typeof Data !== 'undefined' && Data.canonicalGeneId) return Data.canonicalGeneId(id); } catch (e) {}
    return id;
  }
  /* point the engine's data + synonym + label at the gene's genome */
  function applyRef(gene) {
    const fam = familyFromGene(gene);
    const dir = PE.cfg.refDirs[fam] || 'Fgram_ph1';
    PE.state.dataBase = PE.cfg.dataRoot + dir + '/';
    PE.state.refLabel = PE.cfg.refLabels[fam] || '';
    window.__PE_SYNONYM = './synonym/' + dir + '_synonym.tsv';
  }

  /* ---- observed-variant filter: HDF5 missense set for the reference view ---- */
  /* pick the HQ/HC dataset id for the gene's family (graminearum / vert7600 / vertMRC826) */
  function datasetFor(fam, coverage) {
    try {
      const want = coverage === 'hc' ? 'High Coverage' : 'High Quality';
      const ds = Data.datasets().find(d => d.family === fam && d.sub === want);
      return ds || null;
    } catch (e) { return null; }
  }
  /* Query the HDF5 (same path SNPVersity/SNPImpact use) for the gene's variants
     and keep the missense substitutions as a Set of "resi|altAA" keys. Cached
     per gene+coverage. Best-effort: any failure returns an empty set + error. */
  async function getObserved(gene, coverage) {
    const key = gene + '|' + coverage;
    if (observedCache[key]) return observedCache[key];
    const fam = familyFromGene(gene);
    const ds = datasetFor(fam, coverage);
    const datasetName = ds ? (ds.name + ' · ' + ds.sub) : '';
    const set = new Set();
    let error = null;
    if (!ds) {
      error = 'no-dataset';
    } else {
      try {
        /* pass every accession so "observed" spans the whole dataset, not just
           the default founder subset (variant rows are keyed on the genotypes sent) */
        const allIds = Data.accessionsFor(ds.id).map(a => a.id);
        const vars = await Data.queryFoldVariants(gene, ds.id, allIds);
        (vars || []).forEach(v => {
          if (v && v.consClass === 'missense' && v.resi != null &&
              v.alt && v.alt !== '*' && !/fs/i.test(String(v.alt))) {
            set.add(v.resi + '|' + v.alt);
          }
        });
      } catch (e) {
        error = (e && e.message) || 'query-failed';
        console.warn('[paneffect] observed-variant query failed:', error);
      }
    }
    const res = { set, count: set.size, datasetName, coverage, error };
    observedCache[key] = res;
    return res;
  }
  /* Build the config object the engine consumes (initial set + count + the
     HQ/HC switch callback). Default to the observed view when the gene has
     observed missense calls, else fall back to "all" so the view isn't blank. */
  async function buildObservedConfig(gene) {
    const res = await getObserved(gene, PE.state.coverage);
    return {
      set: res.set, count: res.count,
      coverage: PE.state.coverage, coverageLabels: COVERAGE_LABELS,
      datasetName: res.datasetName, error: res.error,
      mode: (res.count > 0 ? 'observed' : 'all'),
      onCoverage: async (cov) => {
        PE.state.coverage = cov;
        const r = await getObserved(gene, cov);
        return { set: r.set, count: r.count, datasetName: r.datasetName, error: r.error };
      },
    };
  }

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
      /* Fusarium stores have a single "all variants" track */
      wgs:     (opts.wgs || null),
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
      .map(m => `<option value="${m}" ${m === st.esm ? 'selected' : ''}>${PE.cfg.esmLabel[m] || m}</option>`).join('');
    const examples = PE.cfg.examples
      .map(g => `<a href="#" class="pe-ex" data-gene="${g}">${g}</a>`)
      .join('<span class="pe-ex-sep">·</span>');
    const vLbl = variantLabel(st.variant);

    return `
    <style>${panelCSS()}</style>

    <div class="sec"><div class="bar"></div><div style="width:100%">
      <div class="n">MISSENSE VARIANT EFFECTS · ESM PROTEIN LANGUAGE MODELS · FUSARIUM</div>
      <h2>Predicted effects of amino-acid substitutions</h2>
      <p>Pick a gene model and a protein language-model, then read every possible
         substitution as a heatmap — in the reference genome and across related
         Fusarium species. Jump here from SNPVersity or SNPFold on a specific
         missense call and it lands pre-highlighted.</p>
    </div></div>

    <div class="card pe-controls">
      <div class="pe-row">
        <div class="field pe-gene">
          <label>Gene model, transcript, or protein</label>
          <input type="text" id="peGene" class="mono-in" placeholder="FGSG_… / FVEG_… / FVERT4_…"
                 value="${escAttr(st.gene)}" autocomplete="off">
        </div>
        <div class="field pe-esm">
          <label>Variant score model</label>
          <select id="peEsm">${esmOpts}</select>
        </div>
        <div class="field pe-view">
          <label>Views</label>
          <div class="pe-segwrap" id="peSeg">
            ${seg('both', 'Both')}${seg('b73', 'Reference only')}${seg('pan', 'Across species')}
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
    /* normalize to the canonical bare gene model, then point at its genome */
    PE.state.gene = canonGene(geneModelOf(PE.state.gene));
    const inp = document.getElementById('peGene'); if (inp && PE.state.gene) inp.value = PE.state.gene;
    applyRef(PE.state.gene);
    let info = null;
    try { info = await resolveGene(PE.state.gene); } catch (e) { info = null; }
    PE.state.resolved = info;
    applyCanonicalPolicy(info);   /* single-isoform → all views enabled */
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
  /* Fusarium genomes are single-isoform: the gene model has exactly one
     transcript (_T001) and one protein (_P001), which is always canonical, so
     the pan-genome (across-species) view is always available. No transcript
     disambiguation or common-name synonym lookup is needed here. */
  async function resolveGene(rawId) {
    const gm = canonGene(geneModelOf(rawId));
    return { rawId: rawId, geneModel: gm,
             transcript: gm ? gm + '_T001' : '', protein: gm ? gm + '_P001' : '',
             canonicalTranscript: gm ? gm + '_T001' : null, canonical: true };
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
  async function loadViews() {
    const host = document.getElementById('peView');
    if (!host) return;

    const nativeReady =
      window.PanEffectEngine && typeof window.PanEffectEngine.render === 'function' &&
      typeof window.runPanEffect === 'function';   /* ported pipeline present */

    if (PE.cfg.engine === 'native' && nativeReady) {
      /* fetch the observed-missense set for the reference-view filter before
         rendering, so the initial paint already reflects the chosen mode */
      let observed = null;
      try { observed = await buildObservedConfig(PE.state.gene); }
      catch (e) { console.warn('[paneffect] observed config failed:', e); }
      renderNative(host, observed);
    } else {
      if (PE.cfg.engine === 'native') {
        status(window.PanEffectEngine
          ? 'PanEffect scripts (./js/pe/*.js) aren’t loaded — showing embedded view. See console.'
          : 'Native engine not loaded — showing embedded view.');
      }
      renderEmbed(host);
    }
  }

  /* Fusarium has no external PanEffect site to embed — the native engine is the
     only renderer. This only shows if the pe/*.js scripts failed to load. */
  function renderEmbed(host) {
    host.innerHTML =
      `<div class="pe-canon-hint" style="display:block">
         The PanEffect engine scripts didn’t load, so the heatmap can’t be drawn.
         Check that <span class="mono">js/pe/*.js</span>, <span class="mono">paneffect-heatmap.js</span>,
         <span class="mono">paneffect-engine.js</span> and d3 are included in <span class="mono">index.html</span>
         (see the browser console / Network tab for 404s).
       </div>`;
    status('PanEffect engine not loaded.');
  }

  function renderNative(host, observed) {
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
        dataBase: PE.state.dataBase,
        refLabel: PE.state.refLabel,
        observed: observed || null,
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
