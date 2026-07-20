/* =====================================================================
 *  snpfold.js — Structure-aware variant interpretation.
 *
 *  Integrates coding variants with predicted protein structure three ways:
 *    1. a linear "protein browser" track that aligns variants, Pfam domains,
 *       secondary structure, and per-residue pLDDT confidence by residue;
 *    2. an interactive 3D viewer (3Dmol.js, loaded on demand) colorable by
 *       confidence / domain / variant impact, with variant residues marked;
 *    3. a per-variant structural-context readout that interprets each change
 *       (domain, local confidence, secondary structure, predicted ΔΔG).
 *
 *  Structure data + PDB come from structure-<gene>.js via Data.structureFor()
 *  / Data.pdbFor(); variants from Data.queryFoldVariants().
 * ===================================================================== */
(function () {
  const THREEDMOL_URL = 'https://3Dmol.org/build/3Dmol-min.js';
  const W = 1000;   // svg user-units width (scales to container)

  const FD = {
    gene: 'Zm00001eb406050',
    struct: null, pdb: null, variants: [],
    selId: null,
    colorMode: 'plddt',     // plddt | domain | impact
    showVar: true,
    viewer: null, libState: 'idle',
    dataset: null, sec: false,  // sec = show PlantCAD2/ESM2/ESM3 (MaizeGDB 2026 only)
    carriers: null, openCarrier: null,   // pos|ref|alt -> {carriersHom,carriersHet,het,hom} (whole-panel, via geneFunction)
    sort: { key: null, dir: 'asc' },     // variant-table sort: column key + direction ('asc'|'desc'); null key = file order
    root: null,          // persistent DOM container — survives navigation to other tools
    loaded: false,       // a gene's heavy content is (being) rendered into root
    loadedGene: null,    // which gene that content is for
  };

  /* ---------- palettes ---------- */
  function plddtHex(b){ return b>=90?'#0053d6': b>=70?'#65cbf3': b>=50?'#ffdb13':'#ff7d45'; }
  function plddtBand(b){ return b>=90?'Very high': b>=70?'Confident': b>=50?'Low':'Very low'; }
  function impactHex(s){
    const n = finiteNumber(s);
    if (n == null) return '#8a93a3';
    const lo=-12,hi=4,t=Math.max(0,Math.min(1,(n-lo)/(hi-lo)));
    const r=t<.5?220:Math.round(220*(1-(t-.5)*2)), g=t<.5?Math.round(180*t*2):170;
    return `rgb(${r},${Math.max(50,g)},55)`;
  }
  const DOM_FILL = ['#cdeccf','#cfe0fb','#efe0fb','#fde9cc'];
  const CONS_FILL = { lof:'#d6322a', missense:'#2f5bbf', lod:'#b54708', splice:'#6d28d9', indel:'#176c3a', syn:'#8a93a3' };
  const SS_LABEL = { H:'α-helix', E:'β-strand', C:'loop / coil' };

  /* language-model score pill — same gradient (red→green) as SNPImpact / SNPFunction */
  function finiteNumber(v){
    if (v == null || v === '' || v === '.' || v === 'NA' || v === 'N/A') return null;
    const n = Number(v);
    return Number.isFinite(n) ? n : null;
  }
  function scoreColor(v){
    const n = finiteNumber(v);
    if (n == null) return '#8a93a3';
    const lo=-12,hi=6,t=Math.max(0,Math.min(1,(n-lo)/(hi-lo)));
    const r=t<.5?255:Math.round(255*(1-(t-.5)*2));
    const g=t<.5?Math.round(255*t*2):200;
    return `rgb(${r},${Math.max(60,g)},60)`;
  }
  function scoreCell(v){
    const n = finiteNumber(v);
    return n==null
      ? '<span style="color:var(--faint)">—</span>'
      : `<span class="imp-score" style="background:${scoreColor(n)}">${n>0?'+':''}${n.toFixed(1)}</span>`;
  }
  function scoreText(v){
    const n = finiteNumber(v);
    return n == null ? '—' : `${n > 0 ? '+' : ''}${n.toFixed(1)}`;
  }

  /* Support both the established fields and common 2026 aliases. */
  const SCORE_KEYS = {
    plantcad:  ['plantcad', 'plantcad1', 'plantcad_1', 'PLANTCAD', 'PLANTCAD1'],
    plantcad2: ['plantcad2', 'plantcad_2', 'PLANTCAD2'],
    esm:       ['esm', 'esm1', 'esm_1', 'ESM', 'ESM1'],
    esm2:      ['esm2', 'esm_2', 'ESM2'],
    esm3:      ['esm3', 'esm_3', 'ESM3'],
  };
  function modelScore(v, model){
    if (!v) return null;
    const keys = SCORE_KEYS[model] || [model];
    for (const key of keys){
      if (Object.prototype.hasOwnProperty.call(v, key) && finiteNumber(v[key]) != null) return finiteNumber(v[key]);
    }
    return null;
  }
  function hasSecondaryVariantScores(variants){
    return (variants || []).some(v =>
      modelScore(v, 'plantcad2') != null ||
      modelScore(v, 'esm2') != null ||
      modelScore(v, 'esm3') != null
    );
  }
  function aiScoreSummary(v){
    const scores = [
      ['PlantCAD', modelScore(v, 'plantcad')],
      ...(FD.sec ? [['PlantCAD2', modelScore(v, 'plantcad2')]] : []),
      ['ESM1', modelScore(v, 'esm')],
      ...(FD.sec ? [
        ['ESM2', modelScore(v, 'esm2')],
        ['ESM3', modelScore(v, 'esm3')],
      ] : []),
    ];
    if (!scores.some(([, value]) => value != null)) return 'n/a';
    //return scores.map(([label, value]) => `${label} ${scoreText(value)}`).join(' · ');
    return scores.map(([label, value]) => `${label} ${scoreText(value)}`).join('<br>');
  }

  function escFold(s){ return String(s==null?'':s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }
  function carrierKey(v){ return v.pos+'|'+v.refNt+'|'+v.altNt; }
  function carrierOf(v){ return FD.carriers ? (FD.carriers[carrierKey(v)] || null) : null; }

  /* ---------- active dataset detection ---------- */
  function nonEmpty(v){ return v != null && v !== ''; }
  function datasetText(dataset){
    if (dataset == null) return '';
    if (typeof dataset === 'string' || typeof dataset === 'number') return String(dataset);
    if (typeof dataset === 'object'){
      return [dataset.id, dataset.key, dataset.value, dataset.name, dataset.label, dataset.title]
        .filter(nonEmpty).join(' ');
    }
    return String(dataset);
  }
  function datasetRepresentations(dataset){
    const out = [];
    const add = v => {
      if (!nonEmpty(v)) return;
      if (!out.some(x => x === v)) out.push(v);
    };
    add(dataset);
    if (dataset && typeof dataset === 'object'){
      add(dataset.id); add(dataset.key); add(dataset.value);
      add(dataset.name); add(dataset.label); add(dataset.title);
    }
    return out;
  }
  function resolveDataset(){
    const candidates = [];
    const add = v => { if (nonEmpty(v)) candidates.push(v); };

    if (typeof S !== 'undefined' && S){
      add(S.dataset);
      add(S.datasetId);
      add(S.datasetKey);
      add(S.selectedDataset);
      add(S.activeDataset);
      add(S.currentDataset);
      add(S.snpDataset);
      add(S.ds);
    }

    if (typeof Data !== 'undefined' && Data){
      if (typeof Data.currentDataset === 'function'){
        try { add(Data.currentDataset()); } catch (e) { /* use other candidates */ }
      } else {
        add(Data.currentDataset);
      }
      add(Data.dataset);
      add(Data.datasetId);
      add(Data.datasetKey);
      add(Data.selectedDataset);
      add(Data.activeDataset);
    }

    if (typeof window !== 'undefined'){
      add(window.CURRENT_DATASET);
      add(window.currentDataset);
      add(window.activeDataset);
    }

    add(FD.dataset);
    return candidates.length ? candidates[0] : null;
  }
  function hasSecondaryScores(dataset){
    /* Prefer the application's own helper, while trying both object and ID/name forms. */
    if (typeof Data !== 'undefined' && typeof Data.hasSecondaryScores === 'function'){
      const reps = datasetRepresentations(dataset);
      if (!reps.length) reps.push(undefined);
      for (const rep of reps){
        try {
          if (Data.hasSecondaryScores(rep) === true) return true;
        } catch (e) { /* fall through to normalized-name check */ }
      }
    }

    /* Fallback accepts labels/IDs such as MaizeGDB2026, maizegdb_2026_hq,
       or "MaizeGDB 2026 (High Coverage)". */
    const normalized = datasetText(dataset).toLowerCase().replace(/[^a-z0-9]+/g, '');
    return normalized.includes('maizegdb2026');
  }

  /* ---------- dataset chooser (compact; shares Data.datasets() with SNPVersity) ---------- */
  function foldDatasets(){
    try {
      if (typeof Data !== 'undefined' && typeof Data.datasets === 'function'){
        const ds = Data.datasets();
        if (Array.isArray(ds)) return ds;
      }
    } catch (e) { /* no dataset catalog available */ }
    return [];
  }
  function datasetId(d){
    if (d && typeof d === 'object'){
      if (nonEmpty(d.id))  return d.id;
      if (nonEmpty(d.key)) return d.key;
      return datasetText(d);
    }
    return d;
  }
  function isCurrentDataset(d){
    const a = String(datasetId(d)), b = String(datasetId(FD.dataset));
    if (a && a === b) return true;
    const na = datasetText(d).toLowerCase().replace(/[^a-z0-9]+/g, '');
    const nb = datasetText(FD.dataset).toLowerCase().replace(/[^a-z0-9]+/g, '');
    return !!na && na === nb;
  }
  function datasetLabel(d){
    const name = (d && nonEmpty(d.name)) ? String(d.name) : datasetText(d);
    let sub = (d && nonEmpty(d.sub)) ? String(d.sub) : '';
    if (sub && name.toLowerCase().includes(sub.toLowerCase())) sub = '';
    return { name: name || 'Dataset', sub };
  }
  function datasetCardsHTML(){
    const list = foldDatasets();
    if (!list.length) return '';
    return list.map(d=>{
      const id  = datasetId(d);
      const sel = isCurrentDataset(d);
      const { name, sub } = datasetLabel(d);
      return `<button type="button" class="fold-ds ${sel?'sel':''}" onclick="FOLD.pickDataset('${escFold(id)}')">
        <span class="fold-ds-dot"></span>
        <span class="fold-ds-txt"><span class="fold-ds-name">${escFold(name)}</span>${sub?`<span class="fold-ds-sub">${escFold(sub)}</span>`:''}</span>
      </button>`;
    }).join('');
  }
  function datasetChooser(){
    const cards = datasetCardsHTML();
    if (!cards) return '';   // if the data layer exposes no catalog, show nothing
    return `<div class="card pad fold-ds-card" style="margin-bottom:16px">
      <div class="fold-ds-head">Dataset</div>
      <div class="fold-ds-grid" id="foldDsGrid">${cards}</div>
    </div>`;
  }
  function syncDatasetChooser(){
    const grid = FD.root && FD.root.querySelector('#foldDsGrid');
    if (grid) grid.innerHTML = datasetCardsHTML();
  }

  /* ---------- structural context for a residue ---------- */
  function domainAt(resi){
    const ds = FD.struct.domains || [];
    return ds.find(d=>d.kind==='domain' && resi>=d.start && resi<=d.end)
        || ds.find(d=>resi>=d.start && resi<=d.end) || null;
  }
  function ctxFor(v){
    const i = v.resi-1;
    const inModel = FD.struct.plddt && i>=0 && i<FD.struct.plddt.length;
    const plddt = inModel ? FD.struct.plddt[i] : null;
    const ss = inModel && FD.struct.ss && FD.struct.ss[i] ? FD.struct.ss[i] : 'C';
    return { plddt, ss, ssLabel:SS_LABEL[ss], domain: inModel ? domainAt(v.resi) : null, inModel };
  }
  function interpret(v, c){
    if (!c.inModel){
      return ['mid', `Residue ${v.resi} is outside the modeled region (1–${FD.struct.length} aa) — shown in the list, but not placed on the structure.`];
    }
    const inDom = c.domain && c.domain.kind==='domain';
    const where = inDom ? (' in the '+c.domain.name+' domain')
                        : (c.domain && c.domain.kind==='region' ? (' in the '+c.domain.name) : ' in an inter-domain region');
    if (v.consClass==='lof'){
      const after = inDom ? ('within / before the '+c.domain.name+' domain') : 'in the C-terminal region';
      const what = v.consequence==='Start lost' ? 'Disrupts the start codon'
                 : v.consequence==='Stop lost'  ? 'Removes the stop codon'
                 : 'Truncates the protein';
      return ['lof', `${what} ${after} — position shown; downstream effect not modeled.`];
    }
    if (v.consClass==='indel'){
      return ['mid', `${v.consequence}${where} — in-frame length change; position shown, folding effect not modeled.`];
    }
    if (v.combined==null){
      return ['mid', `Missense${where} — no language-model score available.`];
    }
    const sev = v.combined<=-6 ? 'Likely damaging' : v.combined<=-3 ? 'Possibly damaging' : 'Likely tolerated';
    const wellStruct = c.plddt>=70;
    const tone = sev==='Likely damaging' ? 'bad' : sev==='Possibly damaging' ? 'mid' : 'ok';
    return [tone, `${sev}: ${wellStruct?'well-structured':'low-confidence'} ${c.ssLabel}${where}.`];
  }

  /* =================== RENDER =================== */
  function searchBar(){
    return `<div class="card pad" style="margin-bottom:16px;display:flex;gap:10px;align-items:center;flex-wrap:wrap">
      <span style="font-size:10.5px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.4px">Gene model</span>
      <input id="foldGeneInput" value="${FD.gene||''}" placeholder="e.g. Zm00001eb406050" spellcheck="false"
        style="flex:1;min-width:280px;border:1px solid var(--line);border-radius:9px;padding:9px 11px;font-family:var(--mono);font-size:13px"
        onkeydown="if(event.key==='Enter')FOLD.loadGene()">
      <button class="btn" onclick="FOLD.loadGene()">Load structure</button>
    </div>`;
  }

  /* ---------- persistent container ----------
     SNPFold's content is expensive to render (structure fetch + 3D viewer), so we
     build it once into FD.root and re-attach that same node whenever the tool is
     shown again — instead of rebuilding from scratch. Because we hold a JS reference
     to FD.root, it survives another tool overwriting #page (it is simply detached,
     not destroyed), which lets the 3D view stay warm across navigation. */
  function ensureRoot(page){
    if (!FD.root){
      FD.root = document.createElement('div');
      FD.root.className = 'snpfold-root';
    }
    if (FD.root.parentNode !== page){
      page.innerHTML = '';
      page.appendChild(FD.root);
    }
  }

  /* Entry point called by the suite shell on navigation.
     - Arriving from another tool (SNPVersity sets S.foldGene) → autoload that gene.
     - Arriving from the side-panel menu → do NOT autoload; show a landing page the
       first time, and simply re-show the already-loaded page on later visits. */
  async function render(page){
    page = page || document.getElementById('page');
    injectCSS();
    ensureRoot(page);

    /* keep FD.dataset in sync with the rest of the app (side panel / other tools) */
    const activeDataset = resolveDataset();
    if (activeDataset != null) FD.dataset = activeDataset;
    FD.sec = hasSecondaryScores(FD.dataset);

    const requestedGene = (typeof S !== 'undefined' && S && S.foldGene) ? S.foldGene : null;

    if (requestedGene){
      /* came from another tool — honor the explicit request and autoload */
      S.foldGene = null;
      if (requestedGene !== FD.loadedGene || !FD.loaded){
        FD.gene = requestedGene;
        FD.selId = null;
        FD.openCarrier = null;
        await loadStructure();
      } else {
        /* same gene already loaded — just re-show it */
        reshowLoaded();
      }
      return;
    }

    /* came from the side-panel menu */
    if (FD.loaded){
      reshowLoaded();            // preserve the rendered page, don't reload
      return;
    }
    renderLanding();             // first visit from the menu: no autoload
  }

  /* re-show the already-rendered content without refetching */
  function reshowLoaded(){
    syncDatasetChooser();
    if (FD.struct && FD.pdb){
      if (!FD.viewer) buildViewer();           // (re)init if it never got a live DOM
      else { try { FD.viewer.resize(); FD.viewer.render(); } catch (e){} }
    }
    if (typeof attachTT==='function') attachTT();
  }

  /* Landing shown when the tool is opened from the menu — dataset chooser + gene
     search, but no structure/variant fetch until the user asks for it. */
  function renderLanding(){
    FD.root.innerHTML = datasetChooser() + searchBar() + `
      <div class="empty-state"><div class="ei">${ICONS.fold}</div>
        <h3>Structure-aware variant interpretation</h3>
        <p>Pick a dataset above and a gene model, then press <b>Load structure</b> to map
        coding variants onto the predicted protein — a linear browser, an interactive 3D
        model, and per-variant structural context.</p>
        <div style="margin-top:14px"><button class="btn primary" onclick="FOLD.loadGene()">Load structure for ${escFold(FD.gene)}</button></div>
      </div>`;
    if (typeof attachTT==='function') attachTT();
  }

  /* The heavy path: fetch the model + variants and render the full SNPFold UI into
     FD.root. Used both for autoload-from-tool and for explicit user loads. */
  async function loadStructure(){
    FD.viewer = null;
    FD.loaded = true;               // commit to the loaded view (keep across navigation)
    FD.loadedGene = FD.gene;

    const header = datasetChooser() + searchBar();

    FD.root.innerHTML = header + `<div class="loading" style="padding:48px;text-align:center">
      <div class="spinner"></div><div>Loading model for ${FD.gene}…</div></div>`;
    try { await Data.ensureStructure(FD.gene); } catch (e){ /* no file for this gene */ }

    FD.struct = Data.structureFor(FD.gene);
    FD.pdb    = Data.pdbFor(FD.gene);

    if (!FD.struct){
      FD.root.innerHTML = datasetChooser() + searchBar() + `<div class="empty-state"><div class="ei">${ICONS.fold}</div>
        <h3>No predicted model for “${FD.gene}”</h3>
        <p>SNPFold loads a protein model from <span class="c-mono">structure-&lt;gene&gt;.js</span>.
        Search another gene model above, or generate its structure file and drop it in
        <span class="c-mono">js/structures/</span>.</p></div>`;
      if (typeof attachTT==='function') attachTT();
      return;
    }
    const s = FD.struct;
    const plddtValues = Array.isArray(s.plddt) ? s.plddt.map(finiteNumber).filter(v => v != null) : [];
    const meanP = plddtValues.length
      ? (plddtValues.reduce((a,b)=>a+b,0)/plddtValues.length).toFixed(0)
      : '—';

    FD.root.innerHTML = datasetChooser() + searchBar() + `<div class="loading" style="padding:48px;text-align:center">
      <div class="spinner"></div><div>Loading coding variants for ${s.gene}…</div></div>`;
    try {
      /* Passing the dataset as a second argument is backward-compatible in JavaScript:
         older one-argument implementations simply ignore it. */
      const variantsPromise = Promise.resolve(Data.queryFoldVariants(FD.gene, FD.dataset));
      const functionPromise = typeof Data.geneFunction === 'function'
        ? Promise.resolve(Data.geneFunction(FD.gene, FD.dataset)).catch(()=>null)
        : Promise.resolve(null);

      const [variants, fn] = await Promise.all([variantsPromise, functionPromise]);
      FD.variants = variants || [];

      /* If the app state did not expose the dataset, actual returned 2026 score fields
         still enable the extra columns. */
      FD.sec = FD.sec || hasSecondaryVariantScores(FD.variants);

      FD.carriers = (fn && fn.variants)
        ? Object.fromEntries(fn.variants.map(v => [v.pos+'|'+v.ref+'|'+v.alt, v]))
        : null;
    }
    catch (e){
      console.error('SNPFold variant loading error', e);
      FD.variants = [];
      FD.carriers = null;
    }

    FD.root.innerHTML = datasetChooser() + searchBar() + `
      <div class="sec"><div class="bar"></div><div>
        <div class="n">STRUCTURE-AWARE INTERPRETATION · AlphaFold + PlantCAD/ESM</div>
        <h2>See where variants, domains, and local structure align</h2>
        <p>Coding variants mapped onto the predicted protein. Read each change against its
        Pfam domain, secondary structure, and local model confidence (pLDDT) — in a linear
        browser and in 3D — to judge whether it strikes a structured core or a flexible loop.</p>
      </div></div>

      <div class="fold-context">
        <span class="g"><b class="mono">${s.gene}</b></span>
        ${s.title?`<span class="dot">·</span><span>${s.title}</span>`:''}
        <span class="dot">·</span><span><b>${s.length}</b> aa</span>
        ${s.uniprot?`<span class="dot">·</span><span>UniProt <a href="https://www.uniprot.org/uniprotkb/${s.uniprot}" target="_blank" rel="noopener">${s.uniprot}</a></span>`:''}
        <span class="dot">·</span><span>mean pLDDT <b>${meanP}</b></span>
        <span class="dot">·</span><span><b>${FD.variants.length}</b> coding variants</span>
      </div>

      <!-- linear protein browser -->
      <div class="sec" style="margin-top:24px"><div class="bar"></div><div><h2 style="font-size:16px">Protein browser</h2>
        <p>Variants (lollipops, height = severity) over domains, secondary structure, and pLDDT. Click a variant.</p></div></div>
      <div class="card pad">
        <div class="fold-track" id="foldTrack">${trackSVG()}</div>
        <div class="fold-legend" style="display:flex;flex-wrap:wrap;gap:10px 34px;align-items:flex-start">
          <div class="lg-group" style="display:flex;flex-direction:column;align-items:flex-start;gap:4px">
            <div class="lg-title" style="font-weight:600;color:var(--ink)">Lollipop:</div>
            <span class="lg"><span class="sw" style="background:${CONS_FILL.missense}"></span>Missense</span>
            <span class="lg"><span class="sw" style="background:${CONS_FILL.lof}"></span>Loss-of-function</span>
            <span class="lg"><span class="sw" style="background:${CONS_FILL.indel}"></span>In-frame indel</span>
          </div>
          <div class="lg-group" style="display:flex;flex-direction:column;align-items:flex-start;gap:4px">
            <div class="lg-title" style="font-weight:600;color:var(--ink)">Secondary structure:</div>
            <span class="lg"><span class="sw" style="background:#d24b6a"></span>helix</span>
            <span class="lg"><span class="sw" style="background:#2f6fd0"></span>strand</span>
            <span class="lg"><span class="sw" style="background:#cbd4e1"></span>coil</span>
          </div>
          <div class="lg-group" style="display:flex;flex-direction:column;align-items:flex-start;gap:4px">
            <div class="lg-title" style="font-weight:600;color:var(--ink)">AlphaFold confidence score (pLDDT):</div>
            <span class="lg"><span class="sw" style="background:#0053d6"></span>≥90</span>
            <span class="lg"><span class="sw" style="background:#65cbf3"></span>70–90</span>
            <span class="lg"><span class="sw" style="background:#ffdb13"></span>50–70</span>
            <span class="lg"><span class="sw" style="background:#ff7d45"></span>&lt;50</span>
          </div>
        </div>
      </div>

      <!-- 3D + context -->
      <div class="fold-main">
        <div class="card" style="overflow:hidden">
          <div class="fold-toolbar">
            <div class="seg">
              <button class="seg-b ${FD.colorMode==='plddt'?'on':''}"  onclick="FOLD.color('plddt')">Confidence</button>
              <button class="seg-b ${FD.colorMode==='domain'?'on':''}" onclick="FOLD.color('domain')">Domain</button>
              <button class="seg-b ${FD.colorMode==='impact'?'on':''}" onclick="FOLD.color('impact')">Variant impact</button>
            </div>
            <label class="chk"><input type="checkbox" ${FD.showVar?'checked':''} onchange="FOLD.toggleVar(this.checked)"> Variant residues</label>
            <button class="btn ghost" onclick="FOLD.reset()">Reset view</button>
          </div>
          <div id="fold3d" class="fold-3d"></div>
        </div>
        <div class="card pad fold-ctx" id="foldCtx">${ctxHTML()}</div>
      </div>

      <!-- variant table -->
      <div class="sec" style="margin-top:24px"><div class="bar"></div><div><h2 style="font-size:16px">Coding variants on this gene</h2></div></div>
      <div class="tbl-wrap" style="max-height:none"><table class="vcf fold-table">
        <thead id="foldTableHead">${tableHeadHTML()}</thead>
        <tbody id="foldTableBody">${sortedVariants().map(rowHTML).join('')}</tbody>
      </table></div>
    `;
    buildViewer();
    if (typeof attachTT==='function') attachTT();
  }

  /* ---------- linear track ---------- */
  function x(r){ return ((r-1)/(FD.struct.length-1))*W; }
  function rleBands(arr, fn){
    const out=[]; let st=0;
    for (let i=1;i<=arr.length;i++){
      if (i===arr.length || fn(arr[i])!==fn(arr[st])){ out.push({from:st+1,to:i,key:fn(arr[st])}); st=i; }
    }
    return out;
  }
  function ssRuns(ss){
    const out=[]; let st=0;
    for (let i=1;i<=ss.length;i++){
      if (i===ss.length || ss[i]!==ss[st]){ out.push({from:st+1,to:i,t:ss[st]}); st=i; }
    }
    return out;
  }
  function trackSVG(){
    const s=FD.struct, N=s.length;
    const H=176, lolliTop=8, lolliH=58, base=lolliTop+lolliH;      // lollipop baseline
    const domY=base+6, domH=18, ssY=domY+domH+22, ssH=12, pY=ssY+ssH+10, pH=12;
    const maxSev=10;
    let g='';

    // pLDDT strip (RLE by band)
    rleBands(s.plddt, b=>plddtHex(b)).forEach(seg=>{
      const x1=x(seg.from), x2=x(seg.to);
      g+=`<rect x="${x1.toFixed(1)}" y="${pY}" width="${Math.max(.6,x2-x1).toFixed(1)}" height="${pH}" fill="${seg.key}"/>`;
    });
    g+=`<text x="0" y="${pY+pH+11}" class="tlab">pLDDT confidence</text>`;

    // secondary structure
    ssRuns(s.ss).forEach(r=>{
      const x1=x(r.from), x2=x(r.to), w=Math.max(.6,x2-x1);
      if (r.t==='H')      g+=`<rect x="${x1.toFixed(1)}" y="${ssY}" width="${w.toFixed(1)}" height="${ssH}" rx="3" fill="#d24b6a"/>`;
      else if (r.t==='E') g+=`<rect x="${x1.toFixed(1)}" y="${ssY+2}" width="${w.toFixed(1)}" height="${ssH-4}" fill="#2f6fd0"/>`;
      else                g+=`<rect x="${x1.toFixed(1)}" y="${ssY+ssH/2-1}" width="${w.toFixed(1)}" height="2" fill="#cbd4e1"/>`;
    });
    g+=`<text x="0" y="${ssY-4}" class="tlab">secondary structure</text>`;

    // domains
    (s.domains||[]).forEach((d,i)=>{
      const x1=x(d.start), x2=x(d.end), w=x2-x1, mid=x1+w/2;
      if (d.kind==='region'){
        g+=`<rect x="${x1.toFixed(1)}" y="${domY+3}" width="${w.toFixed(1)}" height="${domH-6}" rx="3" fill="#eef1f6" stroke="#d3dae6" stroke-dasharray="3 2"/>`;
        g+=`<text x="${mid.toFixed(1)}" y="${domY+domH-4}" class="dlab" text-anchor="middle">${d.name}</text>`;
      } else {
        g+=`<rect x="${x1.toFixed(1)}" y="${domY}" width="${w.toFixed(1)}" height="${domH}" rx="5" fill="${DOM_FILL[i%DOM_FILL.length]}" stroke="#bcd0f5"/>`;
        g+=`<text x="${mid.toFixed(1)}" y="${domY+domH-5}" class="dlab" text-anchor="middle">${d.name} · ${d.pfam}</text>`;
      }
    });

    // ruler ticks
    for (let r=1;r<=N;r+=100){ const xx=x(r);
      g+=`<line x1="${xx.toFixed(1)}" y1="${base}" x2="${xx.toFixed(1)}" y2="${base+3}" stroke="#aab4c4"/>`;
      g+=`<text x="${xx.toFixed(1)}" y="${base-3}" class="rlab">${r}</text>`;
    }
    g+=`<line x1="0" y1="${base}" x2="${W}" y2="${base}" stroke="#dde3ec"/>`;

    // selection guide
    const sel = FD.variants.find(v=>v.id===FD.selId);
    if (sel && sel.resi>=1 && sel.resi<=N){ const xx=x(sel.resi); g+=`<line x1="${xx.toFixed(1)}" y1="${lolliTop}" x2="${xx.toFixed(1)}" y2="${pY+pH}" stroke="#13264a" stroke-dasharray="2 2" opacity=".5"/>`; }

    // lollipops
    FD.variants.forEach(v=>{
      if (v.resi<1 || v.resi>N) return;   // can't place variants outside the model
      const xx=x(v.resi);
      const sev=v.combined!=null?Math.min(maxSev, Math.abs(v.combined)):maxSev*0.5;
      const head=base - 6 - (sev/maxSev)*(lolliH-10);
      const col=CONS_FILL[v.consClass]||'#2f5bbf';
      const on=v.id===FD.selId;
      g+=`<g class="lolli ${on?'on':''}" onclick="FOLD.select('${v.id}')" data-tt="${v.variant} · ${v.consequence} · res ${v.resi}">`;
      g+=`<line x1="${xx.toFixed(1)}" y1="${base}" x2="${xx.toFixed(1)}" y2="${head.toFixed(1)}" stroke="${col}" stroke-width="${on?2:1.4}"/>`;
      g+=`<circle cx="${xx.toFixed(1)}" cy="${head.toFixed(1)}" r="${on?6:4.5}" fill="${col}" stroke="#fff" stroke-width="1.5"/>`;
      if (v.consClass==='lof') g+=`<text x="${xx.toFixed(1)}" y="${(head-8).toFixed(1)}" class="vlab" text-anchor="middle">✱</text>`;
      g+=`</g>`;
    });

    return `<svg viewBox="0 0 ${W} ${H}" class="track-svg" preserveAspectRatio="xMinYMin meet">${g}</svg>`;
  }

  /* ---------- PanEffect jump ---------- */
  function isMissense(v){ return v && v.consClass==='missense' && v.resi; }
  /* internal view switch — highlights the substitution when missense,
     otherwise just opens PanEffect on the gene */
  function panEffect(v){
    if (!FD.gene) return;
    if (typeof goPanEffect !== 'function') return go('paneffect');
    if (isMissense(v)) goPanEffect(FD.gene, {variant:{pos:v.resi, wt:v.ref||'', sub:v.alt||''}});
    else goPanEffect(FD.gene);
  }
  function peJump(v){
    if (!isMissense(v) || !FD.gene) return '';
    const sub = `${v.ref||''}${v.resi}${v.alt||''}`;
    return ` <a class="pe-jump" href="#" title="View ${escFold(sub)} in PanEffect"
      onclick="event.stopPropagation();FOLD.panEffect('${v.id}');return false;">effects ↗</a>`;
  }

  /* ---------- variant table: sortable column model ----------
     One entry per <th>, in display order. `get` returns the value the column is
     sorted on (null/undefined => always sorted to the bottom, either direction).
     `sec` marks the MaizeGDB-2026-only columns, so the header and the row markup
     stay in sync automatically. `desc1` = first click sorts high→low, which reads
     better for counts/ranks; everything else starts low→high. */
  const PRIO_RANK = { high:3, moderate:2, medium:2, low:1, modifier:0 };
  const FOLD_COLS = [
    { key:'variant',     label:'Variant',     type:'str',
      get:v => v.variant },
    { key:'consequence', label:'Consequence', type:'str',
      get:v => v.consequence },
    { key:'resi',        label:'Residue',     type:'num', num:true,
      get:v => finiteNumber(v.resi) },
    { key:'domain',      label:'Domain',      type:'str',
      get:v => { const d = ctxFor(v).domain; return d ? d.name : null; } },
    { key:'plddt',       label:'Local pLDDT', type:'num', num:true, desc1:true,
      get:v => finiteNumber(ctxFor(v).plddt) },
    { key:'ss',          label:'Structure',   type:'str',
      get:v => { const c = ctxFor(v); return c.inModel ? c.ssLabel : null; } },
    { key:'plantcad',    label:'PlantCAD',    type:'num', num:true,
      get:v => modelScore(v, 'plantcad') },
    { key:'plantcad2',   label:'PlantCAD2',   type:'num', num:true, sec:true,
      get:v => modelScore(v, 'plantcad2') },
    { key:'esm',         label:'ESM1',        type:'num', num:true,
      get:v => modelScore(v, 'esm') },
    { key:'esm2',        label:'ESM2',        type:'num', num:true, sec:true,
      get:v => modelScore(v, 'esm2') },
    { key:'esm3',        label:'ESM3',        type:'num', num:true, sec:true,
      get:v => modelScore(v, 'esm3') },
    { key:'priority',    label:'Priority',    type:'num', desc1:true,
      get:v => { const p = v.priority ? PRIO_RANK[String(v.priority).toLowerCase()] : null;
                 return p == null ? null : p; } },
    { key:'carriers',    label:'Carriers',    type:'num', num:true, desc1:true,
      get:v => { const c = carrierOf(v);
                 if (!c) return null;
                 const n = (Number(c.hom) || 0) + (Number(c.het) || 0);
                 return n === 0 ? null : n; } },
  ];
  function foldVisibleCols(){ return FOLD_COLS.filter(c => !c.sec || FD.sec); }
  function foldCol(key){ return FOLD_COLS.find(c => c.key === key) || null; }

  /* Stable sort: ties (and blanks) keep their original order, so repeated sorts
     never shuffle rows arbitrarily. */
  function sortedVariants(){
    const list = FD.variants || [];
    const col  = foldCol(FD.sort.key);
    if (!col) return list.slice();
    const dir = FD.sort.dir === 'desc' ? -1 : 1;
    return list
      .map((v, i) => ({ v, i, k: col.get(v) }))
      .sort((a, b) => {
        const aNull = a.k == null || a.k === '';
        const bNull = b.k == null || b.k === '';
        if (aNull && bNull) return a.i - b.i;
        if (aNull) return 1;          // blanks always last
        if (bNull) return -1;
        const r = col.type === 'num'
          ? (Number(a.k) - Number(b.k))
          : String(a.k).localeCompare(String(b.k), undefined, { numeric:true, sensitivity:'base' });
        return r ? r * dir : a.i - b.i;
      })
      .map(x => x.v);
  }

  function tableHeadHTML(){
    return '<tr>' + foldVisibleCols().map(c => {
      const on = FD.sort.key === c.key;
      const arrow = on ? (FD.sort.dir === 'desc' ? '▼' : '▲') : '↕';
      const tip = on
        ? `Sorted ${FD.sort.dir === 'desc' ? 'high to low' : 'low to high'} — click to reverse`
        : `Sort by ${c.label}`;
      return `<th class="fold-th${c.num ? ' num' : ''}${on ? ' sorted' : ''}" title="${escFold(tip)}"
        aria-sort="${on ? (FD.sort.dir === 'desc' ? 'descending' : 'ascending') : 'none'}"
        onclick="FOLD.sortBy('${c.key}')"><span class="fold-th-in">${escFold(c.label)}<span class="fold-ar">${arrow}</span></span></th>`;
    }).join('') + '</tr>';
  }

  /* ---------- variant table ---------- */
  function rowHTML(v){
    const c = ctxFor(v);
    const on = v.id===FD.selId;
    const cr = carrierOf(v);
    const openC = FD.openCarrier===v.id;
    return `<tr class="fold-row ${on?'sel':''}" onclick="FOLD.select('${v.id}')">
      <td class="c-mono c-alt" style="padding-left:11px">${v.variant}</td>
      <td><span class="cons ${v.consClass}">${v.consequence}</span>${peJump(v)}</td>
      <td class="num">${v.resi}</td>
      <td>${c.domain ? (c.domain.kind==='domain'?`<span class="dom-tag">${c.domain.name}</span>`:`<span style="color:var(--muted);font-size:11px">${c.domain.name}</span>`) : '<span style="color:var(--faint)">—</span>'}</td>
      <td class="num">${c.plddt==null?'<span style="color:var(--faint)">—</span>':`<span class="plddt-chip" style="background:${plddtHex(c.plddt)};color:${c.plddt>=70?'#06294f':'#5c3a06'}">${c.plddt.toFixed(0)}</span>`}</td>
      <td>${c.inModel?`<span class="ss-chip ss-${c.ss}">${c.ssLabel}</span>`:'<span style="color:var(--faint)">—</span>'}</td>
      <td class="num">${scoreCell(modelScore(v, 'plantcad'))}</td>
      ${FD.sec?`<td class="num">${scoreCell(modelScore(v, 'plantcad2'))}</td>`:''}
      <td class="num">${scoreCell(modelScore(v, 'esm'))}</td>
      ${FD.sec?`<td class="num">${scoreCell(modelScore(v, 'esm2'))}</td><td class="num">${scoreCell(modelScore(v, 'esm3'))}</td>`:''}
      <td>${v.priority?`<span class="prio ${v.priority.toLowerCase()}">${v.priority}</span>`:'<span style="color:var(--faint)">—</span>'}</td>
      <td style="text-align:center">${carrierBtn(v, cr, openC)}</td>
    </tr>${openC?carrierRow(v, cr):''}`;
  }
  function foldCols(){ return foldVisibleCols().length; }   // colspan for the expanded carrier row
  function carrierBtn(v, cr, open){
    if (!cr || (cr.hom===0 && cr.het===0)) return '<span style="color:var(--faint)">—</span>';
    return `<button class="fold-cbtn ${open?'on':''}" title="Show carrier accessions (whole panel)"
      onclick="event.stopPropagation();FOLD.carriers('${v.id}')"><span class="cc hom">${cr.hom}</span><span class="cc het">${cr.het}</span></button>`;
  }
  function carrierRow(v, cr){
    if (!cr) return `<tr class="fn-carriers"><td colspan="${foldCols()}"><div class="fn-cwrap muted">No carrier data for this variant.</div></td></tr>`;
    const carriersHom = Array.isArray(cr.carriersHom) ? cr.carriersHom : [];
    const carriersHet = Array.isArray(cr.carriersHet) ? cr.carriersHet : [];
    const chip = (id,cls)=>`<span class="carrier ${cls}">${escFold(id)}</span>`;
    const homs = carriersHom.slice(0,60).map(id=>chip(id,'hom')).join('');
    const hets = carriersHet.slice(0,60).map(id=>chip(id,'het')).join('');
    return `<tr class="fn-carriers"><td colspan="${foldCols()}">
      <div class="fn-cwrap">
        <div><div class="fn-k">Homozygous ${v.consClass==='lof'?'(candidate knockouts)':''} · ${carriersHom.length}</div>
          <div class="fn-chips">${homs||'<span class="muted">none</span>'}${carriersHom.length>60?` <span class="muted">+${carriersHom.length-60} more</span>`:''}</div></div>
        <div style="margin-top:8px"><div class="fn-k">Heterozygous · ${carriersHet.length}</div>
          <div class="fn-chips">${hets||'<span class="muted">none</span>'}${carriersHet.length>60?` <span class="muted">+${carriersHet.length-60} more</span>`:''}</div></div>
      </div></td></tr>`;
  }

  /* ---------- context panel ---------- */
  function ctxHTML(){
    const v = FD.variants.find(x=>x.id===FD.selId);
    if (!v) return `<div class="ctx-empty">
        <div class="ei2">${ICONS.fold}</div>
        <div><b>Select a variant</b><div class="muted">Click a lollipop, a table row, or a residue in 3D to see its structural context.</div></div>
      </div>`;
    const c = ctxFor(v);
    const [tone, msg] = interpret(v, c);
    return `
      <div class="ctx-head">
        <div>
          <div class="ctx-v mono">${v.variant}</div>
          <div><span class="cons ${v.consClass}">${v.consequence}</span> ${v.priority?`<span class="prio ${v.priority.toLowerCase()}">${v.priority}</span>`:''}</div>
        </div>
        <button class="btn ghost" onclick="FOLD.focus()">Focus in 3D</button>
      </div>
      <div class="interp ${tone}">${msg}</div>
      <div class="ctx-grid">
        <div class="ck"><div class="kk">Residue</div><div class="vv mono">${v.ref||''}${v.resi}${v.consClass==='missense'&&v.alt?v.alt:''}</div></div>
        <div class="ck"><div class="kk">Domain</div><div class="vv">${c.domain?c.domain.name:'—'}${c.domain&&c.domain.pfam&&c.domain.pfam!=='region'?` <span class="muted mono">${c.domain.pfam}</span>`:''}</div></div>
        <div class="ck"><div class="kk">Local confidence</div><div class="vv">${c.plddt==null?'<span class="muted">outside model</span>':`<span class="plddt-chip" style="background:${plddtHex(c.plddt)};color:${c.plddt>=70?'#06294f':'#5c3a06'}">${c.plddt.toFixed(0)}</span> <span class="muted">${plddtBand(c.plddt)}</span>`}</div></div>
        <div class="ck"><div class="kk">Secondary structure</div><div class="vv"><span class="ss-chip ss-${c.ss}">${c.ssLabel}</span></div></div>
        <div class="ck"><div class="kk">AI scores</div><div class="vv mono">${aiScoreSummary(v)}</div></div>
      </div>
      <div class="ctx-actions">
        ${v.consClass==='missense' && v.resi ? `<button class="btn" onclick="FOLD.panEffect('${v.id}')">${ICONS.effect||ICONS.star} PanEffect</button>` : ''}
        <button class="btn" onclick="go('snpimpact')">${ICONS.star} SNPImpact</button>
        <button class="btn" onclick="go('snpcompare')">${ICONS.compare} Send to SNPCompare</button>
      </div>`;
  }

  /* =================== 3D VIEWER =================== */
  function ensure3Dmol(cb){
    if (window.$3Dmol) return cb(true);
    if (FD.libState==='failed') return cb(false);
    const done = ()=>cb(!!window.$3Dmol);
    if (FD.libState==='loading'){ const t=setInterval(()=>{ if(window.$3Dmol){clearInterval(t);cb(true);} },120);
      setTimeout(()=>{clearInterval(t); if(!window.$3Dmol) cb(false);},6000); return; }
    FD.libState='loading';
    const sc=document.createElement('script'); sc.src=THREEDMOL_URL;
    sc.onload=()=>{ FD.libState='ready'; done(); };
    sc.onerror=()=>{ FD.libState='failed'; cb(false); };
    document.head.appendChild(sc);
  }
  function buildViewer(){
    const el = document.getElementById('fold3d'); if (!el) return;
    el.innerHTML = `<div class="v-loading"><div class="spinner"></div><div>Loading 3D viewer…</div></div>`;
    ensure3Dmol(ok=>{
      if (!ok || !FD.pdb){
        el.innerHTML = `<div class="fold-3d-fallback">
          <div class="ei2">${ICONS.fold}</div>
          <b>3D viewer unavailable offline</b>
          <p>The interactive structure needs the 3Dmol.js library (loaded from a CDN). The protein
          browser, domains, and structural context above work without it. To enable 3D on a
          restricted network, self-host <span class="mono">3Dmol-min.js</span> and point SNPFold at it.</p>
        </div>`;
        return;
      }
      try{
        el.innerHTML='';
        FD.viewer = $3Dmol.createViewer(el, { backgroundColor:'white' });
        FD.viewer.addModel(FD.pdb, 'pdb');
        applyStyle();
        FD.viewer.zoomTo();
        FD.viewer.render();
        if (FD.selId) focusResidue(false);
      }catch(e){
        console.error('3Dmol error', e);
        el.innerHTML = `<div class="fold-3d-fallback"><b>Could not render structure</b><p class="mono">${String(e.message||e)}</p></div>`;
      }
    });
  }
  function applyStyle(){
    if (!FD.viewer) return;
    const v=FD.viewer; v.setStyle({}, {});
    if (FD.colorMode==='plddt'){
      v.setStyle({}, { cartoon:{ colorfunc:(a)=>plddtHex(a.b) } });
    } else if (FD.colorMode==='domain'){
      v.setStyle({}, { cartoon:{ color:'#dfe5ee' } });
      (FD.struct.domains||[]).filter(d=>d.kind==='domain').forEach((d,i)=>{
        v.setStyle({ resi:d.start+'-'+d.end }, { cartoon:{ color: ['#7fae7f','#6f97d6','#a98fd0','#d6a45b'][i%4] } });
      });
    } else {
      v.setStyle({}, { cartoon:{ color:'#d8dde6' } });
    }
    if (FD.showVar){
      FD.variants.forEach(vr=>{
        const col = FD.colorMode==='impact' ? impactHex(vr.combined) : (CONS_FILL[vr.consClass]||'#2f5bbf');
        v.addStyle({ resi:vr.resi }, { stick:{ radius:.18 } });
        v.addStyle({ resi:vr.resi }, { sphere:{ scale: vr.id===FD.selId ? 0.7 : 0.45, color: col } });
      });
    }
    v.render();
  }
  function focusResidue(animate){
    const v=FD.variants.find(x=>x.id===FD.selId); if(!FD.viewer||!v) return;
    FD.viewer.removeAllLabels();
    FD.viewer.addLabel(v.variant, { position:{ resi:v.resi }, backgroundColor:'#13264a', backgroundOpacity:.9,
      fontColor:'white', fontSize:11, borderThickness:0 }, { resi:v.resi });
    applyStyle();
    FD.viewer.zoomTo({ resi:v.resi }, animate?500:0);
    FD.viewer.render();
  }

  /* =================== handlers =================== */
  function refreshSelection(){
    const t=document.getElementById('foldTrack'); if(t) t.innerHTML=trackSVG();
    const c=document.getElementById('foldCtx'); if(c) c.innerHTML=ctxHTML();
    refreshTable();
  }
  /* re-renders header (sort arrows) and body (sorted rows) together */
  function refreshTable(){
    const h=document.getElementById('foldTableHead'); if(h) h.innerHTML=tableHeadHTML();
    const b=document.getElementById('foldTableBody'); if(b) b.innerHTML=sortedVariants().map(rowHTML).join('');
    if (typeof attachTT==='function') attachTT();
  }
  window.FOLD = {
    select(id){ FD.selId = (FD.selId===id?null:id); refreshSelection(); if(FD.selId) focusResidue(true); else if(FD.viewer){FD.viewer.removeAllLabels();applyStyle();FD.viewer.zoomTo();FD.viewer.render();} },
    carriers(id){ FD.openCarrier = (FD.openCarrier===id?null:id); refreshTable(); },
    /* Click a column header to sort; click the same header again to reverse.
       Blank cells (—) always sort to the bottom, whichever direction is active. */
    sortBy(key){
      const col = foldCol(key); if (!col) return;
      if (FD.sort.key === key) FD.sort = { key, dir: FD.sort.dir === 'asc' ? 'desc' : 'asc' };
      else                     FD.sort = { key, dir: col.desc1 ? 'desc' : 'asc' };
      refreshTable();
    },
    panEffect(id){ panEffect(FD.variants.find(v=>v.id===id)); },
    color(m){ FD.colorMode=m; document.querySelectorAll('.fold-toolbar .seg-b').forEach(b=>b.classList.remove('on'));
      const map={plddt:0,domain:1,impact:2}; const btns=document.querySelectorAll('.fold-toolbar .seg-b'); if(btns[map[m]])btns[map[m]].classList.add('on'); applyStyle(); },
    toggleVar(on){ FD.showVar=on; applyStyle(); },
    reset(){ if(FD.viewer){ FD.viewer.zoomTo(); FD.viewer.render(); } },
    focus(){ focusResidue(true); },
    loadGene(){ const el=document.getElementById('foldGeneInput'); if(!el)return;
      const g=el.value.trim(); if(!g)return; FD.gene=g; FD.selId=null; FD.openCarrier=null; loadStructure(); },
    setDataset(dataset){ FD.dataset=dataset; if(typeof S!=='undefined'&&S)S.dataset=datasetId(dataset); FD.sec=hasSecondaryScores(dataset);
      if(FD.loaded) loadStructure(); else renderLanding(); },
    /* compact dataset chooser on this page — select which dataset to use directly */
    pickDataset(id){
      const ds = foldDatasets().find(d=>String(datasetId(d))===String(id));
      const val = ds ? datasetId(ds) : id;
      if (isCurrentDataset(ds || val)) { syncDatasetChooser(); return; }   // no change
      FD.dataset = val;
      if (typeof S !== 'undefined' && S) S.dataset = val;   // keep the whole app in sync
      FD.sec = hasSecondaryScores(val);
      /* Selecting a dataset only records the choice + moves the highlight. The page
         is NOT reloaded here — the new dataset is applied when the user presses
         "Load structure" (FOLD.loadGene), which reads the current FD.dataset. */
      syncDatasetChooser();
    },
  };

  /* carrier-chip + expandable-row styles — shared look with SNPFunction */
  function injectCSS(){
    if (document.getElementById('snpfold-carrier-css')) return;
    const s=document.createElement('style'); s.id='snpfold-carrier-css';
    s.textContent=`
      /* sortable variant-table headers
         Hover/active states are translucent overlays rather than opaque fills, so they
         tint whatever the .vcf header background already is and never fight its text
         colour. On a light header swap the two rgba(255,255,255,…) values for
         rgba(0,0,0,.05) / rgba(0,0,0,.08) to darken instead of lighten. */
      .fold-table th.fold-th{cursor:pointer;user-select:none;white-space:nowrap;position:relative;
        transition:background-color .12s}
      .fold-table th.fold-th:hover{background-color:rgba(55,55,55,.33)}
      .fold-table th.fold-th .fold-th-in{display:inline-flex;align-items:center;gap:5px}
      .fold-table th.fold-th.num .fold-th-in{justify-content:flex-end}
      .fold-table th.fold-th .fold-ar{font-size:9px;line-height:1;opacity:.35;transition:opacity .12s}
      .fold-table th.fold-th:hover .fold-ar{opacity:.7}
      .fold-table th.fold-th.sorted{background-color:rgba(100,100,100,.75)}
      .fold-table th.fold-th.sorted .fold-ar{opacity:1}
      .fold-cbtn{display:inline-flex;align-items:center;border:1px solid var(--line);border-radius:7px;overflow:hidden;background:#fff;cursor:pointer;padding:0}
      .fold-cbtn:hover{border-color:#c3cee0}
      .fold-cbtn.on{border-color:#9db4dd;box-shadow:0 0 0 2px rgba(47,106,208,.12)}
      .fold-cbtn .cc{font-family:var(--mono);font-size:11px;font-weight:600;padding:2px 7px}
      .fold-cbtn .cc.hom{background:#fdecea;color:#8f281c}
      .fold-cbtn .cc.het{background:#eef4ff;color:#274b8f}
      .fn-carriers td{background:#fbfcfe;border-bottom:1px solid var(--line)}
      .fn-cwrap{padding:10px 12px} .fn-cwrap.muted{color:var(--muted)}
      .fn-chips{display:flex;flex-wrap:wrap;gap:5px;margin-top:5px}
      .fn-k{font-size:10.5px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.4px;margin-bottom:4px}
      .carrier{font-family:var(--mono);font-size:11px;padding:2px 7px;border-radius:6px;border:1px solid var(--line)}
      .carrier.hom{background:#fdecea;border-color:#f0c4bd;color:#8f281c}
      .carrier.het{background:#eef4ff;border-color:#cfe0ff;color:#274b8f}
      /* compact dataset chooser (shares Data.datasets() with SNPVersity) */
      .fold-ds-card{padding:12px 14px}
      .fold-ds-head{font-size:10.5px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.4px;margin-bottom:9px}
      .fold-ds-grid{display:flex;flex-wrap:wrap;gap:10px}
      .fold-ds{display:flex;align-items:center;gap:9px;text-align:left;cursor:pointer;background:#fff;
        border:1px solid var(--line);border-radius:10px;padding:9px 13px;min-width:210px;transition:border-color .12s,box-shadow .12s}
      .fold-ds:hover{border-color:#c3cee0}
      .fold-ds.sel{border-color:#9db4dd;box-shadow:0 0 0 2px rgba(47,106,208,.12);background:#f7faff}
      .fold-ds-dot{flex:0 0 auto;width:10px;height:10px;border-radius:50%;border:2px solid #c3cee0;background:#fff}
      .fold-ds.sel .fold-ds-dot{border-color:#2f6ad0;background:#2f6ad0;box-shadow:inset 0 0 0 2px #fff}
      .fold-ds-txt{display:flex;flex-direction:column;line-height:1.25;min-width:0}
      .fold-ds-name{font-weight:600;font-size:13px;color:var(--ink)}
      .fold-ds-sub{font-size:11.5px;color:var(--muted)}`;
    document.head.appendChild(s);
  }

  window.addEventListener('resize', ()=>{ if(FD.viewer) try{ FD.viewer.resize(); }catch(e){} });

  SNPTools.register('snpfold', { render });
})();
