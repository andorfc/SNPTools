/* =====================================================================
 *  snpversity.js — Visualization & Search tool.
 *  Registers \u2018snpversity\u2019. All render fns are global so the inline
 *  onclick/onchange handlers in the markup resolve against them.
 * ===================================================================== */

/* reference data + chromosome geometry come from the data layer */
const DATASETS    = Data.datasets();
/* PROJECTS + ACCESSIONS depend on the chosen dataset (each .h5 family has
   its own real accession columns), so they are reassigned in selectDataset. */
let   PROJECTS    = Data.projectsFor(S.dataset);
let   ACCESSIONS  = Data.accessionsFor(S.dataset);
const GENE_MODELS = Data.geneModels();
const CHR_LEN     = Data.chromLengths();
const CENTRO      = Data.centromeres();

/* friendly default: preselect up to 12 founders (one run each) for this dataset */
Data.defaultSelectionFor(S.dataset).forEach(id=>S.selected.add(id));

/* default query region */
S.chr='chr10'; S.start=3750832; S.end=3755732;

/* ================= SNPVERSITY PAGE ================= */
/* Other tools (SNPTrait, SNPGeo, SNPFold dataset pills) set S.dataset directly.
   Rebind this page's catalogue when that happened, otherwise an inbound
   hand-off is resolved against the previous dataset's accessions and every
   id is reported missing. Data.accessionsFor() caches one array per family,
   so an identity check is enough. */
function syncDatasetBinding(){
  const acc = Data.accessionsFor(S.dataset);
  if (ACCESSIONS !== acc){
    PROJECTS = Data.projectsFor(S.dataset);
    ACCESSIONS = acc;
    const known = new Set(acc.map(a => a.id));
    [...S.selected].forEach(id => { if (!known.has(id)) S.selected.delete(id); });
  }
}
function renderVersity(){
  injectVersityCSS();
  syncDatasetBinding();
  const inbound = applyPendingRequest();   // e.g. carriers handed over from SNPFunction
  const p=document.getElementById('page');
  p.className='page fade';
  p.innerHTML = `
    ${inbound ? inboundBanner(inbound) : ''}
    <div class="sec"><div class="bar"></div><div style="width:100%">
      <div class="n">VCF BUILDER & VIEWER · B73 v5</div>
      <h2>Build a variant view across maize accessions</h2>
      <p>Choose a dataset, set a genomic interval, and pick the accessions you want. SNPVersity returns a color-coded variant table and a downloadable VCF — alleles, effects, and DNA/protein language-model scores included.</p>
    </div></div>

    <!-- DATASET -->
    <div class="sec"><div class="bar"></div><div><h2 style="font-size:16px">1 · Choose a dataset</h2></div></div>
    <div class="ds-grid" id="dsGrid"></div>

    <!-- REGION -->
    <div class="sec"><div class="bar"></div><div><h2 style="font-size:16px">2 · Select a genomic interval</h2></div></div>
    <div class="card region-card" id="regionCard"></div>

    <!-- ACCESSIONS -->
    <div class="sec"><div class="bar"></div><div><h2 style="font-size:16px">3 · Choose accessions</h2>
      <p>Browse and filter the lines by panel, origin, subpopulation and GRIN traits, take a quick pick, or paste a list.</p></div></div>
    <div class="card acc-card" id="accCard"></div>

    <!-- RUN -->
    <div class="runbar" id="runbar"></div>

    <!-- RESULTS -->
    <div id="resultsAnchor" style="margin-top:30px"></div>
  `;
  renderDatasets(); renderRegion(); renderAccPicker(); renderRunbar();
  if(inbound){
    // prefill the gene box so "Load" re-fetches the same coordinates
    const gi=document.getElementById('geneInput');
    if(gi && inbound.gene){ gi.value=inbound.gene; }
    const st=document.getElementById('geneStatus');
    if(st && inbound.gene){ st.className='status ok'; st.textContent=`Region set from ${inbound.gene} (${S.chr}:${S.start.toLocaleString()}–${S.end.toLocaleString()})`; }
    if(inbound.missing && inbound.missing.length){
      uplSay('warn', `${inbound.missing.length} carrier${inbound.missing.length>1?'s are':' is'} not present in this dataset's accession list and could not be selected.`);
      renderUplReport({matched:inbound.added+inbound.already, considered:inbound.requested,
        unmatched:inbound.missing.map((t,i)=>({line:i+1, text:t})), dupes:[], fanout:[],
        source:inbound.from||'SNPFunction', skippedHeader:null, applied:true});
    }
    document.getElementById('runbar').scrollIntoView({behavior:'smooth',block:'center'});
    return;
  }
  // returning from another tool: restore the last results instead of forcing a rebuild
  if(S.results && S.results.rows && S.results.rows.length) renderResults();
}

/* =====================================================================
 *  INBOUND HANDOFF  —  another tool (SNPFunction) hands us a gene model
 *  and a set of accessions to preselect.
 *    payload = {gene, chr, start, end, dataset, accessions:[…], from, note}
 * ===================================================================== */
window.versityRequest = function(payload){
  if(!payload) return;
  S.pendingVersity = payload;
  S.results = null;                  // the old result no longer matches the new query
  if(typeof go==='function') go('snpversity');
};

/* consume S.pendingVersity: switch dataset, set the region, resolve + select
   the accessions. Honours req.merge ('add' | 'replace'); anything else — including
   a payload from a tool that predates the flag — replaces, which is the old
   behaviour. Returns a summary for the banner, or null. */
function applyPendingRequest(){
  const req = S.pendingVersity; if(!req) return null;
  S.pendingVersity = null;

  // Accession IDs are not portable between dataset families, so a handoff that
  // switches dataset can only ever reset the selection — 'add' is meaningless.
  const crossDataset = !!(req.dataset && req.dataset!==S.dataset && DATASETS.some(d=>d.id===req.dataset));
  const merge = crossDataset ? 'replace' : (req.merge==='add' ? 'add' : 'replace');

  // Snapshot dataset + region + selection *before* touching any of them, so a
  // single Undo reverts the whole query state coherently rather than half of it.
  snapshotQuery(`the handoff from ${req.from||'another tool'}`);

  // 1 · dataset (rebinds this dataset's real accession catalog)
  if(crossDataset){
    S.dataset=req.dataset;
    PROJECTS=Data.projectsFor(S.dataset);
    ACCESSIONS=Data.accessionsFor(S.dataset);
  }

  // 2 · region  (a region has no union semantics — the newest handoff wins)
  if(req.chr){ S.chr=String(req.chr).startsWith('chr')?req.chr:'chr'+req.chr; }
  const flank=+req.flank||0;
  if(req.start!=null && req.end!=null){
    const lo=Math.floor(Math.min(+req.start,+req.end)), hi=Math.ceil(Math.max(+req.start,+req.end));
    S.start=Math.max(0,lo-flank); S.end=hi+flank;
  }

  // 3 · accessions (resolved the same way as an uploaded list, so founder
  //     names / run IDs / composite IDs all work)
  const wanted=[...new Set((req.accessions||[]).filter(Boolean).map(String))];
  const idx=buildAccIndex(); const missing=[];
  const before=new Set(S.selected);          // for the delta report
  const resolved=new Set();
  wanted.forEach(w=>{
    const hit=idxHit(idx,w) || (w.match(RUN_RE)?idxHit(idx,w.match(RUN_RE)[1]):null);
    if(hit) hit.forEach(id=>resolved.add(id)); else missing.push(w);
  });
  if(merge==='replace') S.selected.clear();
  resolved.forEach(id=>S.selected.add(id));

  let added=0, already=0, removed=0;
  resolved.forEach(id=>{ if(before.has(id)) already++; else added++; });
  if(merge==='replace') before.forEach(id=>{ if(!resolved.has(id)) removed++; });

  S.page=1;
  return {gene:req.gene||'', from:req.from||'SNPFunction', note:req.note||'',
          merge, crossDataset, requested:wanted.length,
          added, already, removed,
          selected:S.selected.size, total:S.selected.size, missing};
}
function inboundBanner(i){
  const delta = (typeof Handoff!=='undefined')
    ? Handoff.deltaText({added:i.added, already:i.already, removed:i.removed,
                         missing:i.missing.length, total:i.total})
    : `<b>${i.selected}</b> of ${i.requested} accession${i.requested===1?'':'s'} selected`;
  return `<div class="from-fn" id="inboundBanner">
    <b>From ${escAttr(i.from)}</b>
    <span>${i.gene?`<span class="mono">${escAttr(i.gene)}</span> — `:''}${escAttr(i.note||'carrier accessions')}
      ${i.merge==='add'?'added to your selection':'replaced your selection'}: ${delta}.
      ${i.crossDataset?`<br><span class="ho-hint">Dataset changed — accessions are not shared between datasets, so the previous selection could not be kept.</span>`:''}</span>
    <span class="from-fn-acts">
      ${undoChip()}
      <button class="btn" onclick="runQuery()">${ICONS.dna||''} Build VCF &amp; view</button>
    </span>
  </div>`;
}

/* =====================================================================
 *  UNDO  —  one snapshot of the whole query state (dataset, region,
 *  selection), taken immediately before a handoff or an uploaded list is
 *  applied. It is deliberately shallow: exactly one step, and it expires
 *  the moment the user touches the selection by hand or runs the query,
 *  so a stale Undo button can never sit there waiting to eat real work.
 * ===================================================================== */
function snapshotQuery(what){
  S.undoQuery = { what: what || 'the last change', dataset:S.dataset,
                  chr:S.chr, start:S.start, end:S.end, selected:[...S.selected] };
}
/* call from anything that edits the selection by hand */
function touchSelection(){
  if(!S.undoQuery) return;
  S.undoQuery = null;
  document.querySelectorAll('.ho-undo').forEach(b=>b.remove());
}
function undoChip(){
  if(!S.undoQuery) return '';
  return `<button class="btn ghost ho-undo" onclick="undoLastChange()"
    title="Put the dataset, region and selection back the way they were">Undo</button>`;
}
function undoLastChange(){
  const u=S.undoQuery; if(!u) return;
  S.undoQuery=null;
  if(u.dataset && u.dataset!==S.dataset){
    S.dataset=u.dataset;
    PROJECTS=Data.projectsFor(S.dataset);
    ACCESSIONS=Data.accessionsFor(S.dataset);
  }
  if(u.chr) S.chr=u.chr;
  if(u.start!=null) S.start=u.start;
  if(u.end!=null)   S.end=u.end;
  S.selected.clear(); (u.selected||[]).forEach(id=>S.selected.add(id));
  S.results=null;
  const b=document.getElementById('inboundBanner'); if(b) b.remove();
  renderDatasets(); renderRegion(); renderAccPicker(); renderRunbar();
  const n=S.selected.size;
  accSay(`Undone — ${escAttr(u.what)} was reverted. Selection is back to <b>${n}</b> accession${n===1?'':'s'}.`);
}

function renderDatasets(){
  document.getElementById('dsGrid').innerHTML = DATASETS.map(d=>`
    <div class="ds ${d.id===S.dataset?'sel':''}" onclick="selectDataset('${d.id}')">
      <div class="dot"></div>
      <div class="t">${d.name}</div>
      <div class="ref">${d.sub} · aligned to ${d.ref}</div>
      <div class="stats">
        <div class="stat"><div class="v">${d.acc}</div><div class="k">accessions</div></div>
        <div class="stat"><div class="v">${d.sites}</div><div class="k">variant sites</div></div>
      </div>
      <div class="badges">
        <span class="chiplet ${d.het?'on':''}">${d.het?'✓':'×'} heterozygous</span>
        <span class="chiplet ${d.indel?'on':''}">${d.indel?'✓':'×'} INDELs</span>
        <span class="chiplet ${d.impute?'on':''}">${d.impute?'✓':'×'} imputed</span>
      </div>
    </div>`).join('');
}
function selectDataset(id){
  touchSelection();
  S.dataset=id;
  // switch to this dataset's real accession catalog
  PROJECTS   = Data.projectsFor(id);
  ACCESSIONS = Data.accessionsFor(id);
  S.selected.clear();
  Data.defaultSelectionFor(id).forEach(x=>S.selected.add(x));
  renderDatasets(); renderAccPicker(); renderRunbar();
  const m=document.getElementById('mAcc'); if(m)m.textContent=S.selected.size;
}

/* ---- region card + genome ribbon ---- */
function renderRegion(){
  const len=CHR_LEN[S.chr]||3e8;
  document.getElementById('regionCard').innerHTML=`
    <div class="ribbon-pane">
      <div class="ribbon-head">
        <span class="chr" id="ribChr">Chromosome ${S.chr.replace('chr','')}</span>
        <span class="coord" id="ribCoord"></span>
      </div>
      <div class="ribbon">
        <div class="chrom-track"></div>
        <div class="centro" id="centro"></div>
        <div class="window" id="window"></div>
        <div class="ticks"><span>0</span><span>${(len/1e6).toFixed(0)} Mb</span></div>
      </div>
      <div class="region-meta">
        <div class="m"><div class="v" id="mSpan"></div><div class="k">interval span</div></div>
        <div class="m"><div class="v" id="mAcc"></div><div class="k">accessions selected</div></div>
      </div>
    </div>
    <div class="form-pane">
      <div class="field">
        <label>Chromosome</label>
        <select id="chrInput" onchange="onRegion()">
          ${Object.keys(CHR_LEN).map(c=>`<option value="${c}" ${c===S.chr?'selected':''}>Chromosome ${c.replace('chr','')}</option>`).join('')}
        </select>
      </div>
      <div class="field">
        <label>Interval (bp)</label>
        <div class="row">
          <input type="number" id="startInput" value="${S.start}" min="0" oninput="onRegion()" aria-label="start">
          <span style="align-self:center;color:var(--faint)">–</span>
          <input type="number" id="endInput" value="${S.end}" min="0" oninput="onRegion()" aria-label="end">
        </div>
      </div>
      <div class="gene-row">
        <div class="field" style="margin:0">
          <label>Gene model · auto-fill coordinates</label>
          <input type="text" id="geneInput" placeholder="Zm00001eb…" class="mono-in">
        </div>
        <button class="btn" onclick="loadGene()" style="margin-bottom:1px">Load</button>
      </div>
      <div style="display:flex;gap:14px;align-items:center">
        <button class="link-btn" onclick="exampleGene()">Use an example gene</button>
        <div class="field" style="margin:0;flex:1">
          <label>± flank (bp)</label>
          <input type="number" id="flankInput" value="0" min="0" class="mono-in">
        </div>
        <div class="field" style="margin:0;width:120px">
          <label>Loci / page</label>
          <select id="perPage" onchange="S.perPage=+this.value">
            <option>50</option><option selected>100</option><option>250</option><option>500</option>
          </select>
        </div>
      </div>
      <div class="status" id="geneStatus"></div>
    </div>`;
  drawRibbon();
}
function onRegion(){
  S.chr=document.getElementById('chrInput').value;
  S.start=+document.getElementById('startInput').value||0;
  S.end=+document.getElementById('endInput').value||0;
  drawRibbon(); renderRunbar();
  // update chr label
  document.getElementById('ribChr').textContent='Chromosome '+S.chr.replace('chr','');
}
function drawRibbon(){
  const len=CHR_LEN[S.chr]||3e8;
  const a=Math.max(0,Math.min(S.start,len)), b=Math.max(0,Math.min(S.end,len));
  const lo=Math.min(a,b), hi=Math.max(a,b);
  const w=document.getElementById('window');
  const left=(lo/len)*100, width=Math.max(((hi-lo)/len)*100,0.4);
  if(w){w.style.left=left+'%';w.style.width=width+'%';}
  const ce=document.getElementById('centro'); if(ce)ce.style.left=((CENTRO[S.chr]||.4)*100)+'%';
  const span=hi-lo;
  const fmt=n=>n>=1e6?(n/1e6).toFixed(2)+' Mb':n>=1e3?(n/1e3).toFixed(1)+' kb':n+' bp';
  document.getElementById('ribCoord').textContent=`${S.chr}:${lo.toLocaleString()}–${hi.toLocaleString()}`;
  document.getElementById('mSpan').textContent=fmt(span);
  document.getElementById('mAcc').textContent=S.selected.size;
}
async function loadGene(){
  const id=document.getElementById('geneInput').value.trim();
  const st=document.getElementById('geneStatus');
  if(!id){st.className='status err';st.textContent='Enter a gene model ID first.';return;}
  st.className='status'; st.textContent='Looking up '+id+'…';
  let g;
  try{
    g=await Data.lookupGene(id);
  }catch(err){
    st.className='status err'; st.textContent='Lookup failed: '+((err&&err.message)?err.message:err);
    console.error('Gene lookup failed:', err);
    return;
  }
  if(!g){st.className='status err';st.textContent='Gene model “'+id+'” not found.';return;}
  const flank=+document.getElementById('flankInput').value||0;
  S.chr=g.chr; S.start=Math.max(0,g.start-flank); S.end=g.end+flank;
  const chrSel=document.getElementById('chrInput');
  if(chrSel){
    if(![...chrSel.options].some(o=>o.value===g.chr)){
      chrSel.add(new Option('Chromosome '+g.chr.replace('chr',''), g.chr));
    }
    chrSel.value=g.chr;
  }
  document.getElementById('startInput').value=S.start;
  document.getElementById('endInput').value=S.end;
  document.getElementById('ribChr').textContent='Chromosome '+g.chr.replace('chr','');
  st.className='status ok'; st.textContent='Loaded '+id+' ('+g.chr+':'+g.start.toLocaleString()+'–'+g.end.toLocaleString()+')';
  drawRibbon(); renderRunbar();
}
function exampleGene(){
  const ids=Data.exampleGenes(); const id=pick(ids);
  document.getElementById('geneInput').value=id; loadGene();
}

/* ---- accession picker ---- */
function escAttr(s){return String(s==null?'':s).replace(/&/g,'&amp;').replace(/"/g,'&quot;').replace(/</g,'&lt;').replace(/>/g,'&gt;');}
/* REF/ALT sequences can be long indels; truncate the visible text to keep the
   table columns from stretching, while the full sequence stays available on
   hover via the data-tt tooltip. */
function alleleDisp(seq){
  const s=String(seq==null?'':seq);
  return s.length>10 ? s.slice(0,10)+'…' : s;
}
function alleleTT(seq,label){
  const s=String(seq==null?'':seq);
  return s.length>10 ? `${label}: ${s}` : label;
}
/* ---- step 3: the selection, quick picks, and the full selector in a drawer ----
   SNPVersity no longer carries its own browser of the catalogue (it was a list of projects of
   accession chips, a lesser copy of SNPTrait). "Browse & filter lines…" opens SNPTrait over the
   page (TraitDrawer, js/snptrait.js) on a draft of this selection; Apply writes it back. Step 3
   keeps the quick jobs: what is selected (count, panel mix, chips), one-click picks, and a
   pasted or uploaded list (below). */
let accShowAll=false;               // show every selected chip, not the first ACC_CHIPS
const ACC_CHIPS=40;
function renderAccPicker(){
  document.getElementById('accCard').innerHTML=`
    <div class="acc-body">
      <div class="acc2">
        <div class="acc2-top">
          <div class="acc2-count"><span class="acc2-n" id="selCount">0</span> <span class="acc2-of" id="selOf"></span></div>
          <button class="btn primary" id="browseLinesBtn" onclick="openLineSelector()" aria-haspopup="dialog">
            ${ICONS.search||''} Browse &amp; filter lines…</button>
        </div>
        <div class="acc2-mix" id="selMix"></div>
        <div class="sel-chips acc2-chips" id="selChips"></div>
        <div class="acc2-picks">
          <div class="acc2-k">Quick picks <span>replace the selection · Undo is offered</span></div>
          <div class="acc2-pickrow" id="accPicks"></div>
          <div class="acc2-say" id="accSay" role="status" aria-live="polite"></div>
        </div>
      </div>
      <div class="acc-side">
        <div class="upl">
          <label>Or paste / upload a list (one accession per line)</label>
          <div class="file-row">
            <input type="file" id="fileUpload" accept=".txt,.tsv,.csv,.list,text/plain" onchange="onAccFilePicked(this)">
          </div>
          <div class="upl-name" id="uplName">No file chosen — you can also paste a list below.</div>
          <textarea id="accPaste" class="upl-paste" spellcheck="false"
            placeholder="ZmG_B73&#10;PI 550473&#10;Mo17&#10;…"
            oninput="refreshUplButtons()"></textarea>
          <div class="upl-actions">
            <button class="btn primary" id="uplLoadBtn" onclick="loadAccList()" disabled>Load accessions</button>
            <button class="btn" onclick="clearAccUpload()">Clear</button>
          </div>
          <div class="upl-opt" data-ho-mount data-ho-id="uplReplace"></div>
          <div class="upl-status" id="uplStatus"></div>
          <div class="upl-report" id="uplReport"></div>
        </div>
      </div>
    </div>`;
  accShowAll=false;
  // same control, same wording, same session state as the handoff checkbox in
  // SNPFunction / SNPFold — setting it in one place sets it in all of them
  if(typeof Handoff!=='undefined') Handoff.sync();
  renderAccList(); renderSelected();
}
/* The drawer: SNPTrait on a draft of the selection. */
function openLineSelector(){
  if(typeof TraitDrawer==='undefined'){ go('snptrait'); return; }
  TraitDrawer.open({dataset:S.dataset, selected:[...S.selected], title:'Choose accessions for SNPVersity',
    onApply:ids=>{
      const before=new Set(S.selected);
      const same=ids.length===before.size && ids.every(id=>before.has(id));
      if(same){ accSay('No change to the selection.'); return; }
      snapshotQuery('the selection from Browse & filter');
      S.selected=new Set(ids);
      let added=0, removed=0;
      ids.forEach(id=>{ if(!before.has(id)) added++; });
      before.forEach(id=>{ if(!S.selected.has(id)) removed++; });
      renderAccList(); renderSelected(); renderRunbar();
      if(typeof refreshSelChip==='function') refreshSelChip(true);
      accSay(`Applied from Browse &amp; filter: ${added?`<b>+${added}</b>`:''}${added&&removed?' · ':''}${removed?`<b>−${removed}</b>`:''}`+
             ` · selection is now <b>${S.selected.size}</b>. ${undoChip()}`);
    }});
}
function accSay(html){ const el=document.getElementById('accSay'); if(el) el.innerHTML=html||''; }
/* Quick picks. A release with panel-membership flags (GRIN-linked: inNAM, inAmes282, inWiDiv;
   a line can be in several) offers each panel as a whole; otherwise each project of the
   catalogue, and one run per NAM founder where founders are tagged. */
function accPicks(){
  const A=ACCESSIONS, picks=[];
  const flags=[['inNAM','NAM founders + B73'],['inAmes282','Ames 282'],['inWiDiv','WiDiv']];
  if(A.some(a=>a.inNAM!=null || a.inAmes282!=null || a.inWiDiv!=null)){
    flags.forEach(([k,label])=>{ const ids=A.filter(a=>a[k]==='yes').map(a=>a.id); if(ids.length) picks.push({key:k, label, ids}); });
    const other=A.filter(a=>a.panel==='Other GRIN').map(a=>a.id);
    if(other.length) picks.push({key:'other', label:'Other GRIN', ids:other});
  } else {
    PROJECTS.forEach(p=>{ const ids=[]; p.groups.forEach(g=>g.accessions.forEach(a=>ids.push(a.id))); picks.push({key:'p:'+p.id, label:p.title, ids}); });
    if(Data.namFoundersFor(S.dataset).length){
      const seen=new Set(), ids=[];
      A.forEach(a=>{ if(a.namFounder && !seen.has(a.namFounder)){ seen.add(a.namFounder); ids.push(a.id); } });
      picks.push({key:'nam1', label:'One per NAM founder', ids});
    }
  }
  [.02,.05,.10,.25].forEach(f=>picks.push({key:'r'+f, label:`Random ${Math.round(f*100)}%`, random:f, n:Math.ceil(A.length*f)}));
  picks.push({key:'all', label:'All', ids:A.map(a=>a.id)});
  return picks;
}
function sameSet(ids){ return ids.length===S.selected.size && ids.every(id=>S.selected.has(id)); }
/* renderAccList (the name the upload and hand-off code call) draws the quick picks */
function renderAccList(){
  const box=document.getElementById('accPicks'); if(!box) return;
  box.innerHTML=accPicks().map(pk=>{
    const n=pk.ids ? pk.ids.length : pk.n;
    const on=pk.ids && pk.ids.length && sameSet(pk.ids);
    return `<button class="qbtn${on?' on':''}" onclick="applyPick('${pk.key}')" ${on?'aria-pressed="true"':''}>${escAttr(pk.label)} <span class="qn">${n}</span></button>`;
  }).join('')+`<button class="qbtn" onclick="applyPick('clear')" ${S.selected.size?'':'disabled'}>Clear</button>`;
}
function applyPick(key){
  let ids, label;
  if(key==='clear'){ ids=[]; label='nothing (cleared)'; }
  else {
    const pk=accPicks().find(x=>x.key===key); if(!pk) return;
    label=pk.label;
    ids = pk.random ? randomSample(ACCESSIONS.map(a=>a.id), pk.n) : pk.ids;
  }
  if(sameSet(ids)){ accSay(`The selection already is ${escAttr(label)}.`); return; }
  snapshotQuery(`the quick pick “${label}”`);
  S.selected=new Set(ids);
  renderAccList(); renderSelected(); renderRunbar();
  accSay(`Selection replaced with <b>${escAttr(label)}</b> — <b>${ids.length}</b> accession${ids.length===1?'':'s'}. ${undoChip()}`);
}
/* the selection: count, a bar of its make-up by panel / project, then the chips */
function renderSelected(){
  const arr=[...S.selected];
  const c=document.getElementById('selCount'); if(c) c.textContent=arr.length.toLocaleString();
  const of=document.getElementById('selOf'); if(of) of.textContent=`of ${ACCESSIONS.length.toLocaleString()} accessions selected`;
  const m=document.getElementById('mAcc'); if(m) m.textContent=arr.length;
  if(typeof refreshSelChip==='function') refreshSelChip();
  const byId=new Map(ACCESSIONS.map(a=>[a.id,a]));
  const mixBox=document.getElementById('selMix');
  if(mixBox){
    const mix=new Map();
    arr.forEach(id=>{ const a=byId.get(id); const k=a?(a.projTitle||a.proj):'not in this dataset';
      const cur=mix.get(k)||{n:0,color:a?a.projColor:'#999'}; cur.n++; mix.set(k,cur); });
    const sub=new Map();
    arr.forEach(id=>{ const a=byId.get(id); if(a && a.subpop) sub.set(a.subpop,(sub.get(a.subpop)||0)+1); });
    const subTop=[...sub.entries()].sort((x,y)=>y[1]-x[1]);
    mixBox.innerHTML = arr.length ? `
      <div class="acc2-bar" role="img" aria-label="Selection by panel: ${escAttr([...mix.entries()].map(([k,v])=>k+' '+v.n).join(', '))}">
        ${[...mix.entries()].map(([k,v])=>`<span style="flex:${v.n};background:${v.color}" title="${escAttr(k)}: ${v.n}"></span>`).join('')}</div>
      <div class="acc2-legend">${[...mix.entries()].map(([k,v])=>`<span><i style="background:${v.color}"></i>${escAttr(k)} <b>${v.n}</b></span>`).join('')}</div>
      ${subTop.length?`<div class="acc2-sub">Subpopulation: ${subTop.slice(0,6).map(([k,n])=>`${escAttr(k)} <b>${n}</b>`).join(' · ')}${subTop.length>6?` · +${subTop.length-6} more`:''}</div>`:''}` : '';
  }
  const box=document.getElementById('selChips'); if(!box) return;
  if(!arr.length){ box.innerHTML='<div class="empty">Nothing selected yet — browse &amp; filter the lines, take a quick pick, or paste a list.</div>'; return; }
  const show=accShowAll ? arr : arr.slice(0, ACC_CHIPS);
  box.innerHTML=show.map(id=>{ const a=byId.get(id);
      return `<span class="sel-chip" title="${escAttr(a?(a.label||a.id)+' · '+(a.projTitle||''):id)}"><span class="dotc" style="display:inline-block;width:7px;height:7px;border-radius:2px;background:${a?a.projColor:'#999'}"></span>${escAttr(a?a.run:id)}${a?sampleHetMark(a):''}<button onclick="toggleAcc('${escAttr(id)}')" aria-label="Remove ${escAttr(a?a.run:id)}">×</button></span>`;
    }).join('')
    + (arr.length>ACC_CHIPS ? `<button class="link-more" onclick="accShowAll=!accShowAll;renderSelected()">${accShowAll?'Show fewer':`+${arr.length-ACC_CHIPS} more — show all`}</button>` : '');
}
/* =====================================================================
 *  ACCESSION LIST UPLOAD  (file or paste)
 *  Reads a plain list, resolves each entry against the current dataset's
 *  accession catalog, and reports exactly what matched / what did not.
 * ===================================================================== */
let uplFileText = null;      // text of the last successfully read file
let uplFileName = '';

const MAX_UPL_BYTES = 5 * 1024 * 1024;
const RUN_RE = /\b([SEDC]RR\d{4,})\b/i;                 // SRR / ERR / DRR / CRR run accessions
const RUN_SUFFIX_RE = /[_\s,\t-]+[SEDC]RR\d{4,}\s*$/i;  // trailing "_SRR12460455"

function normKey(s){ return String(s==null?'':s).replace(/\uFEFF/g,'').trim().replace(/^["']+|["']+$/g,'').toLowerCase(); }
function looseKey(s){ return normKey(s).replace(/[^a-z0-9]+/g,''); }

/* key → Set(accession id). Every accession is indexed under its id, run,
   founder, label and the founder_run composite (both orders), in exact and
   punctuation-insensitive form. Rebuilt on demand so it follows the dataset. */
function buildAccIndex(){
  const idx=new Map();
  const add=(k,id)=>{ if(!k) return; if(!idx.has(k)) idx.set(k,new Set()); idx.get(k).add(id); };
  ACCESSIONS.forEach(a=>{
    // GRIN-linked lines also answer to their GRIN accession (PI 550473), line name and genotype id
    const keys=[a.id, a.run, a.founder, a.label, a.grin, a.strain, a.genotypeId];
    if(a.founder&&a.run){ keys.push(a.founder+'_'+a.run, a.run+'_'+a.founder); }
    keys.filter(Boolean).forEach(v=>{ add(normKey(v),a.id); add(looseKey(v),a.id); });
  });
  return idx;
}
function idxHit(idx,key){ return idx.get(normKey(key)) || idx.get(looseKey(key)) || null; }

/* resolve one line to accession id(s); returns {ok, ids, via} or {ok:false, text} */
function resolveAccEntry(raw, idx){
  const line=String(raw).replace(/\uFEFF/g,'').trim();
  if(!line || /^[#!]/.test(line)) return null;                 // blank or comment → skipped
  // whole line first, then each delimited field (handles csv/tsv exports)
  const toks=[line].concat(line.split(/[,\t;|]+/)).map(t=>t.trim().replace(/^["']+|["']+$/g,'')).filter(Boolean);
  for(const t of toks){ const hit=idxHit(idx,t); if(hit) return {ok:true, ids:[...hit], via:t, text:line}; }
  // fall back to the embedded run accession …
  const m=line.match(RUN_RE);
  if(m){ const hit=idxHit(idx,m[1]); if(hit) return {ok:true, ids:[...hit], via:m[1], text:line}; }
  // … then to the founder part with the run suffix stripped
  const f=line.replace(RUN_SUFFIX_RE,'').trim();
  if(f && f!==line){ const hit=idxHit(idx,f); if(hit) return {ok:true, ids:[...hit], via:f, text:line}; }
  return {ok:false, text:line};
}

function uplSay(cls,msg){ const el=document.getElementById('uplStatus'); if(el){ el.className='upl-status '+(cls||''); el.innerHTML=msg||''; } }
function refreshUplButtons(){
  const btn=document.getElementById('uplLoadBtn'); if(!btn) return;
  const pasted=(document.getElementById('accPaste')||{}).value||'';
  btn.disabled = !(uplFileText && uplFileText.trim()) && !pasted.trim();
}
function clearAccUpload(){
  uplFileText=null; uplFileName='';
  const f=document.getElementById('fileUpload'); if(f) f.value='';
  const t=document.getElementById('accPaste');  if(t) t.value='';
  const n=document.getElementById('uplName');   if(n) n.textContent='No file chosen — you can also paste a list below.';
  const r=document.getElementById('uplReport'); if(r) r.innerHTML='';
  uplSay('',''); refreshUplButtons();
}

/* file chosen → validate + read into memory (nothing is selected until "Load") */
function onAccFilePicked(input){
  const file=input.files && input.files[0];
  const nameEl=document.getElementById('uplName');
  const rep=document.getElementById('uplReport'); if(rep) rep.innerHTML='';
  uplFileText=null; uplFileName='';
  if(!file){ if(nameEl) nameEl.textContent='No file chosen — you can also paste a list below.'; uplSay('',''); refreshUplButtons(); return; }
  if(file.size>MAX_UPL_BYTES){
    if(nameEl) nameEl.textContent=file.name;
    uplSay('err',`That file is ${(file.size/1048576).toFixed(1)} MB. Please upload a plain list under 5 MB.`);
    refreshUplButtons(); return;
  }
  if(file.size===0){
    if(nameEl) nameEl.textContent=file.name;
    uplSay('err','That file is empty.'); refreshUplButtons(); return;
  }
  if(/\.(xlsx|xls|pdf|docx?|zip|gz|bam|h5|vcf)$/i.test(file.name)){
    if(nameEl) nameEl.textContent=file.name;
    uplSay('err','Unsupported file type. Use a plain text, .csv or .tsv list with one accession per line.');
    refreshUplButtons(); return;
  }
  const fr=new FileReader();
  fr.onerror=()=>{ uplSay('err','Could not read that file.'); refreshUplButtons(); };
  fr.onload=()=>{
    const txt=String(fr.result||'');
    if(/\u0000/.test(txt)){ uplSay('err','That looks like a binary file, not a text list.'); refreshUplButtons(); return; }
    uplFileText=txt; uplFileName=file.name;
    if(nameEl) nameEl.textContent=`${file.name} · ${(file.size/1024).toFixed(1)} KB · ${txt.split(/\r\n|\r|\n/).filter(l=>l.trim()).length} non-empty lines`;
    uplSay('ok','File read. Press <b>Load accessions</b> to match it against this dataset.');
    refreshUplButtons();
  };
  fr.readAsText(file);
}

/* main entry: parse whatever is available (file wins, else the textarea) */
function loadAccList(){
  const pasted=((document.getElementById('accPaste')||{}).value||'');
  const text = (uplFileText && uplFileText.trim()) ? uplFileText : pasted;
  const src  = (uplFileText && uplFileText.trim()) ? (uplFileName||'uploaded file') : 'pasted list';
  if(!text || !text.trim()){ uplSay('err','Nothing to load — choose a file or paste a list first.'); return; }
  applyAccList(text, src);
}

function applyAccList(text, source){
  const idx=buildAccIndex();
  const lines=String(text).split(/\r\n|\r|\n/);
  const matched=new Map();      // accession id → the entry that matched it
  const unmatched=[];           // {line, text}
  const dupes=[];               // entries that resolved to something already matched
  const fanout=[];              // entries that resolved to >1 accession (e.g. a founder with reps)
  let considered=0, skippedHeader=null;

  lines.forEach((raw,i)=>{
    const r=resolveAccEntry(raw, idx);
    if(r===null) return;                                   // blank / comment
    // tolerate a single header row at the top
    if(!r.ok && considered===0 && skippedHeader===null &&
       /^(accession|accessions|id|sample|sample_?id|name|run|taxa|line)\b/i.test(r.text)){
      skippedHeader=r.text; return;
    }
    considered++;
    if(!r.ok){ unmatched.push({line:i+1, text:r.text}); return; }
    if(r.ids.length>1) fanout.push({text:r.text, n:r.ids.length});
    r.ids.forEach(id=>{ if(matched.has(id)) dupes.push(r.text); else matched.set(id, r.text); });
  });

  if(!considered){ uplSay('err','No usable entries found — the list looks empty or contains only comments.'); return; }
  if(!matched.size){
    uplSay('err',`None of the ${considered} entries in <b>${escAttr(source)}</b> matched an accession in this dataset. `+
                 `Check that you picked the right dataset above, or that the file is one accession per line.`);
    renderUplReport({matched:0, considered, unmatched, dupes, fanout, source, skippedHeader, applied:false});
    return;
  }

  const merge = (typeof Handoff!=='undefined') ? Handoff.mode()
              : (((document.getElementById('uplReplace')||{}).checked) ? 'replace' : 'add');
  snapshotQuery(`the list loaded from ${source}`);
  const before=new Set(S.selected);
  if(merge==='replace') S.selected.clear();
  matched.forEach((_,id)=>S.selected.add(id));

  let added=0, already=0, removed=0;
  matched.forEach((_,id)=>{ if(before.has(id)) already++; else added++; });
  if(merge==='replace') before.forEach(id=>{ if(!matched.has(id)) removed++; });

  const cls = unmatched.length ? 'warn' : 'ok';
  const delta = (typeof Handoff!=='undefined')
    ? Handoff.deltaText({added, already, removed, missing:unmatched.length, total:S.selected.size})
    : `selection is now <b>${S.selected.size}</b> accessions`;
  uplSay(cls, `Loaded <b>${matched.size}</b> of ${considered} entries from <b>${escAttr(source)}</b> · ${delta}. ${undoChip()}`);
  renderUplReport({matched:matched.size, considered, unmatched, dupes, fanout, source, skippedHeader, applied:true});

  renderAccList(); renderSelected(); renderRunbar();
}

function renderUplReport(r){
  const box=document.getElementById('uplReport'); if(!box) return;
  const bits=[];
  if(r.skippedHeader) bits.push(`<div class="upl-note">Skipped header row: <span class="mono">${escAttr(r.skippedHeader)}</span></div>`);
  if(r.dupes.length)  bits.push(`<div class="upl-note">${r.dupes.length} duplicate entr${r.dupes.length>1?'ies were':'y was'} ignored.</div>`);
  if(r.fanout.length) bits.push(`<div class="upl-note">${r.fanout.length} entr${r.fanout.length>1?'ies':'y'} matched more than one run (all replicates were selected).</div>`);
  if(r.unmatched.length){
    const show=r.unmatched.slice(0,25);
    bits.push(`<details class="upl-bad" open><summary>${r.unmatched.length} entr${r.unmatched.length>1?'ies':'y'} not found in this dataset</summary>
      <ul>${show.map(u=>`<li><span class="ln">line ${u.line}</span> <span class="mono">${escAttr(u.text)}</span></li>`).join('')}</ul>
      ${r.unmatched.length>show.length?`<div class="upl-note">…and ${r.unmatched.length-show.length} more.</div>`:''}
      <div class="upl-note">Accepted forms: SNPVersity ID, run accession (SRR/ERR/DRR/CRR…), founder or line name, GRIN accession (PI…), or <span class="mono">FOUNDER_RUN</span>.</div>
    </details>`);
  }
  box.innerHTML=bits.join('');
}

/* minimal styling for the uploader (kept local so it can't clash with the suite CSS) */
function injectVersityCSS(){
  if(document.getElementById('snpversity-upl-css')) return;
  const s=document.createElement('style'); s.id='snpversity-upl-css';
  s.textContent=`
    .upl-name{font-size:11.5px;color:var(--muted);margin:6px 0 6px;word-break:break-all}
    .upl-paste{width:100%;box-sizing:border-box;min-height:64px;resize:vertical;border:1px solid var(--line);
      border-radius:8px;padding:7px 9px;font-family:var(--mono);font-size:11.5px;color:var(--ink);background:#fff}
    .upl-actions{display:flex;gap:6px;flex-wrap:wrap;margin-top:7px}
    .upl-opt{margin-top:9px}
    .upl-status{font-size:11.5px;margin-top:8px;line-height:1.45}
    .upl-status.ok{color:#176c3a} .upl-status.err{color:#c0362c} .upl-status.warn{color:#8a6d1e}
    .upl-report{margin-top:8px}
    .upl-note{font-size:11px;color:var(--muted);margin-top:5px}
    .upl-bad{margin-top:7px;border:1px solid #f0c4bd;background:#fdf6f5;border-radius:8px;padding:7px 9px}
    .upl-bad summary{cursor:pointer;font-size:11.5px;font-weight:600;color:#8f281c}
    .upl-bad ul{margin:7px 0 0;padding-left:16px;max-height:190px;overflow:auto}
    .upl-bad li{font-size:11px;margin-bottom:3px}
    .upl-bad .ln{color:var(--faint);margin-right:5px}
    .upl-bad .mono,.upl-note .mono,.upl-name .mono{font-family:var(--mono)}
    .from-fn{background:#eef4ff;border:1px solid #cfe0ff;color:#274b8f;border-radius:9px;padding:9px 12px;
      font-size:12.5px;margin-bottom:12px;display:flex;gap:9px;align-items:center;flex-wrap:wrap}
    .from-fn .mono{font-family:var(--mono)}
    .from-fn-acts{margin-left:auto;display:flex;gap:8px;align-items:center;flex-wrap:wrap}
    .from-fn .ho-hint{color:#4a6ca8}
    /* step 3 (selection summary + quick picks; the full selector is SNPTrait's drawer) */
    .acc2{padding:16px 20px;display:flex;flex-direction:column;gap:13px;min-width:0}
    .acc2-top{display:flex;align-items:center;gap:14px;flex-wrap:wrap}
    .acc2-count{flex:1;min-width:180px}
    .acc2-n{font-family:var(--mono);font-size:24px;font-weight:700;color:var(--ink)}
    .acc2-of{font-size:13px;color:var(--muted)}
    .acc2-bar{display:flex;height:10px;border-radius:6px;overflow:hidden;background:#eef1f5}
    .acc2-bar span{display:block;min-width:3px}
    .acc2-legend{display:flex;flex-wrap:wrap;gap:4px 14px;margin-top:7px;font-size:12px;color:var(--muted)}
    .acc2-legend i{display:inline-block;width:9px;height:9px;border-radius:2px;margin-right:5px;vertical-align:-1px}
    .acc2-legend b,.acc2-sub b{color:var(--ink);font-family:var(--mono);font-weight:600}
    .acc2-sub{font-size:12px;color:var(--muted);margin-top:4px}
    .acc2-chips{max-height:150px;padding:0}
    .acc2-chips .empty{text-align:left;padding:6px 0}
    .link-more{border:0;background:none;color:var(--blue-600);font:600 12px var(--body);cursor:pointer;padding:4px 6px}
    .acc2-picks{border-top:1px solid var(--line);padding-top:12px}
    .acc2-k{font-size:11px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.4px;margin-bottom:8px}
    .acc2-k span{text-transform:none;letter-spacing:0;font-weight:400;color:var(--faint);margin-left:6px}
    .acc2-pickrow{display:flex;flex-wrap:wrap;gap:6px}
    .acc2-pickrow .qn{font-family:var(--mono);color:var(--faint);margin-left:3px;font-weight:500}
    .acc2-pickrow .qbtn.on{background:var(--blue-50);border-color:#a9c5fb;color:var(--blue-600)}
    .acc2-pickrow .qbtn:disabled{opacity:.45;cursor:not-allowed}
    .acc2-say{font-size:12px;color:var(--muted);margin-top:8px;min-height:16px}
    .gene-more{display:inline-block;margin-left:4px;padding:0 5px;border-radius:5px;background:#eef2f8;
      color:var(--muted);font-size:10.5px;font-weight:600;cursor:help;vertical-align:1px}
    /* per-row jump links in the Effect column (matches .pe-jump) */
    .effect-cell .fold-jump{font-size:11px;white-space:nowrap;margin-left:2px;
      color:#176c3a;text-decoration:none;border-bottom:1px dotted #9ecdb1}
    .effect-cell .fold-jump:hover{color:#0f4f2a;border-bottom-style:solid}
    table.vcf thead th.annot-off{opacity:.62;font-style:italic}
    table.vcf tbody td.annot-off{background:repeating-linear-gradient(135deg,#f6f8fb 0 6px,#eef1f6 6px 12px)}
    .annot-note{margin:6px 0 10px;padding:8px 11px;border:1px solid var(--line);border-left:3px solid #c9a227;
      border-radius:6px;background:#fffdf4;font-size:12px;color:var(--ink);line-height:1.5}
    .annot-note .annot-k{font-weight:600}
    .annot-note summary{cursor:pointer;color:var(--muted);margin-top:3px}
    .annot-note .annot-why{margin-top:4px;color:var(--muted)}`;
  document.head.appendChild(s);
}

/* every one of these is a hand edit, so it expires the pending Undo */
function toggleAcc(id){touchSelection();accSay('');S.selected.has(id)?S.selected.delete(id):S.selected.add(id);renderAccList();renderSelected();renderRunbar();}

/* ---- run bar ---- */
function renderRunbar(){
  const rb=document.getElementById('runbar'); if(!rb)return;
  const d=DATASETS.find(x=>x.id===S.dataset);
  const lo=Math.min(S.start,S.end),hi=Math.max(S.start,S.end);
  const span=hi-lo;
  const ready=S.selected.size>0 && span>0;
  // rough pre-run prediction (the real table-vs-download call uses the exact
  // variant count the server returns).
  const est=ready ? Data.estimateResult(S.dataset, span, S.selected.size) : null;
  const mode = est ? (est.overLimit
    ? '<span style="color:#b42318">likely over the server’s build limit — choose a smaller interval or fewer accessions</span>'
    : est.willDownload
    ? '<span style="color:var(--faint)">likely too large for the table — returns a downloadable VCF</span>'
    : 'opens as an interactive table') : '';
  const spanStr = span>=1e6 ? (span/1e6).toFixed(2)+' Mb' : (span/1e3).toFixed(1)+' kb';
  rb.innerHTML=`
    <div class="summ">
      <span class="kick${ready?'':' wait'}">${ready?'Ready to run':'Finish steps 1–3 to run'}</span>
      <b>${d.name} ${d.sub}</b> · <span class="mono">${S.chr}:${lo.toLocaleString()}–${hi.toLocaleString()}</span><br>
      <span>${S.selected.size} accessions · ${spanStr} interval${mode?' · '+mode:''}</span>
    </div>
    <div class="spacer"></div>
    <button class="btn primary" ${ready?'':'disabled'} onclick="runQuery()">
      ${ICONS.dna} Build VCF & view
    </button>`;
}

/* ---- run query (LIVE: processForm.php -> h5_to_vcf.py -> VCF) ---- */
function noticeCard(title, body){
  return `<div class="card pad fade" style="text-align:center">
    <h3 style="font-family:var(--disp);margin:0 0 8px">${title}</h3>
    <p style="color:var(--muted);margin:0 auto;max-width:560px">${body}</p></div>`;
}
/* The query a result answers. The region form and the dataset (picked here, or in SNPTrait,
   SNPFold or SNPFunction) can move on after a run; the table, its links and every hand-off
   follow what was queried, not what the form says now. A result built outside runQuery
   (the test harness) has no record and reads the current form. */
function resultQuery(res){
  res=res||S.results||{};
  return res.q || {dataset:S.dataset, chr:S.chr, lo:Math.min(S.start,S.end), hi:Math.max(S.start,S.end)};
}
/* A result as the other tools take it (SNPTree, SNPMatrix, SNPCompare, SNPImpact, SNPGeo). */
function resultHandoff(res){
  res=res||S.results;
  const q=resultQuery(res);
  return {rows:res.rows, accs:res.accs, chr:q.chr, start:q.lo, end:q.hi, dataset:q.dataset,
    datasetName:(Data.datasets().find(d=>d.id===q.dataset)||{}).name||q.dataset, vcfUrl:res.vcfUrl};
}
let RUN_SEQ=0;               // a second Build overtakes the first; the first's answer is dropped
async function runQuery(){
  touchSelection();          // the selection is committed — the Undo no longer applies
  const anchor=document.getElementById('resultsAnchor');
  anchor.innerHTML=`<div class="card" style="overflow:hidden"><div class="loading">
    <div class="spinner"></div>
    <div>Reading HDF5 store and assembling VCF…</div>
    <div class="load-hint">
      Large regions and big accession sets can take a while. A smaller region, fewer accessions,
      or a lower-density SNP set will load faster.
    </div>
  </div></div>`;
  anchor.scrollIntoView({behavior:'smooth',block:'start'});
  const lo=Math.min(S.start,S.end), hi=Math.max(S.start,S.end);
  const q={dataset:S.dataset, chr:S.chr, lo, hi}, seq=++RUN_SEQ;
  try{
    const res = await Data.queryVariants(q.dataset, q.chr, lo, hi, [...S.selected]);
    if(seq!==RUN_SEQ) return;
    res.q=q; S.results=res;
    if(S.results.wide){
      anchor.innerHTML=wideResultCard(S.results, lo, hi);
      return;
    }
    if(S.results.empty || !S.results.rows.length){
      anchor.innerHTML=noticeCard('No variants in this range',
        'The query returned no variant sites for the selected accessions. Try a wider interval or a different region.');
      return;
    }
    S.page=1; renderResults();
  }catch(err){
    if(seq!==RUN_SEQ) return;
    const url=err&&err.detail&&err.detail.vcfUrl;   // set when the VCF was built but the fetch failed
    anchor.innerHTML=noticeCard('The query failed',
      `<span style="color:var(--muted)">${(err&&err.message?err.message:err)}</span>`+
      (url?`<br><br><button class="btn" onclick="downloadVcfUrl('${url}')">${ICONS.download} Download the VCF</button>`:''));
    console.error('SNPVersity query failed:', err);
  }
}
/* Result is too large for the in-browser table — VCF download only. */
function wideResultCard(res, lo, hi){
  const V=res.variants, A=(res.accs||[]).length;
  const size = V!=null ? `<b>${V.toLocaleString()}</b> variants x <b>${A}</b> accessions`
                       : `a <b>${((hi-lo)/1e6).toFixed(2)} Mb</b> interval`;
  return `<div class="card pad fade" style="text-align:center">
    <h3 style="font-family:var(--disp);margin:0 0 8px">Result too large for the table view</h3>
    <p style="color:var(--muted);margin:0 auto 14px;max-width:560px">This query returns ${size}.
    Download the VCF to work with it directly, or else choose a smaller region, fewer accessions, or lower-density SNP set.</p>
    <button class="btn primary" onclick="downloadVCF()">${ICONS.download} Download the VCF</button>
  </div>`;
}
/* download a VCF by explicit URL */
function downloadVcfUrl(url){
  if(!url)return;
  const a=document.createElement('a');
  a.href=url; a.download=url.split('/').pop();
  document.body.appendChild(a); a.click(); a.remove();
}
/* download the VCF produced by the most recent query */
function downloadVCF(){ downloadVcfUrl(S.results&&S.results.vcfUrl); }
function gColor(s){
  if(s===null||s===undefined) return '#e7ebf1';
  const min=-15,max=10, c=Math.max(min,Math.min(s,max));
  if(c<=1&&c>=-1) return '#aab4be';
  const r=(c-min)/(max-min);
  return `rgb(${Math.round(255*(1-r))},${Math.round(255*r)},0)`;
}

/* Hard "you can't send this there" gate for the result toolbar — greys a
   button only when the target tool would almost certainly kill the tab
   (memory ceiling, or a ~minute-long synchronous compute). Slower-but-doable
   sizes stay enabled and are caught by each tool's own "build anyway" notice.
   The header "Send selection to…" menu isn't greyed; the tool guards catch it. */
function sendGate(tool){
  const V = (S.results && S.results.rows) ? S.results.rows.length : 0;
  const A = (S.results && S.results.accs) ? S.results.accs.length : 0;   // the matrix the tool receives
  if(tool==='impact')
    return V > IBS_COST.impactBlock
      ? {ok:false, msg:'Selected data too large for SNPImpact - choose a smaller region or lower-density SNP set'}
      : {ok:true};
  if(tool==='tree')
    return ibsWork(V,A) > IBS_COST.workBlock
      ? {ok:false, msg:'Selected data too large for SNPTree - choose a smaller region, a lower-density SNP set, or fewer accessions'}
      : {ok:true};
  if(tool==='compare')
    return (ibsWork(V,A) > IBS_COST.cmpWorkBlock || V*A > IBS_COST.cmpMemBlock)
      ? {ok:false, msg:'Selected data too large for SNPCompare - choose a smaller region, a lower-density SNP set, or fewer accessions'}
      : {ok:true};
  if(tool==='geo')
    return !(typeof SNPTools!=='undefined' && SNPTools.registry && SNPTools.registry.snpgeo)
      ? {ok:false, msg:'SNPGeo is not loaded on this deployment'}
      : (V > 400000 ? {ok:false, msg:'Too many variants for SNPGeo - choose a smaller region'} : {ok:true});
  return {ok:true};
}
function sendBtnHTML(fn, gate, label, icon){
  return `<button class="btn${gate.ok?'':' off'}" ${gate.ok?'':`aria-disabled="true" data-tt="${escAttr(gate.msg)}"`} onclick="${fn}">${icon||''} ${label}</button>`;
}

function renderResults(){
  const a=document.getElementById('resultsAnchor');
  const {chr,lo,hi}=resultQuery();
  const pgOK=(hi-lo)<=PANGENOME_MAX_SPAN;
  a.innerHTML=`
    <div class="sec"><div class="bar"></div>
      <div style="display:flex;align-items:center;gap:14px;flex-wrap:wrap;width:100%">
        <h2 style="font-size:16px;margin:0">Variant view · <span class="c-mono" style="color:var(--blue-600)">${chr}:${lo.toLocaleString()}–${hi.toLocaleString()}</span></h2>
        <div style="margin-left:auto;display:flex;gap:8px;flex-wrap:wrap">
          ${sendBtnHTML('sendToImpact()',  sendGate('impact'),  'Send to SNPImpact',   ICONS.star)}
          ${sendBtnHTML('sendToCompare()', sendGate('compare'), 'Send to SNPCompare',  ICONS.compare||ICONS.grid)}
          ${sendBtnHTML('sendToTree()',    sendGate('tree'),    'Send data to SNPTree', ICONS.tree)}
          ${sendBtnHTML('sendToGeo()',     sendGate('geo'),     'Send to SNPGeo',      ICONS.map||'')}
          <button class="btn${pgOK?'':' off'}" ${pgOK?'':'aria-disabled="true" data-tt="Selected region too large for Pangenome viewer - choose a region under 10 kb"'} onclick="openPangenomeRegion()">${ICONS.grid||''} Pangenome viewer ↗</button>
          <button class="btn" onclick="downloadVCF()">${ICONS.download} Download VCF</button>
        </div>
      </div>
    </div>
    <div class="result-tabs">
      <button class="rtab active" data-rt="table" onclick="switchRT('table')">${ICONS.table} Table view</button>
    </div>
    <div id="rtBody"></div>`;
  switchRT('table');
}
function switchRT(rt){
  document.querySelectorAll('.rtab').forEach(t=>t.classList.toggle('active',t.dataset.rt===rt));
  return renderTable();
}
/* hand the generated VCF matrix (+ metadata) to SNPTree for a local IBS phylogeny */
/* "Load data from SNPVersity": SNPTree, SNPImpact and SNPCompare open empty until a result is sent
   to them; when SNPVersity holds a table result, their empty page offers it in one click (the
   button runs the same hand-off as SNPVersity's own Send button). Nothing when there is no result
   or it was a VCF-only (wide) one. */
function versityResultReady(){ return !!(S.results && S.results.rows && S.results.rows.length && !S.results.wide && S.results.accs); }
function loadFromVersityHTML(sendFn, what){
  if(!versityResultReady()) return '';
  const q=resultQuery(), r=S.results;
  return `<div class="load-versity" id="loadVersity">
    <div><b>SNPVersity has a result</b> <span class="lv-sub">${escAttr(q.chr)}:${(+q.lo).toLocaleString()}–${(+q.hi).toLocaleString()} ·
      ${r.accs.length.toLocaleString()} lines · ${r.rows.length.toLocaleString()} variants</span>
      <div class="lv-sub">${what}</div></div>
    <button class="btn primary" onclick="${sendFn}()">${ICONS.dna||''} Load data from SNPVersity</button></div>`;
}
function sendToTree(){
  if(S.results&&S.results.rows&&S.results.rows.length) S.treeInput=resultHandoff();
  go('snptree');   // navigates even with no result yet (SNPTree shows a guided empty state)
}
/* hand the region's genotype matrix to SNPMatrix (IBS distance heatmap) */
function sendToMatrix(){
  if(S.results&&S.results.rows&&S.results.rows.length) S.matrixInput=resultHandoff();
  go('snpmatrix');
}
/* hand the generated VCF matrix (+ metadata) to SNPCompare for local vs global IBS */
function sendToCompare(){
  if(S.results&&S.results.rows&&S.results.rows.length) S.compareInput=resultHandoff();
  go('snpcompare');
}
/* hand the region's variants to SNPGeo for the geographic view.
   SNPGeo's carrier percentages are over EVERY sample of the dataset known from a
   country, so a partial selection would understate them. When the current result
   covers only part of the dataset, re-query the same region for all of the
   dataset's samples (as SNPGeo's own gene search does), provided the table
   budget allows it; otherwise hand off what we have and SNPGeo shows a
   partial-selection warning. */
const GEO_HANDOFF_WORK_MAX = 4e7;   // = Data TABLE_WORK_MAX (variants x samples)
async function sendToGeo(){
  if(!(S.results&&S.results.rows&&S.results.rows.length)){ go('snpgeo'); return; }
  const q=resultQuery();
  const allIds=(Data.accessionsFor(q.dataset)||[]).map(a=>a.id);
  let res=S.results, requeried=false;
  const partial=(S.results.accs||[]).length < allIds.length;
  if(partial && S.results.rows.length*allIds.length <= GEO_HANDOFF_WORK_MAX){
    const btns=[...document.querySelectorAll('button')].filter(b=>/Send to SNPGeo/.test(b.textContent));
    btns.forEach(b=>{b.dataset.lbl=b.innerHTML; b.innerHTML='Querying all '+allIds.length+' samples…'; b.disabled=true;});
    try{
      const r=await Data.queryVariants(q.dataset, q.chr, q.lo, q.hi, allIds);
      if(r && !r.wide && r.rows && r.rows.length){ r.q=q; res=r; requeried=true; }
    }catch(e){ console.warn('[sendToGeo] full-dataset re-query failed; handing off the selection only', e); }
    btns.forEach(b=>{b.innerHTML=b.dataset.lbl; b.disabled=false;});
  }
  S.geoInput=Object.assign(resultHandoff(res), {requeried});
  go('snpgeo');
}
/* hand the region's variants to SNPImpact for ranking (accessions irrelevant there) */
function sendToImpact(){
  if(S.results&&S.results.rows&&S.results.rows.length) S.impactInput=resultHandoff();
  go('snpimpact');
}

function setMaf(el){
  let v=parseFloat(el.value);
  if(isNaN(v)) v=0;
  v=Math.max(0, Math.min(0.5, v));      // clamp to valid MAF range
  v=Math.round(v*100)/100;              // snap to 0.01
  S.fMaf=v; S.page=1; renderTable();
}
/* Cache the O(rows) work (full-table filter + distinct-effects scan) so that
   paging — which doesn't change the filter — doesn't re-scan every variant on
   each click. Keyed on the result set + the four filter controls. The Site QC
   counts are taken once per result. */
let _tblCache=null, _qcCache={res:null, val:null};
function versityQc(){ return S.fQc || Data.SITE_QC_DEFAULTS.versity; }
function qcSummary(rows){
  const res=S.results||null;
  if(_qcCache.res===res && _qcCache.val) return _qcCache.val;
  const n={flagged:0, nocarrier:0, any:false};
  for(const r of rows){ if(r.qc==null) continue; n.any=true;
    const g=Data.qcGroup(r.qc); if(g==='flagged') n.flagged++; else if(g==='nocarrier') n.nocarrier++; }
  _qcCache={res, val:n};
  return n;
}
function tableView(rows){
  const fq=versityQc();
  const key=S.results&&S.results.vcfUrl+'|'+S.fImpact+'|'+S.fEffect+'|'+S.fMaf+'|'+fq;
  if(_tblCache&&_tblCache.key===key) return _tblCache;
  const qc=qcSummary(rows);
  const fr=rows.filter(r=>
    (S.fImpact==='all'||r.impact===S.fImpact) &&
    (S.fEffect==='all'||r.effect===S.fEffect) &&
    (r.maf>=S.fMaf) && Data.qcKeep(r.qc, fq));
  const effects=['all',...new Set(rows.map(r=>r.effect))];
  _tblCache={key,fr,effects,qc};
  return _tblCache;
}

function renderTable(){
  const {rows,accs}=S.results;
  // accession header height scales to the longest full ID so it isn't clipped
  const maxIdLen=accs.length?Math.max(...accs.map(a=>String(a.id).length)):8;
  const hetMarks=accs.some(a=>a.sampleQC==='Heterozygous sample');     // room for the "het" marker
  const thH=Math.max(118, Math.min(300, Math.round(maxIdLen*6.4)+30+(hetMarks?28:0)));
  const {fr,effects,qc}=tableView(rows);
  const perPage=S.perPage, pages=Math.max(1,Math.ceil(fr.length/perPage));
  if(S.page>pages)S.page=1;
  const slice=fr.slice((S.page-1)*perPage, S.page*perPage);

  const b=document.getElementById('rtBody');
  b.innerHTML=`
    <div class="filterbar">
      <div class="fld"><label>Effect impact</label>
        <select onchange="S.fImpact=this.value;S.page=1;renderTable()">
          ${['all','HIGH','MODERATE','LOW','MODIFIER'].map(v=>`<option value="${v}" ${v===S.fImpact?'selected':''}>${v==='all'?'All impacts':v}</option>`).join('')}
        </select></div>
      <div class="fld"><label>Effect type</label>
        <select onchange="S.fEffect=this.value;S.page=1;renderTable()">
          ${effects.map(v=>`<option value="${v}" ${v===S.fEffect?'selected':''}>${v==='all'?'All effects':v}</option>`).join('')}
        </select></div>
      <div class="fld"><label>Min. MAF (0.0–0.5)</label>
        <input type="number" min="0" max="0.5" step="0.01" value="${S.fMaf.toFixed(2)}" inputmode="decimal"
          title="Enter a value between 0.0 and 0.5" onchange="setMaf(this)"
          onkeydown="if(event.key==='Enter'){this.blur();}">
      </div>
      ${qc.any?`<div class="fld"><label data-tt="Site QC classes each site by its heterozygous and homozygous carriers among all lines of the release. Flagged: het only or het excess. The downloaded VCF keeps every site.">Site QC</label>
        <select id="versityQc" onchange="S.fQc=this.value;S.page=1;renderTable()">
          ${Data.SITE_QC_CHOICES.map(([v,l])=>`<option value="${v}" ${v===versityQc()?'selected':''}>${l}</option>`).join('')}
        </select></div>`:''}
    </div>
    ${qc.any?`<div class="qc-note" id="versityQcNote">${qc.flagged.toLocaleString()} flagged, ${qc.nocarrier.toLocaleString()} with no carrier in this release, of ${rows.length.toLocaleString()} sites.</div>`:''}
    <div class="legend">
      <span style="font-weight:600;color:var(--ink)">Genotype</span>
      <span class="li"><span class="sw" style="background:#e9f4ec"></span>Reference allele (0)</span>
      <span class="li"><span class="sw" style="background:#fdf3d2"></span>Alternative heterozygous allele (1)</span>
      <span class="li"><span class="sw" style="background:#f4a259"></span>Alternative homozygous allele (2)</span>
      <span class="li"><span class="sw" style="background:#eaeef4"></span>Missing (.)</span>
      <span style="margin-left:14px;font-weight:600;color:var(--ink)">LM score</span>
      <span class="li"><span class="sw" style="width:54px;background:linear-gradient(90deg,rgb(255,0,0),#aab4be,rgb(0,255,0))"></span>deleterious → tolerated</span>
    </div>
    ${annotNoteHTML()}
    ${genesPanel()}
    ${curatedLineHTML()}
    <div class="tbl-wrap"><table class="vcf">
      <thead><tr>
        <th data-tt="Chromosome — reference chromosome containing the variant.">CHR</th><th data-tt="Position — 1-based coordinate on B73 v5.">POS</th><th class="num" data-tt="Reference allele in B73 v5.">REF</th><th class="num" data-tt="Alternate allele represented by this row.">ALT</th>
        ${annotHeaderHTML()}
        ${accs.map(a=>`<th class="acc-th" style="height:${thH}px" title="${escAttr((a.projTitle||a.proj||'')+' — '+a.id)}"><span class="proj-bar" style="background:${a.projColor};height:8px" title="${escAttr(a.projTitle||a.proj||'')}"></span><span class="v">${a.id}${sampleHetMark(a)}</span></th>`).join('')}
      </tr></thead>
      <tbody>
        ${slice.map(r=>rowHTML(r)).join('')}
      </tbody>
    </table></div>
    <div class="pager">
      <div class="info">Showing <span class="mono">${slice.length}</span> of <span class="mono">${fr.length}</span> loci · page <span class="mono">${S.page}/${pages}</span></div>
      <div class="pages">
        <button onclick="setPage(${Math.max(1,S.page-1)})">‹</button>
        ${pageBtns(pages)}
        <button onclick="setPage(${Math.min(pages,S.page+1)})">›</button>
      </div>
    </div>`;
  attachTT();
}
/* Annotation columns, in this order. Per-set availability comes from
   Data.annotationFields(): a column a set lacks ('na') stays in the table with empty
   cells, a dashed header and the reason in its tooltip; a column the set never records
   ('hidden': MQ and COMP for the GRIN-linked 2026 call set) is left out. */
const ANNOT_TT={
  gene:'B73 v5 gene model overlapping the variant; where a site has consequences in several genes, the most severe is shown and +N lists the others.',
  effect:'Predicted consequence — e.g. missense, synonymous, intron, frameshift.',
  impact:'Predicted severity — HIGH, MODERATE, LOW, or MODIFIER.',
  domain:'Pfam protein domain overlapping the affected residue, when available.',
  mq:'Mapping quality — confidence that reads aligned to the correct location; higher is better.',
  comp:'Completeness — fraction of accessions with a non-missing call at this site.',
  r2:'Maximum LD r² — linkage-disequilibrium correlation; closer to 1 is stronger.',
  maf:'Minor-allele frequency — frequency of the less common allele (0 to 0.5).',
  qc:'Site QC — the class of the site from its heterozygous and homozygous carriers among all lines of the release. Pass: fewer than a quarter of carriers heterozygous. Het elevated: a quarter or more. Het excess: more heterozygous than homozygous carriers. Het only: no homozygous carrier. No carrier: no line carries the allele. Hover a cell for its counts.',
  pc1:'PlantCAD DNA language-model score; more extreme values are more disruptive.',
  pc2:'Second-generation PlantCAD DNA score.',
  evo2:'Evo2 DNA language-model score (log-likelihood ratio), for SNPs within 1 kb of a gene; more negative is more disruptive.',
  esm1:'ESM protein language-model score for the amino-acid change.',
  esm2:'ESM2 protein language-model score.',
  esm3:'ESM3 protein language-model score.',
  esmc:'ESM-C protein language-model score for the amino-acid change (missense sites).'};
const ANNOT_NUM={mq:1,comp:1,r2:1,maf:1,pc1:1,pc2:1,evo2:1,esm1:1,esm2:1,esm3:1,esmc:1};
let _annotCache={ds:null, res:null, val:null};
function annotFields(){
  // once per dataset + result: rowHTML asks for every row, and a 'pending' column scans the rows
  const res=S.results||null, ds=resultQuery().dataset;
  if(_annotCache.val && _annotCache.ds===ds && _annotCache.res===res) return _annotCache.val;
  const F=(Data.annotationFields ? Data.annotationFields(ds) : Object.keys(ANNOT_TT).map(k=>({key:k,label:k,status:'ok',note:''})));
  // A 'pending' column is shown as soon as the queried store carries it (e.g. PlantCAD merged
  // for one chromosome first); it stays pending (empty + note) where the store has no value.
  const rows=(res&&res.rows)||[];
  // An 'auto' column (Site QC) is shown only when a row of the result carries a value, and
  // left out otherwise, so a store without a sidecar shows the table as before.
  const val=F.filter(f=>f.status!=='hidden' && (f.status!=='auto' || rows.some(r=>r[f.key]!=null)))
    .map(f=>f.status==='auto' ? Object.assign({}, f, {status:'ok'})
      : (f.status==='pending' && rows.some(r=>r[f.key]!=null))
      ? Object.assign({}, f, {status:'ok', note:'Merged for this chromosome (rounded to 0.1); sites without a score show N/A.'}) : f);
  _annotCache={ds, res, val};
  return val;
}
function annotHeaderHTML(){
  return annotFields().map(f=>{
    const tt=ANNOT_TT[f.key]+(f.status==='pending'?' Pending for this set: '+f.note
      :(f.status==='na'?' Not available for this set: '+f.note:(f.note?' '+f.note:'')));
    const cls=[ANNOT_NUM[f.key]?'num':'', f.status!=='ok'?'annot-off':''].filter(Boolean).join(' ');
    return `<th${cls?` class="${cls}"`:''} data-col="${f.key}" data-status="${f.status}" data-tt="${escAttr(tt)}">${scoreHeadHTML(f.label)}</th>`;
  }).join('');
}
function annotNoteHTML(){
  const F=annotFields(), na=F.filter(f=>f.status==='na'), pend=F.filter(f=>f.status==='pending');
  if(!na.length && !pend.length) return '';
  const dsId=resultQuery().dataset, ds=(Data.datasets().find(d=>d.id===dsId)||{});
  const nm=escAttr((ds.name||dsId)+(ds.sub?' · '+ds.sub:''));
  const li=L=>L.map(f=>`<b>${f.label}</b> — ${escAttr(f.note)}`).join('<br>');
  return `<div class="annot-note" id="annotNote">`+
    (na.length?`<div><span class="annot-k">Not available for this set (${nm}):</span> ${na.map(f=>f.label).join(', ')}. The columns stay in the table with empty cells.</div>`:'')+
    (pend.length?`<div><span class="annot-k">Pending for this set:</span> ${pend.map(f=>f.label).join(', ')} — empty until the merge lands.</div>`:'')+
    `<details><summary>Why</summary><div class="annot-why">${li(na.concat(pend))}</div></details></div>`;
}
/* SNPCurate: the mark of a curated allele beside its position, and a line above the table naming
   the curated alleles in the queried interval, grey ones (not genotyped here) included. */
function curatedLineHTML(){
  const q=resultQuery(), L=Data.curatedInInterval ? Data.curatedInInterval(q.chr, q.lo, q.hi) : [];
  if(!L.length) return '';
  return `<div class="qc-note cur-line" id="versityCurated">Curated alleles in this interval: ${L.map(e=>
    `<span class="cur-item">${curateBadge(e)} ${escAttr(e.symbol)} ${escAttr(e.label)}</span>`).join(' · ')}</div>`;
}
function rowHTML(r){
  const chr=resultQuery().chr;               // the queried chromosome, even if the form has moved on
  const lo=Math.max(1,r.pos-10000),hi=r.pos+10000;
  const link=`https://jbrowse.maizegdb.org/?data=B73&loc=${chr}:${lo}..${hi}&highlight=${chr}:${r.pos}..${r.pos}`;
  const effectText = String(r.effect||'');
  const isMis = /missense/i.test(effectText);
  const isSyn = /synonymous/i.test(effectText) && !/non[-_ ]?synonymous/i.test(effectText);
  const isFoldCoding = /missense|protein[_ -]?altering|non[_ -]?synonymous|stop[_ -]?gained|nonsense|frameshift|start[_ -]?lost|initiator[_ -]?codon|stop[_ -]?lost|inframe[_ -]?(insertion|deletion)/i.test(effectText);
  const hasGene = r.gene && r.gene!=='—';
  const peJump = (isMis && hasGene)
    ? ` <a class="pe-jump" href="#" title="View this substitution in PanEffect" onclick="goPanEffect('${r.gene}',{variant:'${escAttr(r.sub||'')}'});return false;">effects ↗</a>`
    : '';
  /* SNPFold: residue-level coding changes and severe LOF changes get a structure
     link. The genomic position/alleles are also passed so stop-gained records can
     still be selected when the VCF SUB field lacks a protein substitution. */
  const foldJump = (isFoldCoding && hasGene)
    ? ` <a class="fold-jump" href="#" title="Show this variant on the predicted protein structure in SNPFold" onclick="goFold('${r.gene}',{chr:'${escAttr(chr)}',pos:${r.pos},ref:'${escAttr(r.ref)}',alt:'${escAttr(r.alt)}',sub:'${escAttr(r.sub||'')}',effect:'${escAttr(r.effect||'')}'});return false;">fold ↗</a>`
    : '';
  const eff = (r.sub?`<span class="sub">(${r.sub})</span> ${r.effect}`:r.effect) + peJump + foldJump;
  const fields=annotFields();
  const F={}; fields.forEach(f=>{F[f.key]=f;});
  const cur=Data.curatedAt ? Data.curatedAt(chr, r.pos, r.ref, r.alt) : null;   // SNPCurate
  const off=k=>F[k] && F[k].status!=='ok';           // column not available / pending for this set
  const blank=k=>`<td class="num annot-off" data-tt="${escAttr(F[k].note)}"></td>`;
  const sc=(v,k)=>off(k)?blank(k):`<td class="score ${v===null?'na':''}" style="${v===null?'':'background:'+gColor(v)}">${v===null?'N/A':v}</td>`;
  // one cell per annotation column, in the header's order (a hidden column has neither)
  const cell={
    gene:  ()=>`<td>${geneCell(r.gene)}${moreGenes(r)}</td>`,
    effect:()=>`<td class="effect-cell">${eff}</td>`,
    impact:()=>`<td><span class="pill ${r.impact.toLowerCase()}">${r.impact}</span></td>`,
    domain:()=>`<td>${domTag(r.domain)}</td>`,
    mq:    ()=>off('mq')?blank('mq'):`<td class="num">${r.mq}</td>`,
    comp:  ()=>off('comp')?blank('comp'):`<td class="num">${r.comp}</td>`,
    r2:    ()=>off('r2')?blank('r2'):`<td class="num">${r.r2==null?'<span style="color:var(--faint)">NA</span>':(+r.r2).toFixed(2)}</td>`,
    maf:   ()=>off('maf')?blank('maf'):`<td class="num">${r.maf==null?'<span style="color:var(--faint)">—</span>':r.maf}</td>`,
    qc:    ()=>`<td class="qc-cell">${siteQcPill(r.qc, r.nHet, r.nHom)}</td>`,
    pc1:()=>sc(r.pc1,'pc1'), pc2:()=>sc(r.pc2,'pc2'), evo2:()=>sc(r.evo2,'evo2'),
    esm1:()=>sc(r.esm1,'esm1'), esm2:()=>sc(r.esm2,'esm2'), esm3:()=>sc(r.esm3,'esm3'), esmc:()=>sc(r.esmc,'esmc'),
  };
  /* a curated allele's row: its site and annotation cells (not the score or genotype cells) on light
     gold, with a gold bar on the left (curated-row CSS in main.css) */
  return `<tr${cur?` class="cur-row cur-row-${cur.mark}" data-curate="${escAttr(cur.id)}"`:''}>
    <td class="c-mono" style="padding-left:11px">${chr.replace('chr','')}</td>
    <td class="c-pos"><a class="gene-link" href="${link}" target="_blank" rel="noopener">${r.pos.toLocaleString()}</a>${cur?curateBadge(cur):''}</td>
    <td class="c-allele c-ref" data-tt="${escAttr(alleleTT(r.ref,'REF allele'))}">${escAttr(alleleDisp(r.ref))}</td>
    <td class="c-allele c-alt" data-tt="${escAttr(alleleTT(r.alt,'ALT allele'))}">${escAttr(alleleDisp(r.alt))}</td>
    ${fields.map(f=>cell[f.key]?cell[f.key]():'<td></td>').join('')}
    ${gtCells(r.gts)}
  </tr>`;
}
/* r.gts is an Int8Array of codes (0 ref, 1 het, 2 alt, 3 missing). Build the
   genotype cells with an indexed loop — Int8Array.map would coerce back to
   numbers, and this is the hottest string-building path in the table. */
const GT_CELL_CLASS=['gt-00','gt-01','gt-11','gt-na'];
const GT_CELL_TEXT =['0','1','2','·'];
const GT_CELL_TT   =['0/0','0/1','1/1','./.'];
function gtCells(gts){
  let s='';
  for(let k=0;k<gts.length;k++){ const c=gts[k];
    s+=`<td class="gt ${GT_CELL_CLASS[c]}" data-tt="${GT_CELL_TT[c]}">${GT_CELL_TEXT[c]}</td>`; }
  return s;
}
/* "Gene model" table cell: link real gene models to their MaizeGDB page.
   maizegdbGeneURL / isSingleGeneModel are shared helpers defined in core.js.
   Intergenic spans and the "—" placeholder stay as plain text. */
/* A site SnpEff annotates in more than one gene (1,464 of 3,744 GRIN-linked test-window sites) shows its
   first, most severe consequence in the Gene model / Effect columns; "+N" names the others
   on hover. SNPFunction, SNPFold and SNPGeo read each gene's own entry (Data.rowForGene). */
function moreGenes(r){
  if(!r.anns || !Data.annotationsOf) return '';
  const seen=new Set(), lines=[];
  Data.annotationsOf(r).forEach((a,i)=>{
    const k=a.gene+'|'+a.effect+'|'+a.sub;
    if(seen.has(k)) return; seen.add(k);
    if(i>0) lines.push(`${a.gene.replace(/\s+/g,'\u2013')}: ${a.effect}${a.sub?` (${a.sub})`:''}`);   // intergenic A_B -> A\u2013B
  });
  if(!lines.length) return '';
  return ` <span class="gene-more" data-tt="${escAttr('Also annotated at this site: '+lines.join('; '))}">+${lines.length}</span>`;
}
function geneCell(g){
  if(!g || g==='—') return '<span style="color:var(--faint)">—</span>';
  return isSingleGeneModel(g)
    ? `<a class="gene-link" href="${maizegdbGeneURL(g)}" target="_blank" rel="noopener" title="MaizeGDB gene page for ${escAttr(g)}">${g}</a>`
    : g;
}
function genesPanel(){
  const rows=(S.results&&S.results.rows)||[];
  // Only single, real gene models — intergenic / boundary spans dropped (isSingleGeneModel, core.js).
  // Every gene a site has a consequence in, not only the first-listed one (Data.genesOf).
  const genes=[...new Set(rows.flatMap(r=>Data.genesOf?Data.genesOf(r):[r.gene]).filter(isSingleGeneModel))].sort();
  if(!genes.length) return '';
  const jb=g=>`https://jbrowse.maizegdb.org/index.html?data=B73&loc=${encodeURIComponent(g)}`;
  const items=genes.map(g=>`<div style="display:flex;flex-wrap:wrap;justify-content:space-between;align-items:center;gap:4px 12px;padding:5px 2px;border-bottom:1px solid #eef1f5">
      <a class="c-mono" style="font-size:12.5px" href="${maizegdbGeneURL(g)}" target="_blank" rel="noopener" title="MaizeGDB gene model page for ${escAttr(g)}">${g}</a>
      <span style="display:flex;flex-wrap:wrap;gap:4px 12px;font-size:12px;white-space:nowrap">
        <a href="${maizegdbGeneURL(g)}" target="_blank" rel="noopener" title="MaizeGDB gene model page">MaizeGDB ↗</a>
        <a href="${jb(g)}" target="_blank" rel="noopener">JBrowse ↗</a>
        <a href="${pangenomeGeneURL(g)}" target="_blank" rel="noopener">Pangenome ↗</a>
        <a href="#" onclick="goPanEffect('${g}');return false;">PanEffect →</a>
        <a href="#" onclick="goFold('${g}');return false;">SNPFold →</a>
        <a href="#" onclick="goFunction('${g}');return false;">SNPFunction →</a>
      </span>
    </div>`).join('');
  return `<details class="card pad" style="margin-bottom:14px">
    <summary style="cursor:pointer;font-weight:600">${genes.length} gene model${genes.length>1?'s':''} in this region</summary>
    <div style="display:grid;grid-template-columns:1fr;gap:0;margin-top:10px">${items}</div>
  </details>`;
}
/* ---- MaizeGDB Pangenome viewer links (B73 v5 coordinates / gene models) ---- */
const PANGENOME_BASE = 'https://pangenome-viewer.maizegdb.org/';
/* The Pangenome viewer only accepts interval spans up to 10 kb; a wider
   region stalls it, so the region hand-off is greyed out past this. */
const PANGENOME_MAX_SPAN = 10000;
/* link for a single B73 gene model */
function pangenomeGeneURL(gene, set){
  return `${PANGENOME_BASE}?set=${encodeURIComponent(set||'NAM')}&geneID=${encodeURIComponent(gene)}`;
}
/* link for a genomic interval */
function pangenomeRegionURL(chr, start, end, set){
  const lo=Math.min(start,end), hi=Math.max(start,end);
  const c=String(chr).startsWith('chr')?chr:`chr${chr}`;
  return `${PANGENOME_BASE}?set=${encodeURIComponent(set||'NAM')}&chr=${encodeURIComponent(c)}&start=${lo}&end=${hi}`;
}
/* open the current query region in the Pangenome viewer */
function openPangenomeRegion(){
  const {chr,lo,hi}=resultQuery();          // the button sits on the result, so its region
  if(hi-lo > PANGENOME_MAX_SPAN) return;   // button is greyed out in this case
  window.open(pangenomeRegionURL(chr,lo,hi),'_blank','noopener');
}
/* Jump to SNPFold for a specific gene model.
   `variant` is optional: {chr,pos,ref,alt,sub,effect} (or a plain "A123T" string).
   When present SNPFold selects that site — lollipop, table row, 3D label — after
   the gene is shown. Passing it through state (not a reload) means clicking fold
   links repeatedly on the same gene never re-queries the HDF5 store. */
function goFold(gene, variant){
  S.foldGene=gene;
  S.foldVariant = variant
    ? (typeof variant==='string' ? {sub:variant} : variant)
    : null;
  go('snpfold');
}
/* jump to SNPFunction (gene function & allele mining) for a specific gene */
function goFunction(gene, dataset){ S.functionGene=gene; S.functionDataset=dataset||S.dataset; go('snpfunction'); }
function pageBtns(pages){
  let out='';const cur=S.page;
  const show=new Set([1,pages,cur,cur-1,cur+1,cur-2,cur+2]);
  let last=0;
  for(let i=1;i<=pages;i++){
    if(!show.has(i))continue;
    if(i-last>1)out+='<span style="align-self:center;color:var(--faint);padding:0 2px">…</span>';
    out+=`<button class="${i===cur?'on':''}" onclick="setPage(${i})">${i}</button>`;
    last=i;
  }
  return out;
}
function setPage(p){S.page=p;renderTable();document.querySelector('.tbl-wrap').scrollTop=0;}


/* register with the suite shell */
SNPTools.register('snpversity', { render(){ renderVersity(); } });
