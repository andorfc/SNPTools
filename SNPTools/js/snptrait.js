/* =====================================================================
 *  snptrait.js — SNPTrait = metadata selector ("Line Selector" on maize).
 *  Registers 'snptrait'. Lets the user filter and select samples by
 *  metadata, then hand the selected set to SNPVersity to build a VCF.
 *
 *  Maize (main): the schema for a family comes from
 *  window.SNPTRAIT_SCHEMA[family], written by build_maize_trait_catalog.py
 *  into js/zmgrin.catalog.js (family zmgrin2026: panel, subpopulation,
 *  GRIN origin, kernel type, binned GRIN evaluation traits). Families
 *  without a generated or built-in schema get a neutral project-only
 *  schema. The Fusarium schemas below are kept so the module stays
 *  identical in behaviour on the fusarium branch.
 *
 *  SCHEMA-DRIVEN: each reference has its own metadata shape, so the facets,
 *  grouping and columns are configured per family (TRAIT_SCHEMA):
 *    - F. graminearum : population structure (NA1/NA2/…), species, host,
 *      chemotype, country. Grouped by population.
 *    - F. verticillioides (7600 + MRC826, same isolates): geography +
 *      substrate — country, substrate group, host species, region. Grouped
 *      by country. (from Fvert_MetaData xlsx, keyed by SNP_ID = HDF5 id.)
 *
 *  Data comes from Data.accessionsFor(dataset). All render fns are global so
 *  the inline handlers resolve.
 * ===================================================================== */

const TRAIT_SCHEMA = {
  graminearum: {
    groupBy:'population', groupLabel:'population structure',
    groupOrder:['NA1','NA2','NA3','Admixture','Outgroups','Unknown'],
    groupColors:{NA1:'#56B4E9',NA2:'#e0c200',NA3:'#009E73',Admixture:'#E69F00',Outgroups:'#CC79A7',Unknown:'#999999'},
    groupNames:{NA1:'North American 1 (NA1)',NA2:'North American 2 (NA2)',NA3:'North American 3 (NA3)',
                Admixture:'Admixture',Outgroups:'Outgroups',Unknown:'Unknown'},
    facets:[['population','Population'],['species','Species'],['host','Host'],['chemotype','Chemotype'],['country','Country']],
    columns:[['species','Species'],['host','Host'],['chemotype','Chemotype'],['country','Country']],
    search:['id','population','species','host','chemotype','country'],
  },
  vert: {   // shared by F. verticillioides 7600 + MRC826 (same isolate set)
    groupBy:'country', groupLabel:'country of origin',
    groupOrder:['USA','CAN','LATAM','ETH','NPL','ZAF','ZMB','FIN','HUN','GER','Unknown'],
    groupColors:{USA:'#2563eb',CAN:'#0e7490',LATAM:'#e0952b',ETH:'#1f8a4c',NPL:'#7c3aed',
                 ZAF:'#cc79a7',ZMB:'#b45309',FIN:'#0891b2',HUN:'#9333ea',GER:'#dc2626',Unknown:'#999999'},
    groupNames:{USA:'United States (USA)',CAN:'Canada (CAN)',LATAM:'Latin America (LATAM)',ETH:'Ethiopia (ETH)',
                NPL:'Nepal (NPL)',ZAF:'South Africa (ZAF)',ZMB:'Zambia (ZMB)',FIN:'Finland (FIN)',
                HUN:'Hungary (HUN)',GER:'Germany (GER)',Unknown:'Unknown origin'},
    facets:[['country','Country'],['substrate','Substrate group'],['host','Host species'],['region','Region']],
    columns:[['strain','Strain'],['substrate','Substrate'],['host','Host species'],['region','Region']],
    search:['id','strain','country','region','substrate','substrateDetail','host','geo'],
  },
};
TRAIT_SCHEMA.vert7600 = TRAIT_SCHEMA.vertMRC826 = TRAIT_SCHEMA.vert;   // fusarium family ids

/* Schema lookup: a generated schema (window.SNPTRAIT_SCHEMA[family]) first,
   then a built-in TRAIT_SCHEMA entry, then a neutral schema — never silently
   another organism's schema. */
const TRAIT_NEUTRAL_SCHEMA = {groupBy:'projTitle', groupLabel:'project', groupOrder:[], groupColors:{}, groupNames:{},
  facets:[['projTitle','Project']], columns:[], search:['id','label','projTitle']};
function traitFamily(dataset){
  return (typeof Data!=='undefined' && Data.familyOf) ? Data.familyOf(dataset) : null;
}
function traitSchema(dataset){
  const fam = traitFamily(dataset);
  const gen = (typeof window!=='undefined' && window.SNPTRAIT_SCHEMA) || {};
  if (fam && gen[fam]) return gen[fam];
  if (fam && TRAIT_SCHEMA[fam]) return TRAIT_SCHEMA[fam];
  return TRAIT_NEUTRAL_SCHEMA;
}
function traitHasSchema(dataset){
  const fam = traitFamily(dataset);
  return !!(fam && (((typeof window!=='undefined' && window.SNPTRAIT_SCHEMA) || {})[fam] || TRAIT_SCHEMA[fam]));
}

/* Organism wording. Maize is the default on this branch; Fusarium families keep
   their isolate wording. A deployment can override per family through
   window.SNPTRAIT_ORGANISM[family] = {kicker, title, one, many, exportBase}. */
const TRAIT_ORGANISM = {
  _default:    {kicker:'LINE SELECTOR',   title:'Select maize lines by metadata',       one:'line',    many:'lines',    exportBase:'selected_lines'},
  graminearum: {kicker:'STRAIN SELECTOR', title:'Select Fusarium isolates by metadata', one:'isolate', many:'isolates', exportBase:'selected_strains'},
};
TRAIT_ORGANISM.vert7600 = TRAIT_ORGANISM.vertMRC826 = TRAIT_ORGANISM.graminearum;
function traitOrganism(dataset){
  const fam = traitFamily(dataset);
  const ov = (typeof window!=='undefined' && window.SNPTRAIT_ORGANISM) || {};
  return Object.assign({}, TRAIT_ORGANISM._default, TRAIT_ORGANISM[fam] || {}, ov[fam] || {});
}
function traitUnit(n){ const o = TRAIT.org || TRAIT_ORGANISM._default; return n === 1 ? o.one : o.many; }

const TRAIT = {
  dataset:null, schema:TRAIT_NEUTRAL_SCHEMA, org:TRAIT_ORGANISM._default,
  rows:[], selected:new Set(), q:'', facets:{}, openGroups:new Set(), hasMeta:false,
  ranges:[],          // active numeric trait ranges [{code, min, max}] (design change 10)
  traitData:{},       // family -> parsed data/traits/<family>.traits.json | null (absent) | 'loading'
};

/* ---- numeric trait ranges from the trait side-file ----
   data/traits/<family>.traits.json (build_maize_trait_catalog.py) carries, per
   sample, GRIN evaluation summaries: numeric traits as [n_obs, mean, sd, min, max],
   coded traits as [n_obs, mode]. A range filter keeps lines whose per-accession
   MEAN lies in [min, max]; lines without a record for that trait are excluded
   while the range is active. Loaded lazily, only for families that have a file. */
function traitSideFile(){ const f = traitFamily(TRAIT.dataset); return f ? TRAIT.traitData[f] : null; }
function traitEnsureSideFile(dataset){
  const fam = traitFamily(dataset);
  if (!fam || TRAIT.traitData[fam] !== undefined || typeof fetch !== 'function') return;
  if (!traitHasSchema(dataset) || TRAIT_SCHEMA[fam]) return;      // only generated (maize) schemas ship a side-file
  TRAIT.traitData[fam] = 'loading';
  fetch(`./data/traits/${encodeURIComponent(fam)}.traits.json`)
    .then(r => { if (!r.ok) throw new Error('HTTP ' + r.status); return r.json(); })
    .then(j => { TRAIT.traitData[fam] = j; if (TRAIT.dataset === dataset && document.getElementById('traitFacets')) traitRenderFacets(); })
    .catch(() => { TRAIT.traitData[fam] = null; });
}
function traitNumericTraits(){
  const tf = traitSideFile(); if (!tf || tf === 'loading') return [];
  const n = {};
  Object.values(tf.samples || {}).forEach(s => Object.entries(s.traits || {}).forEach(([k, v]) => {
    if (Array.isArray(v) && v.length === 5 && typeof v[1] === 'number') n[k] = (n[k] || 0) + 1; }));
  return Object.entries(n).filter(([k]) => ((tf.dictionary || {})[k] || {}).scale === 'numeric')
    .map(([k, c]) => ({code: k, n: c, name: ((tf.dictionary || {})[k] || {}).name || k, unit: ((tf.dictionary || {})[k] || {}).unit || ''}))
    .sort((a, b) => a.name.localeCompare(b.name));
}
function traitValue(id, code){
  const tf = traitSideFile(); if (!tf || tf === 'loading') return null;
  const s = (tf.samples || {})[id]; const v = s && s.traits && s.traits[code];
  return (Array.isArray(v) && v.length === 5 && typeof v[1] === 'number') ? v[1] : null;
}
function traitAddRange(code, min, max){
  if (code == null){
    const sel = document.getElementById('trRangeTrait'); if (!sel) return;
    code = sel.value;
    min = document.getElementById('trRangeMin').value; max = document.getElementById('trRangeMax').value;
  }
  const lo = (min === '' || min == null) ? -Infinity : +min, hi = (max === '' || max == null) ? Infinity : +max;
  if (!code || Number.isNaN(lo) || Number.isNaN(hi)) return;
  TRAIT.ranges = TRAIT.ranges.filter(r => r.code !== code).concat([{code, min: lo, max: hi}]);
  traitRenderFacets(); traitRenderTable(); traitStatus();
}
function traitRemoveRange(code){
  TRAIT.ranges = TRAIT.ranges.filter(r => r.code !== code);
  traitRenderFacets(); traitRenderTable(); traitStatus();
}
function traitRangeStats(code){
  const vals = TRAIT.rows.map(r => traitValue(r.id, code)).filter(v => v != null).sort((a, b) => a - b);
  if (!vals.length) return null;
  const q = f => vals[Math.min(vals.length - 1, Math.floor(f * (vals.length - 1)))];
  return {n: vals.length, min: vals[0], q1: q(.25), med: q(.5), q3: q(.75), max: vals[vals.length - 1]};
}
function traitFmt(v){ return (v == null || !isFinite(v)) ? '' : (Math.abs(v) >= 100 ? v.toFixed(0) : +v.toPrecision(3)); }
function traitRangeBlockHTML(){
  const tf = traitSideFile();
  if (!tf) return '';
  if (tf === 'loading') return `<div class="facet-block"><h4>GRIN trait ranges</h4><em>Loading trait records…</em></div>`;
  const list = traitNumericTraits(); if (!list.length) return '';
  const opts = list.map(t => `<option value="${escAttrT(t.code)}">${escT(t.name)}${t.unit ? ' (' + escT(t.unit) + ')' : ''} · ${t.n}</option>`).join('');
  const active = TRAIT.ranges.map(r => {
    const t = list.find(x => x.code === r.code) || {name: r.code, unit: ''};
    const lab = `${isFinite(r.min) ? traitFmt(r.min) : '−∞'} – ${isFinite(r.max) ? traitFmt(r.max) : '∞'}${t.unit ? ' ' + t.unit : ''}`;
    return `<div class="tr-chip"><span><b>${escT(t.name)}</b> ${escT(lab)}</span>
      <button class="link-btn" onclick="traitRemoveRange('${escAttrT(r.code)}')" title="Remove this range">×</button></div>`;
  }).join('');
  return `<div class="facet-block tr-range"><h4>GRIN trait ranges</h4>
    <select id="trRangeTrait" class="tr-sel" onchange="traitRangeHint()">${opts}</select>
    <div class="tr-mm"><input id="trRangeMin" type="number" step="any" placeholder="min">
      <span>–</span><input id="trRangeMax" type="number" step="any" placeholder="max">
      <button class="qbtn" onclick="traitAddRange()">Add</button></div>
    <div class="tr-hint" id="trRangeHint"></div>
    ${active}
    <div class="tr-note">Per-accession mean of GRIN evaluation records. Lines without a record for a trait are excluded while its range is active.</div>
  </div>`;
}
function traitRangeHint(){
  const el = document.getElementById('trRangeHint'), sel = document.getElementById('trRangeTrait');
  if (!el || !sel) return;
  const st = traitRangeStats(sel.value);
  el.textContent = st ? `${st.n} lines with data · min ${traitFmt(st.min)} · median ${traitFmt(st.med)} · max ${traitFmt(st.max)}` : 'No lines with data';
}

function traitLoad(dataset){
  TRAIT.dataset = dataset;
  const sc = TRAIT.schema = traitSchema(dataset);
  TRAIT.org = traitOrganism(dataset);
  const list = Data.accessionsFor(dataset) || [];
  // fields the current schema needs on each row
  const fields = new Set([sc.groupBy, 'strain'].concat(sc.facets.map(f=>f[0]), sc.columns.map(c=>c[0]), sc.search));
  TRAIT.rows = list.map(a => {
    const row = {id:a.id, strain:a.strain || a.label || a.id};
    fields.forEach(k => { const v=a[k]; row[k] = (v!=null && v!=='') ? v : (row[k]!=null?row[k]:'Unknown'); });
    return row;
  });
  TRAIT.facets = {}; sc.facets.forEach(([k]) => TRAIT.facets[k] = new Set());
  TRAIT.q = ''; TRAIT.openGroups = new Set(); TRAIT.ranges = [];
  traitEnsureSideFile(dataset);
  TRAIT.hasMeta = TRAIT.rows.some(r => sc.facets.some(([k]) => r[k] && r[k]!=='Unknown' && r[k]!=='unknown'));
  /* Starts empty: it was seeded with SNPVersity's selection, so "replace" in the Send dialog
     sent those lines back unless they were unticked here. Kept between visits (traitLoad runs
     only when the dataset changes). */
  TRAIT.selected = new Set();
}

function traitMatch(r){
  const sc = TRAIT.schema, f = TRAIT.facets;
  for (const [k] of sc.facets){ if (f[k] && f[k].size && !f[k].has(r[k])) return false; }
  if (TRAIT.q){
    const hay = sc.search.map(k=>r[k]||'').join(' ').toLowerCase();
    if (!hay.includes(TRAIT.q)) return false;
  }
  for (const rg of TRAIT.ranges){
    const v = traitValue(r.id, rg.code);
    if (v == null || v < rg.min || v > rg.max) return false;
  }
  return true;
}
function traitVisible(){ return TRAIT.rows.filter(traitMatch); }

/* ================= PAGE ================= */
SNPTools.register('snptrait', { render: renderTrait });

function renderTrait(){
  injectTraitCSS();
  if (TRAIT.dataset !== S.dataset) traitLoad(S.dataset);
  const sc = TRAIT.schema;
  const p = document.getElementById('page');
  p.className = 'page fade';
  const dsName = (Data.datasets().find(d=>d.id===S.dataset)||{}).name || S.dataset;
  p.innerHTML = `
    <div class="sec"><div class="bar"></div><div style="width:100%">
      <div class="n">${escT(TRAIT.org.kicker)}</div>
      <h2>${escT(TRAIT.org.title)}</h2>
      <p>An interactive catalogue of the ${escT(TRAIT.org.many)} in <b>${escT(dsName)}</b>, grouped by ${escT(sc.groupLabel)}.
         Filter with the facets or the search box, batch-select, then send the selection to SNPVersity to build a VCF.</p>
    </div></div>

    <div class="ds-picker" id="traitDsPicker"></div>

    ${TRAIT.hasMeta ? '' : `<div class="trait-note">No metadata is available for this dataset yet — showing ${escT(TRAIT.org.one)} IDs only.</div>`}

    <div class="trait-shell">
      <aside class="trait-facets" id="traitFacets"></aside>
      <div class="trait-main">
        <div class="trait-toolbar">
          <div class="trait-search">${ICONS.search}
            <input id="traitSearch" type="search" placeholder="Search ${escAttrT(sc.search.slice(1,5).join(', '))}…"
                   oninput="TRAIT.q=this.value.trim().toLowerCase();traitRenderTable();traitRenderFacets();traitStatus()">
          </div>
          <div class="trait-qbtns">
            <button class="qbtn" onclick="traitRandom(.05)">5%</button>
            <button class="qbtn" onclick="traitRandom(.10)">10%</button>
            <button class="qbtn" onclick="traitRandom(.25)">25%</button>
            <button class="qbtn" onclick="traitRandom(.50)">50%</button>
            <button class="qbtn solid" onclick="traitSelectVisible(true)">Select visible</button>
            <button class="qbtn" onclick="traitSelectVisible(false)">Clear</button>
          </div>
        </div>
        <div class="trait-grid" id="traitGrid"></div>
      </div>
    </div>

    <div class="runbar" id="traitRunbar"></div>
  `;
  renderTraitDsPicker();
  traitRenderFacets();
  traitRenderTable();
  traitRenderRunbar();
}

function renderTraitDsPicker(){
  const el=document.getElementById('traitDsPicker'); if(!el) return;
  el.innerHTML = Data.datasets().map(d=>`
    <button class="ds-pill ${d.id===S.dataset?'on':''}${traitHasSchema(d.id)?'':' nometa'}" onclick="traitPickDataset('${d.id}')"
      title="${traitHasSchema(d.id)?'':'No trait / passport metadata for this dataset yet'}">
      <span class="dot"></span>${escT(d.name)} · ${d.sub}</button>`).join('');
}
function traitPickDataset(id){
  S.dataset=id; traitLoad(id);
  renderTrait();
}

/* ---- facets ---- */
function facetCounts(rows){
  const sc = TRAIT.schema, z = {};
  sc.facets.forEach(([k]) => z[k] = {});
  rows.forEach(r => sc.facets.forEach(([k]) => { const v = r[k]||'Unknown'; z[k][v] = (z[k][v]||0)+1; }));
  return z;
}
function facetLabel(k, v){
  const sc = TRAIT.schema;
  if (k===sc.groupBy && sc.groupNames[v]) return sc.groupNames[v];
  return v;
}
function traitRenderFacets(){
  const box=document.getElementById('traitFacets'); if(!box) return;
  const sc=TRAIT.schema, C=facetCounts(traitVisible());
  const blocks = sc.facets.map(([k,label])=>{
    const entries=Object.entries(C[k]).sort((a,b)=>b[1]-a[1]);
    if(!entries.length) return '';
    const sel=TRAIT.facets[k];
    const items=entries.map(([v,n])=>`
      <label class="facet-row">
        <input type="checkbox" ${sel.has(v)?'checked':''}
          onchange="traitToggleFacet('${escAttrT(k)}', this.value, this.checked)" value="${escAttrT(v)}">
        <span class="fv" title="${escAttrT(facetLabel(k,v))}">${escT(facetLabel(k,v))}</span>
        <span class="fn">${n}</span></label>`).join('');
    return `<div class="facet-block"><h4>${escT(label)}</h4>${items}</div>`;
  }).filter(Boolean).join('');
  box.innerHTML = `
    <div class="facet-head">
      <span>Filters</span>
      <button class="link-btn" onclick="traitClearFacets()">Clear all</button>
    </div>
    ${traitRangeBlockHTML()}
    ${blocks || '<div class="facet-empty">No metadata to filter.</div>'}`;
  traitRangeHint();
}
function traitToggleFacet(k,v,on){
  if(!TRAIT.facets[k]) return;
  on ? TRAIT.facets[k].add(v) : TRAIT.facets[k].delete(v);
  traitRenderTable(); traitRenderFacets(); traitStatus();
}
function traitClearFacets(){
  Object.values(TRAIT.facets).forEach(s=>s.clear()); TRAIT.ranges = [];
  traitRenderTable(); traitRenderFacets(); traitStatus();
}

/* ---- table grouped by schema.groupBy ---- */
function traitRenderTable(){
  const grid=document.getElementById('traitGrid'); if(!grid) return;
  const sc=TRAIT.schema, gkey=sc.groupBy, rows=traitVisible();
  const groups={};
  rows.forEach(r=>{ const g=r[gkey]||'Unknown'; (groups[g]=groups[g]||[]).push(r); });
  const order=sc.groupOrder||[];
  const keys=Object.keys(groups).sort((a,b)=>{
    const ia=order.indexOf(a), ib=order.indexOf(b);
    return (ia<0?99:ia)-(ib<0?99:ib) || a.localeCompare(b);
  });
  if(!rows.length){ grid.innerHTML=`<div class="trait-empty">No ${escT(TRAIT.org.many)} match the current filters.</div>`; traitStatus(); return; }
  const showMeta=TRAIT.hasMeta, cols=sc.columns;
  const rcols = TRAIT.ranges.map(rg => { const t = traitNumericTraits().find(x => x.code === rg.code); return [rg.code, t ? t.name : rg.code]; });
  grid.innerHTML = keys.map(gv=>{
    const items=groups[gv];
    const open = TRAIT.openGroups.has(gv) || TRAIT.q || TRAIT.ranges.length || Object.values(TRAIT.facets).some(s=>s.size);
    const selN = items.reduce((n,r)=>n+(TRAIT.selected.has(r.id)?1:0),0);
    const color = (sc.groupColors&&sc.groupColors[gv])||'#888';
    const gname = (sc.groupNames&&sc.groupNames[gv])||gv;
    const head = `<div class="tg-head" onclick="traitToggleGroup('${escAttrT(gv)}')" style="border-left:4px solid ${color}">
        <span class="caret ${open?'op':''}">${ICONS.caret}</span>
        <b>${escT(gname)}</b>
        <span class="tg-count">${items.length} ${escT(traitUnit(items.length))}</span>
        <span class="tg-sel">${selN} selected</span>
        <span class="tg-acts">
          <button class="qbtn" onclick="event.stopPropagation();traitSelectGroup('${escAttrT(gv)}',true)">all</button>
          <button class="qbtn" onclick="event.stopPropagation();traitSelectGroup('${escAttrT(gv)}',false)">none</button>
        </span></div>`;
    if(!open) return `<div class="tg">${head}</div>`;
    const rowsHtml = items.map(r=>`
      <tr class="${TRAIT.selected.has(r.id)?'on':''}" onclick="traitToggle('${escAttrT(r.id)}')">
        <td class="tcb"><span class="cb">${ICONS.check}</span></td>
        <td class="mono">${escT(r.id)}</td>
        ${showMeta?cols.map(([k])=>`<td>${escT(r[k]==null?'':r[k])}</td>`).join(''):''}
        ${rcols.map(([c])=>`<td class="num">${escT(traitFmt(traitValue(r.id, c)))}</td>`).join('')}
      </tr>`).join('');
    const thead = showMeta
      ? `<tr><th></th><th>ID</th>${cols.map(([,lbl])=>`<th>${escT(lbl)}</th>`).join('')}${rcols.map(([,lbl])=>`<th>${escT(lbl)} (mean)</th>`).join('')}</tr>`
      : `<tr><th></th><th>ID</th></tr>`;
    return `<div class="tg">${head}<table class="tg-table"><thead>${thead}</thead><tbody>${rowsHtml}</tbody></table></div>`;
  }).join('');
  traitStatus();
}
function traitToggleGroup(gv){
  TRAIT.openGroups.has(gv) ? TRAIT.openGroups.delete(gv) : TRAIT.openGroups.add(gv);
  traitRenderTable();
}
function traitToggle(id){
  TRAIT.selected.has(id) ? TRAIT.selected.delete(id) : TRAIT.selected.add(id);
  traitRenderTable(); traitRenderRunbar();
}
function traitSelectGroup(gv,on){
  const gkey=TRAIT.schema.groupBy;
  TRAIT.rows.filter(r=>(r[gkey]||'Unknown')===gv && traitMatch(r)).forEach(r=> on?TRAIT.selected.add(r.id):TRAIT.selected.delete(r.id));
  traitRenderTable(); traitRenderFacets(); traitRenderRunbar();
}
function traitSelectVisible(on){
  traitVisible().forEach(r=> on?TRAIT.selected.add(r.id):TRAIT.selected.delete(r.id));
  traitRenderTable(); traitRenderRunbar();
}
function traitRandom(frac){
  const vis=traitVisible().slice().sort(()=>Math.random()-0.5);
  const n=Math.max(1, Math.round(vis.length*frac));
  vis.forEach(r=>TRAIT.selected.delete(r.id));
  for(let i=0;i<n;i++) TRAIT.selected.add(vis[i].id);
  traitRenderTable(); traitRenderRunbar();
}

function traitStatus(){
  const el=document.getElementById('traitStatus'); if(!el) return;
  el.textContent = `Visible ${traitVisible().length} / ${TRAIT.rows.length} · Selected ${TRAIT.selected.size}`;
}

/* ---- run bar: export + hand off to SNPVersity ---- */
function traitRenderRunbar(){
  const rb=document.getElementById('traitRunbar'); if(!rb) return;
  const n=TRAIT.selected.size;
  rb.innerHTML=`
    <div class="summ">
      <span class="kick${n?'':' wait'}">${n?'Ready':'Select '+escT(TRAIT.org.many)+' to continue'}</span>
      <b id="traitStatus"></b>
    </div>
    <div class="spacer"></div>
    <button class="btn" onclick="traitExport('csv')">${ICONS.download} CSV</button>
    <button class="btn" onclick="traitExport('json')">${ICONS.download} JSON</button>
    <button class="btn primary" ${n?'':'disabled'} onclick="traitOpenSend()" aria-haspopup="dialog">
      ${ICONS.dna} Send ${n} ${escT(traitUnit(n))} to SNPVersity…</button>`;
  traitStatus();
}
/* ---- hand-off dialog: replace SNPVersity's selection or add to it ----
   The same choice GWAS Explorer's "Send to SNPVersity" offers: two mutually exclusive
   boxes, "replace" pre-ticked, Send disabled until one is ticked. The add box says how
   many of the lines are new to SNPVersity's selection. */
function traitSendCounts(){
  const ids=[...TRAIT.selected];
  const cur=(typeof S!=='undefined' && S.selected instanceof Set) ? S.selected : new Set();
  const fresh=ids.filter(id=>!cur.has(id)).length;
  return {n:ids.length, current:cur.size, fresh, already:ids.length-fresh};
}
function traitOpenSend(){
  if(!TRAIT.selected.size) return;
  traitCloseSend();
  const c=traitSendCounts(), unit=escT(traitUnit(c.n));
  const accs=k=>`accession${k===1?'':'s'}`;
  const host=document.createElement('div'); host.id='traitSendDlg';
  host.innerHTML=`<div class="trait-send-backdrop" onclick="traitCloseSend()"></div>
    <div class="trait-send" role="dialog" aria-modal="true" aria-labelledby="traitSendTitle">
      <div class="trait-send-head"><h3 id="traitSendTitle">Send to SNPVersity</h3>
        <button type="button" class="trait-send-x" onclick="traitCloseSend()" aria-label="Close">&times;</button></div>
      <div class="trait-send-body">
        <div class="trait-send-what"><b>${c.n}</b> ${unit} selected in SNPTrait</div>
        <label class="trait-send-check"><input type="checkbox" id="traitSendReplace" checked onchange="traitSendPick(this)">
          <span>Send ${unit} — replace the <b>${c.current}</b> ${accs(c.current)} currently selected in SNPVersity</span></label>
        <label class="trait-send-check"><input type="checkbox" id="traitSendAdd" onchange="traitSendPick(this)">
          <span>Send ${unit} — keep the <b>${c.current}</b> ${accs(c.current)} currently selected in SNPVersity and add to them
            <span class="trait-send-n">(${c.fresh} new${c.already?`, ${c.already} already selected`:''})</span></span></label>
        <div class="trait-send-hint" id="traitSendHint"></div>
      </div>
      <div class="trait-send-foot">
        <button type="button" class="btn" onclick="traitCloseSend()">Cancel</button>
        <button type="button" class="btn primary" id="traitSendGo" onclick="traitSendConfirm()">Send</button>
      </div>
    </div>`;
  document.body.appendChild(host);
  document.addEventListener('keydown', traitSendKey);
  const first=document.getElementById('traitSendReplace'); if(first) first.focus();
}
function traitSendKey(e){ if(e.key==='Escape') traitCloseSend(); }
function traitCloseSend(){
  const el=document.getElementById('traitSendDlg'); if(el) el.remove();
  document.removeEventListener('keydown', traitSendKey);
}
/* the two boxes exclude each other; Send needs one of them */
function traitSendPick(box){
  const rep=document.getElementById('traitSendReplace'), add=document.getElementById('traitSendAdd');
  if(box && box.checked){ (box===rep?add:rep).checked=false; }
  const any=rep.checked||add.checked;
  document.getElementById('traitSendGo').disabled=!any;
  document.getElementById('traitSendHint').textContent=any?'':'Choose whether to replace SNPVersity’s selection or add to it.';
}
function traitSendConfirm(){
  const rep=document.getElementById('traitSendReplace'), add=document.getElementById('traitSendAdd');
  const mode=rep&&rep.checked ? 'replace' : (add&&add.checked ? 'add' : null);
  if(!mode) return;
  traitCloseSend();
  traitSendToVersity(mode);
}
/* mode: 'replace' (the default, as before the dialog) or 'add' */
function traitSendToVersity(mode){
  const ids=[...TRAIT.selected];
  if(!ids.length) return;
  const merge = mode==='add' ? 'add' : 'replace';
  const payload={dataset:TRAIT.dataset||S.dataset, accessions:ids, merge,
                 from:'SNPTrait', note:`${ids.length} ${traitUnit(ids.length)} selected in SNPTrait`};
  if(typeof Handoff!=='undefined' && Handoff.toVersity){ Handoff.toVersity(payload); return; }
  if(typeof window.versityRequest==='function'){ window.versityRequest(payload); return; }
  S.selected = merge==='add' ? new Set([...S.selected, ...ids]) : new Set(ids);
  go('snpversity');
}
function traitExport(kind){
  const sel=TRAIT.rows.filter(r=>TRAIT.selected.has(r.id));
  if(!sel.length){ return; }
  const sc=TRAIT.schema;
  const hdr=[...new Set(['id', sc.groupBy].concat(sc.columns.map(c=>c[0]), sc.facets.map(f=>f[0])))];
  TRAIT.ranges.forEach(rg => { hdr.push(rg.code); sel.forEach(r => { r[rg.code] = traitValue(r.id, rg.code); }); });
  let blob, name;
  if(kind==='json'){
    blob=new Blob([JSON.stringify(sel.map(r=>{const o={};hdr.forEach(h=>o[h]=r[h]);return o;}),null,2)],{type:'application/json'});
    name=TRAIT.org.exportBase+'.json';
  } else {
    const esc=v=>{v=String(v==null?'':v); return /[",\n]/.test(v)?`"${v.replace(/"/g,'""')}"`:v;};
    const csv=[hdr.join(',')].concat(sel.map(r=>hdr.map(h=>esc(r[h])).join(','))).join('\n');
    blob=new Blob([csv],{type:'text/csv'}); name=TRAIT.org.exportBase+'.csv';
  }
  const a=document.createElement('a'); a.href=URL.createObjectURL(blob); a.download=name;
  document.body.appendChild(a); a.click(); a.remove(); setTimeout(()=>URL.revokeObjectURL(a.href),1500);
}

/* ---- helpers ---- */
function escT(s){return String(s==null?'':s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');}
function escAttrT(s){return escT(s).replace(/"/g,'&quot;').replace(/'/g,'&#39;');}

function injectTraitCSS(){
  if(document.getElementById('snptrait-css')) return;
  const s=document.createElement('style'); s.id='snptrait-css';
  s.textContent=`
  .ds-picker{display:flex;gap:8px;flex-wrap:wrap;margin:10px 0 4px}
  .ds-pill{display:inline-flex;align-items:center;gap:7px;border:1px solid var(--line);background:#fff;border-radius:20px;
    padding:7px 13px;font:600 12.5px/1 var(--body);color:var(--ink);cursor:pointer}
  .ds-pill .dot{width:8px;height:8px;border-radius:50%;background:#c3cbd8}
  .ds-pill.on{border-color:var(--blue-600);color:var(--blue-600);background:#eef4ff}
  .ds-pill.on .dot{background:var(--blue-600)}
  .trait-note{background:#fff8ec;border:1px solid #f0dcae;color:#8a6d1e;border-radius:9px;padding:9px 12px;font-size:12.5px;margin:10px 0}
  .trait-shell{display:grid;grid-template-columns:230px 1fr;gap:16px;margin-top:14px;align-items:start}
  @media(max-width:820px){.trait-shell{grid-template-columns:1fr}}
  .trait-facets{border:1px solid var(--line);border-radius:12px;background:#fff;padding:10px 12px;position:sticky;top:10px;max-height:80vh;overflow:auto}
  .facet-head{display:flex;align-items:center;justify-content:space-between;font-weight:700;font-size:12.5px;margin-bottom:6px}
  .facet-block{margin:10px 0}
  .facet-block h4{margin:0 0 5px;font-size:11px;letter-spacing:.04em;text-transform:uppercase;color:var(--muted)}
  .facet-row{display:flex;align-items:center;gap:7px;font-size:12px;padding:2px 0;cursor:pointer}
  .facet-row .fv{flex:1;min-width:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
  .facet-row .fn{color:var(--faint);font-family:var(--mono);font-size:11px}
  .facet-empty,.facet-block em{color:var(--faint);font-size:12px}
  .trait-toolbar{display:flex;gap:10px;flex-wrap:wrap;align-items:center;margin-bottom:10px}
  .trait-search{flex:1;min-width:200px;display:flex;align-items:center;gap:7px;border:1px solid var(--line);border-radius:9px;padding:7px 10px;background:#fff}
  .trait-search svg{width:16px;height:16px;color:var(--muted)}
  .trait-search input{border:0;outline:0;flex:1;font-size:13px;background:transparent;color:var(--ink)}
  .trait-qbtns{display:flex;gap:6px;flex-wrap:wrap}
  .qbtn{border:1px solid var(--line);background:#fff;border-radius:7px;padding:6px 10px;font:600 12px/1 var(--body);color:var(--ink);cursor:pointer}
  .qbtn.solid{background:var(--blue-600,#2563eb);border-color:var(--blue-600,#2563eb);color:#fff}
  .trait-grid{display:flex;flex-direction:column;gap:8px}
  .ds-pill.nometa{opacity:.55}
  .tr-range .tr-sel{width:100%;font-size:12px;padding:5px;border:1px solid var(--line);border-radius:7px;background:#fff}
  .tr-range .tr-mm{display:flex;gap:5px;align-items:center;margin-top:6px}
  .tr-range .tr-mm input{width:70px;font-size:12px;padding:4px 6px;border:1px solid var(--line);border-radius:6px}
  .tr-range .tr-hint,.tr-range .tr-note{font-size:11px;color:var(--muted);margin-top:5px;line-height:1.35}
  .tr-chip{display:flex;justify-content:space-between;align-items:center;gap:6px;margin-top:6px;padding:4px 8px;border:1px solid #cfe3d6;background:#f3faf5;border-radius:8px;font-size:12px}
  .tg-table td.num{text-align:right;font-family:var(--mono)}
  .tg{border:1px solid var(--line);border-radius:10px;overflow:hidden;background:#fff}
  .tg-head{display:flex;align-items:center;gap:9px;padding:9px 11px;cursor:pointer;background:#fbfcfe}
  .tg-head .caret{transition:transform .15s;display:inline-flex;color:var(--muted)}
  .tg-head .caret.op{transform:rotate(90deg)}
  .tg-head .caret svg{width:14px;height:14px}
  .tg-count{color:var(--muted);font-size:12px}
  .tg-sel{color:var(--blue-600);font-size:11.5px;font-weight:600}
  .tg-acts{margin-left:auto;display:flex;gap:6px}
  .tg-table{width:100%;border-collapse:collapse;font-size:12.5px}
  .tg-table thead th{text-align:left;padding:6px 10px;color:var(--muted);font-size:11px;text-transform:uppercase;letter-spacing:.03em;border-top:1px solid var(--line);background:#fafbfd}
  .tg-table tbody td{padding:6px 10px;border-top:1px solid #eef1f6}
  .tg-table tbody tr{cursor:pointer}
  .tg-table tbody tr:hover{background:#f5f8ff}
  .tg-table tbody tr.on{background:#eaf2ff}
  .tg-table .tcb{width:26px}
  .tg-table .cb{display:inline-flex;width:15px;height:15px;border:1.5px solid #c3ccda;border-radius:4px;color:#fff}
  .tg-table .cb svg{width:12px;height:12px;opacity:0}
  .tg-table tr.on .cb{background:var(--blue-600);border-color:var(--blue-600)}
  .tg-table tr.on .cb svg{opacity:1}
  .trait-empty{padding:26px;text-align:center;color:var(--faint);border:1px dashed var(--line);border-radius:10px}
  /* hand-off dialog (same shape as GWAS Explorer's Send to SNPVersity) */
  .trait-send-backdrop{position:fixed;inset:0;background:rgba(10,15,28,.35);z-index:43}
  .trait-send{position:fixed;top:50%;left:50%;transform:translate(-50%,-50%);width:min(480px,calc(100vw - 32px));max-height:85vh;
    background:#fff;border:1px solid var(--line);border-radius:var(--rad);box-shadow:var(--shadow-lg);z-index:44;
    display:flex;flex-direction:column;overflow:hidden}
  .trait-send-head{display:flex;align-items:center;justify-content:space-between;padding:16px 18px 12px;border-bottom:1px solid var(--line)}
  .trait-send-head h3{font-family:var(--disp);font-size:15px;font-weight:700;margin:0;color:var(--ink)}
  .trait-send-x{border:0;background:none;font-size:22px;line-height:1;color:var(--muted);cursor:pointer;padding:0 4px}
  .trait-send-body{padding:16px 18px;display:flex;flex-direction:column;gap:14px;overflow:auto}
  .trait-send-what{font-size:13px;color:var(--ink)} .trait-send-what b{font-family:var(--mono)}
  .trait-send-check{display:flex;gap:9px;align-items:flex-start;font-size:13px;line-height:1.45;cursor:pointer;color:var(--ink)}
  .trait-send-check input{margin:0;position:relative;top:3px;cursor:pointer;flex:0 0 auto}
  .trait-send-check b{font-family:var(--mono)}
  .trait-send-n{color:var(--muted);white-space:nowrap}
  .trait-send-hint{font-size:11.5px;color:#b45309;min-height:14px}
  .trait-send-foot{padding:12px 18px 16px;border-top:1px solid var(--line);display:flex;justify-content:flex-end;gap:8px}
  `;
  document.head.appendChild(s);
}
