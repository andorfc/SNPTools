/* =====================================================================
 *  snptrait.js — SNPTrait = metadata "Strain Selector".
 *  Registers 'snptrait'. Lets the user filter and select isolates by
 *  metadata, then hand the selected set to SNPVersity to build a VCF.
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
function traitSchema(dataset){
  const fam = (typeof Data!=='undefined' && Data.familyOf) ? Data.familyOf(dataset) : 'graminearum';
  return fam==='graminearum' ? TRAIT_SCHEMA.graminearum : TRAIT_SCHEMA.vert;
}

const TRAIT = {
  dataset:null, schema:TRAIT_SCHEMA.graminearum,
  rows:[], selected:new Set(), q:'', facets:{}, openGroups:new Set(), hasMeta:false,
};

function traitLoad(dataset){
  TRAIT.dataset = dataset;
  const sc = TRAIT.schema = traitSchema(dataset);
  const list = Data.accessionsFor(dataset) || [];
  // fields the current schema needs on each row
  const fields = new Set([sc.groupBy, 'strain'].concat(sc.facets.map(f=>f[0]), sc.columns.map(c=>c[0]), sc.search));
  TRAIT.rows = list.map(a => {
    const row = {id:a.id, strain:a.strain || a.label || a.id};
    fields.forEach(k => { const v=a[k]; row[k] = (v!=null && v!=='') ? v : (row[k]!=null?row[k]:'Unknown'); });
    return row;
  });
  TRAIT.facets = {}; sc.facets.forEach(([k]) => TRAIT.facets[k] = new Set());
  TRAIT.q = ''; TRAIT.openGroups = new Set();
  TRAIT.hasMeta = TRAIT.rows.some(r => sc.facets.some(([k]) => r[k] && r[k]!=='Unknown' && r[k]!=='unknown'));
  TRAIT.selected = new Set([...S.selected].filter(id => TRAIT.rows.some(r => r.id === id)));
}

function traitMatch(r){
  const sc = TRAIT.schema, f = TRAIT.facets;
  for (const [k] of sc.facets){ if (f[k] && f[k].size && !f[k].has(r[k])) return false; }
  if (TRAIT.q){
    const hay = sc.search.map(k=>r[k]||'').join(' ').toLowerCase();
    if (!hay.includes(TRAIT.q)) return false;
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
      <div class="n">STRAIN SELECTOR</div>
      <h2>Select Fusarium isolates by metadata</h2>
      <p>An interactive catalogue of the isolates in <b>${escT(dsName)}</b>, grouped by ${escT(sc.groupLabel)}.
         Filter with the facets or the search box, batch-select, then send the selection to SNPVersity to build a VCF.</p>
    </div></div>

    <div class="ds-picker" id="traitDsPicker"></div>

    ${TRAIT.hasMeta ? '' : `<div class="trait-note">No metadata is available for this reference yet — showing isolate IDs only.</div>`}

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
    <button class="ds-pill ${d.id===S.dataset?'on':''}" onclick="traitPickDataset('${d.id}')">
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
    ${blocks || '<div class="facet-empty">No metadata to filter.</div>'}`;
}
function traitToggleFacet(k,v,on){
  if(!TRAIT.facets[k]) return;
  on ? TRAIT.facets[k].add(v) : TRAIT.facets[k].delete(v);
  traitRenderTable(); traitRenderFacets(); traitStatus();
}
function traitClearFacets(){
  Object.values(TRAIT.facets).forEach(s=>s.clear());
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
  if(!rows.length){ grid.innerHTML='<div class="trait-empty">No isolates match the current filters.</div>'; traitStatus(); return; }
  const showMeta=TRAIT.hasMeta, cols=sc.columns;
  grid.innerHTML = keys.map(gv=>{
    const items=groups[gv];
    const open = TRAIT.openGroups.has(gv) || TRAIT.q || Object.values(TRAIT.facets).some(s=>s.size);
    const selN = items.reduce((n,r)=>n+(TRAIT.selected.has(r.id)?1:0),0);
    const color = (sc.groupColors&&sc.groupColors[gv])||'#888';
    const gname = (sc.groupNames&&sc.groupNames[gv])||gv;
    const head = `<div class="tg-head" onclick="traitToggleGroup('${escAttrT(gv)}')" style="border-left:4px solid ${color}">
        <span class="caret ${open?'op':''}">${ICONS.caret}</span>
        <b>${escT(gname)}</b>
        <span class="tg-count">${items.length} isolate${items.length!==1?'s':''}</span>
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
      </tr>`).join('');
    const thead = showMeta
      ? `<tr><th></th><th>ID</th>${cols.map(([,lbl])=>`<th>${escT(lbl)}</th>`).join('')}</tr>`
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
      <span class="kick${n?'':' wait'}">${n?'Ready':'Select isolates to continue'}</span>
      <b id="traitStatus"></b>
    </div>
    <div class="spacer"></div>
    <button class="btn" onclick="traitExport('csv')">${ICONS.download} CSV</button>
    <button class="btn" onclick="traitExport('json')">${ICONS.download} JSON</button>
    <button class="btn primary" ${n?'':'disabled'} onclick="traitSendToVersity()">
      ${ICONS.dna} Send ${n} isolate${n===1?'':'s'} to SNPVersity</button>`;
  traitStatus();
}
function traitSendToVersity(){
  const ids=[...TRAIT.selected];
  if(!ids.length) return;
  if(typeof window.versityRequest==='function'){
    window.versityRequest({
      dataset:S.dataset, accessions:ids, merge:'replace',
      from:'SNPTrait Strain Selector', note:`${ids.length} isolate${ids.length===1?'':'s'} from the Strain Selector`,
    });
  } else {
    S.selected=new Set(ids); go('snpversity');
  }
}
function traitExport(kind){
  const sel=TRAIT.rows.filter(r=>TRAIT.selected.has(r.id));
  if(!sel.length){ return; }
  const sc=TRAIT.schema;
  const hdr=[...new Set(['id', sc.groupBy].concat(sc.columns.map(c=>c[0]), sc.facets.map(f=>f[0])))];
  let blob, name;
  if(kind==='json'){
    blob=new Blob([JSON.stringify(sel.map(r=>{const o={};hdr.forEach(h=>o[h]=r[h]);return o;}),null,2)],{type:'application/json'});
    name='selected_strains.json';
  } else {
    const esc=v=>{v=String(v==null?'':v); return /[",\n]/.test(v)?`"${v.replace(/"/g,'""')}"`:v;};
    const csv=[hdr.join(',')].concat(sel.map(r=>hdr.map(h=>esc(r[h])).join(','))).join('\n');
    blob=new Blob([csv],{type:'text/csv'}); name='selected_strains.csv';
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
  `;
  document.head.appendChild(s);
}
