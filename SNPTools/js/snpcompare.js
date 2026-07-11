/* =====================================================================
 *  snpcompare.js — SNPCompare: identity-by-state similarity viewer.
 *
 *  Ranks every MaizeGDB accession by similarity to a chosen focal
 *  accession, in two scopes:
 *    • Genome-wide (global): precomputed, served per-focal-accession by
 *      ibsCompare.php  (only the focal row is sent to the browser).
 *    • This region (local): computed in-browser from a SNPVersity result
 *      (the generated VCF matrix), using the same definitions.
 *    • Both: global vs local side-by-side, with Δ = local − global, to
 *      surface region-specific relatedness (introgression / selection).
 *
 *  Similarity  = fraction of co-called SNP sites (both non-missing) that
 *                have the same genotype.
 *  Missing%    = percent of SNP sites where at least one of the pair is
 *                missing a genotype.
 *
 *  Metadata columns (Project / SRA ID / Accession Name) are joined from
 *  the accession catalog by ID, so the precomputed file only needs the
 *  numbers (id, similarity, missing).
 *
 *  Depends on core.js (SNPTools, S, go, ICONS) and data.js (Data).
 * ===================================================================== */
const SNPCompare = (function () {

  const CFG = {
    globalEndpoint : 'ibsCompare.php', // ?focal=<ID> -> {rows:[{id,similarity,missing}]}
    useDemoGlobal  : false,            // real data via ibsCompare.php (set true to preview without a backend)
    defaultDataset : 'mgdb2026_hq',
  };

  const ST = {
    mode:'global', dataset:CFG.defaultDataset, focal:null,
    sortKey:'sim', sortDir:-1,
    fSimMin:null, fSimMax:null, fMissMax:null, fProj:'all',
    input:null,            // local hand-off {rows, accs, chr, start, end, dataset,...}
    globalCache:{},        // focalId -> {rows, demo}
    allRows:[], gdemo:false, _meta:null, _metaDs:null,
  };

  /* ---------------- genotype + local IBS ---------------- */
  function dose(g){ if(g==null)return null; if(g==='0/0')return 0; if(g==='1/1')return 2;
    if(g==='./.'||g==='.'||g==='')return null; return 1; }
  function localCompute(input, focalId){
    if(!input) return null;
    const accs=input.accs, n=accs.length, fi=accs.findIndex(a=>a.id===focalId);
    if(fi<0) return null;
    const total=input.rows.length, out=[];
    for(let o=0;o<n;o++){ let both=0,match=0,either=0;
      for(const r of input.rows){ const gf=dose(r.gts[fi]), go=dose(r.gts[o]);
        if(gf==null||go==null) either++; else { both++; if(gf===go) match++; } }
      out.push({id:accs[o].id, sim: both? match/both : null, miss: total?100*either/total:0, both});
    }
    return out;
  }

  /* ---------------- metadata join ---------------- */
  function metaMap(ds){
    if(ST._meta && ST._metaDs===ds) return ST._meta;
    const projBio={};
    (Data.projectsFor(ds)||[]).forEach(p=>{ projBio[p.id]=(p.bioprojects&&p.bioprojects.length)?p.bioprojects.join(', '):(p.title||''); });
    const m={};
    (Data.accessionsFor(ds)||[]).forEach(a=>{ m[a.id]={id:a.id, name:a.founder, run:a.run, proj:a.proj, projColor:a.projColor, bio:projBio[a.proj]||''}; });
    ST._meta=m; ST._metaDs=ds; return m;
  }
  function projectOptions(ds){
    const seen={}, out=[];
    (Data.projectsFor(ds)||[]).forEach(p=>{ const b=(p.bioprojects&&p.bioprojects.length)?p.bioprojects.join(', '):(p.title||p.id);
      if(!seen[b]){seen[b]=1; out.push({id:p.id, label:b});} });
    return out;
  }

  /* ---------------- global source (endpoint + demo fallback) ---------------- */
  function hash(s){ let h=2166136261>>>0; for(let i=0;i<s.length;i++){ h^=s.charCodeAt(i); h=Math.imul(h,16777619);} return h>>>0; }
  function globalDemo(ds, focalId){
    return (Data.accessionsFor(ds)||[]).map(a=>{
      if(a.id===focalId) return {id:a.id, sim:1, miss:0};
      const h=hash(focalId+'|'+a.id);
      return {id:a.id,
        sim:+(0.85+(h%1500)/1500*0.149).toFixed(4),   // 0.85–0.999
        miss:+(((h>>>7)%3000)/100).toFixed(2)};        // 0–30, unsigned shift (no negatives)
    });
  }
  async function getGlobal(ds, focalId){
    if(ST.globalCache[focalId]) return ST.globalCache[focalId];
    let res;
    if(CFG.useDemoGlobal){ res={rows:globalDemo(ds,focalId), demo:true}; }
    else {
      const resp=await fetch(`${CFG.globalEndpoint}?focal=${encodeURIComponent(focalId)}`,{cache:'no-store'});
      if(!resp.ok) throw new Error('ibsCompare.php failed (HTTP '+resp.status+')');
      const raw=await resp.text(); let j;
      try{ j=JSON.parse(raw); }catch(e){ throw new Error('ibsCompare.php did not return JSON:\n'+raw.slice(0,500)); }
      const rows=(j.rows|| (j.ids? j.ids.map((id,i)=>({id, similarity:j.sim[i], missing:j.miss[i]})):[]))
        .map(r=>({id:r.id, sim:+(r.similarity!=null?r.similarity:r.sim), miss:Math.max(0,+(r.missing!=null?r.missing:r.miss))}));
      res={rows, demo:false};
    }
    ST.globalCache[focalId]=res; return res;
  }

  /* ---------------- build the combined row set ---------------- */
  async function buildRows(){
    const ds=ST.dataset, meta=metaMap(ds);
    let g=null, l=null;
    if(ST.mode==='global'||ST.mode==='both') g=await getGlobal(ds, ST.focal);
    if(ST.mode==='local' ||ST.mode==='both') l=localCompute(ST.input, ST.focal);
    const byId={};
    const add=id=>{ if(!byId[id]) byId[id]=Object.assign({id}, meta[id]||{name:id,run:id,bio:''}); return byId[id]; };
    if(g) g.rows.forEach(r=>{ const x=add(r.id); x.gsim=r.sim; x.gmiss=r.miss; });
    if(l) l.forEach(r=>{ const x=add(r.id); x.lsim=r.sim; x.lmiss=r.miss; x.lboth=r.both; });
    let rows=Object.values(byId);
    if(ST.mode==='local') rows=rows.filter(r=>r.lsim!=null || r.id===ST.focal);
    if(ST.mode==='both')  rows=rows.filter(r=>r.lsim!=null); // intersection = region accessions
    rows.forEach(r=>{ r.dsim=(r.lsim!=null&&r.gsim!=null)? r.lsim-r.gsim : null;
      r.sim = ST.mode==='local'? r.lsim : r.gsim; r.miss = ST.mode==='local'? r.lmiss : r.gmiss; });
    ST.allRows=rows; ST.gdemo=!!(g&&g.demo);
  }

  /* ---------------- sort + filter ---------------- */
  function viewRows(){
    let rows=ST.allRows.slice();
    const num=v=>v==null?-Infinity:v;
    if(ST.fSimMin!=null) rows=rows.filter(r=>num(r.sim)>=ST.fSimMin);
    if(ST.fSimMax!=null) rows=rows.filter(r=>num(r.sim)<=ST.fSimMax);
    if(ST.fMissMax!=null)rows=rows.filter(r=>(r.miss==null?0:r.miss)<=ST.fMissMax);
    if(ST.fProj!=='all') rows=rows.filter(r=>r.proj===ST.fProj);
    const k=ST.sortKey, d=ST.sortDir;
    rows.sort((a,b)=>{
      let va=a[k], vb=b[k];
      if(k==='id'||k==='name'||k==='run'||k==='bio'){ va=(va||'').toString(); vb=(vb||'').toString(); return d*va.localeCompare(vb,undefined,{numeric:true}); }
      va=va==null?-Infinity:va; vb=vb==null?-Infinity:vb; return d*(va-vb);
    });
    return rows;
  }

  /* ---------------- render ---------------- */
  // Genome-wide precomputed matrices exist only for the MaizeGDB 2026 accession set.
  function globalAvailable(){ try{ return Data.familyOf(ST.dataset)==='mgdb2026'; }catch(e){ return true; } }

  function render(){
    injectCSS();
    const crumb=document.getElementById('crumbTool'); if(crumb) crumb.innerHTML='<b>SNPCompare</b>';
    const page=document.getElementById('page');
    const prevDs=ST.dataset;
    // pick up a hand-off from SNPVersity
    if(S.compareInput){ ST.input=S.compareInput; }
    if(ST.input){ ST.dataset=ST.input.dataset||ST.dataset; }
    if(ST.dataset!==prevDs){ ST.focal=null; ST._meta=null; ST.globalCache={}; ST.allRows=[]; }
    const gAvail=globalAvailable();
    // sensible defaults
    if(!ST.focal){
      if(ST.input && ST.input.accs && ST.input.accs.length){
        ST.focal=ST.input.accs[0].id; ST.mode = gAvail ? 'both' : 'local';
      } else {
        const a=(Data.accessionsFor(ST.dataset)||[])[0]; ST.focal=a?a.id:null; ST.mode = gAvail ? 'global' : 'local';
      }
    }
    if(!gAvail && ST.mode!=='local') ST.mode='local';   // no matrix -> region scope only
    // when arriving from another tool, ensure the focal accession exists in this dataset;
    // if not, fall back to the first accession selected in SNPVersity (or first in dataset)
    if(ST.focal){
      const inRegion = ST.input && ST.input.accs && ST.input.accs.some(a=>a.id===ST.focal);
      const inDataset = !!(Data.accessionById && Data.accessionById(ST.focal));
      const ok = ST.input ? (inRegion || (gAvail && inDataset)) : inDataset;
      if(!ok){
        ST.focal = (ST.input && ST.input.accs && ST.input.accs.length) ? ST.input.accs[0].id
                 : ((Data.accessionsFor(ST.dataset)||[])[0]||{}).id || null;
      }
    }
    page.className='page fade';
    page.innerHTML=shell();
    recompute();
  }

  function shell(){
    const ds=ST.dataset;
    const gAvail=globalAvailable();
    const hasLocal=!!(ST.input&&ST.input.accs&&ST.input.accs.length);
    // focal picker lists genome-wide accessions when a matrix exists, else the region's accessions
    const ids=gAvail ? (Data.accessionsFor(ds)||[]).map(a=>a.id)
                     : (hasLocal ? ST.input.accs.map(a=>a.id) : []);
    const region=hasLocal?`${ST.input.chr}:${(+ST.input.start).toLocaleString()}–${(+ST.input.end).toLocaleString()}`:null;
    const projOpts=projectOptions(ds);
    const dsName=ST.input?(ST.input.datasetName||ST.dataset):ST.dataset;
    return `
    <section class="sec"><div class="bar"></div><div style="width:100%">
      <h1>SNPCompare</h1>
      <p>Rank accessions by identity-by-state similarity to a focal accession — genome-wide and within a queried region.</p>
    </div></section>

    <div class="card pad" style="margin-bottom:16px">
      <div style="display:flex;gap:22px;flex-wrap:wrap;align-items:flex-end">
        <div style="min-width:480px;flex:1 1 480px">
          <div class="fl-lbl">Focal accession</div>
          <div style="display:flex;gap:8px">
            <input id="cmpFocal" list="cmpFocalList" value="${esc(ST.focal||'')}" placeholder="type a SNPVersity ID…"
              style="flex:1;min-width:340px;border:1px solid var(--line);border-radius:9px;padding:9px 11px;font-family:var(--mono);font-size:13px">
            <datalist id="cmpFocalList">${ids.slice(0,4000).map(i=>`<option value="${esc(i)}">`).join('')}</datalist>
            <button class="btn solid" onclick="SNPCompare.setFocalFromInput()">Show similarities</button>
          </div>
        </div>
        <div>
          <div class="fl-lbl">Scope</div>
          ${scopeBtn('global','Genome-wide',gAvail,'Genome-wide matrix available for MaizeGDB 2026 only')}
          ${scopeBtn('local','This region',hasLocal,'Send a result from SNPVersity to enable')}
          ${scopeBtn('both','Both (Δ)',gAvail&&hasLocal, !gAvail?'Genome-wide matrix available for MaizeGDB 2026 only':'Send a result from SNPVersity to enable')}
        </div>
        ${region?`<div><div class="fl-lbl">Region</div><div class="c-mono" style="color:var(--blue-600);font-size:13px;padding:8px 0">${region}</div></div>`:''}
      </div>
      ${!gAvail?`<div class="mtx-note" style="margin-top:12px">Genome-wide precomputed IBS is available for <b>MaizeGDB 2026</b> only. For <b>${esc(dsName)}</b>, use <b>This region</b> scope — SNPCompare computes identity-by-state live from your SNPVersity result (no prebuilt matrix needed).</div>`:''}

      <div style="display:flex;gap:14px;flex-wrap:wrap;align-items:flex-end;margin-top:14px">
        <div><div class="fl-lbl">Similarity ≥</div><input class="cmp-f" id="fSimMin" type="number" step="0.01" min="0" max="1" value="${ST.fSimMin??''}" oninput="SNPCompare.setF('fSimMin',this.value)"></div>
        <div><div class="fl-lbl">Similarity ≤</div><input class="cmp-f" id="fSimMax" type="number" step="0.01" min="0" max="1" value="${ST.fSimMax??''}" oninput="SNPCompare.setF('fSimMax',this.value)"></div>
        <div><div class="fl-lbl">Missing% ≤</div><input class="cmp-f" id="fMissMax" type="number" step="1" min="0" max="100" value="${ST.fMissMax??''}" oninput="SNPCompare.setF('fMissMax',this.value)"></div>
        <div><div class="fl-lbl">Project</div>
          <select class="cmp-f" style="width:auto;min-width:230px;max-width:360px" onchange="SNPCompare.setF('fProj',this.value)">
            <option value="all">All projects</option>
            ${projOpts.map(p=>`<option value="${p.id}" ${ST.fProj===p.id?'selected':''}>${esc(p.label)}</option>`).join('')}
          </select></div>
        <button class="qbtn" onclick="SNPCompare.clearFilters()">Clear filters</button>
        <span style="margin-left:auto;display:flex;gap:8px">
          <button class="qbtn" onclick="SNPCompare.exportCSV()">Export CSV</button>
          ${hasLocal?`<button class="qbtn" onclick="SNPCompare.toTree()">Open region in SNPTree</button>`:''}
          ${hasLocal?`<button class="qbtn" onclick="SNPCompare.toMatrix()">Open in SNPMatrix</button>`:''}
        </span>
      </div>
    </div>

    <div id="cmpBanner"></div>
    <div id="cmpCount" style="color:var(--muted);font-size:13px;margin:0 0 8px"></div>
    <div class="card pad fade" id="cmpTableWrap" style="overflow:auto;max-height:70vh"></div>`;
  }
  function scopeBtn(m,label,enabled,reason){
    return `<button class="qbtn ${ST.mode===m?'solid':''}" ${enabled?'':`disabled title="${esc(reason||'Unavailable')}"`}
      style="${enabled?'':'opacity:.45;cursor:not-allowed'}" onclick="SNPCompare.setMode('${m}')">${label}</button>`;
  }

  /* ---------------- table ---------------- */
  function cols(){
    const c=[{k:'rank',t:'#'},{k:'id',t:'Final_ID'}];
    if(ST.mode==='both'){
      c.push({k:'gsim',t:'Global sim'},{k:'lsim',t:'Local sim'},{k:'dsim',t:'Δ (local−global)'},
             {k:'gmiss',t:'Global miss%'},{k:'lmiss',t:'Local miss%'});
    } else if(ST.mode==='local'){
      c.push({k:'lsim',t:'Similarity'},{k:'lmiss',t:'Missing%'},{k:'lboth',t:'Co-called sites'});
    } else {
      c.push({k:'gsim',t:'Similarity'},{k:'gmiss',t:'Missing%'});
    }
    c.push({k:'bio',t:'Project'},{k:'run',t:'SRA ID'},{k:'name',t:'Accession Name'});
    return c;
  }
  const fmtSim=v=>v==null?'—':v.toFixed(4);
  const fmtMiss=v=>v==null?'—':v.toFixed(2);
  const fmtD=v=>v==null?'—':(v>=0?'+':'')+v.toFixed(4);
  function cellVal(r,k){
    if(k==='gsim'||k==='lsim') return fmtSim(r[k]);
    if(k==='gmiss'||k==='lmiss') return fmtMiss(r[k]);
    if(k==='dsim') return fmtD(r.dsim);
    if(k==='lboth') return r.lboth==null?'—':r.lboth.toLocaleString();
    return esc(r[k]==null?'':r[k]);
  }
  function renderTable(){
    const rows=viewRows();
    const wrap=document.getElementById('cmpTableWrap'); if(!wrap) return;
    const cs=cols();
    const sortKey = ST.mode==='local' && ST.sortKey==='sim' ? 'lsim' : (ST.mode!=='local'&&ST.sortKey==='sim'?'gsim':ST.sortKey);
    const head=cs.map(c=>{
      if(c.k==='rank') return `<th>#</th>`;
      const active = (c.k===sortKey) || (c.k==='gsim'&&ST.sortKey==='sim'&&ST.mode!=='local') || (c.k==='lsim'&&ST.sortKey==='sim'&&ST.mode==='local');
      const arrow = active ? (ST.sortDir<0?' ▾':' ▴') : '';
      return `<th onclick="SNPCompare.sortBy('${c.k}')" style="cursor:pointer;white-space:nowrap">${esc(c.t)}${arrow}</th>`;
    }).join('');
    const body=rows.map((r,i)=>{
      const focal=r.id===ST.focal;
      const tds=cs.map(c=> c.k==='rank'
        ? `<td class="num">${i+1}</td>`
        : `<td class="${(c.k.endsWith('sim')||c.k.endsWith('miss')||c.k==='dsim'||c.k==='lboth')?'num':''}">${cellVal(r,c.k)}</td>`).join('');
      return `<tr class="${focal?'cmp-focal':''}">${tds}</tr>`;
    }).join('');
    wrap.innerHTML=`<table class="cmp-table"><thead><tr>${head}</tr></thead><tbody>${body}</tbody></table>`;
    const cnt=document.getElementById('cmpCount');
    if(cnt) cnt.innerHTML=`Showing <b>${rows.length.toLocaleString()}</b> of <b>${ST.allRows.length.toLocaleString()}</b> accessions for “${esc(ST.focal||'')}”.`;
    const ban=document.getElementById('cmpBanner');
    if(ban) ban.innerHTML = (ST.gdemo && ST.mode!=='local')
      ? `<div class="mtx-note">Genome-wide values shown here are <b>demonstration data</b>. Connect <span class="c-mono">ibsCompare.php</span> and set <span class="c-mono">useDemoGlobal=false</span> in snpcompare.js to load your precomputed IBS scores.</div>` : '';
  }

  /* ---------------- orchestration ---------------- */
  async function recompute(){
    const wrap=document.getElementById('cmpTableWrap');
    if(ST.mode!=='global' && !ST.input){ if(wrap) wrap.innerHTML=notice('Send a result from SNPVersity to enable region (local) comparison.'); return; }
    if(!ST.focal){ if(wrap) wrap.innerHTML=notice('Pick a focal accession to begin.'); return; }
    if(ST.mode!=='global' && ST.input && !ST.input.accs.some(a=>a.id===ST.focal)){
      if(wrap) wrap.innerHTML=notice(`Focal accession “${esc(ST.focal)}” isn’t in the region result. Pick one of the accessions you queried, or switch scope to Genome-wide.`); return;
    }
    if(wrap) wrap.innerHTML=`<div class="loading"><div class="spinner"></div><div>Computing similarities…</div></div>`;
    try{ await buildRows(); renderTable(); }
    catch(err){ if(wrap) wrap.innerHTML=notice('Could not load similarities: '+esc(err.message||err)); }
  }

  /* ---------------- public setters ---------------- */
  function setFocalFromInput(){ const el=document.getElementById('cmpFocal'); if(!el)return;
    const v=el.value.trim(); if(!v)return; ST.focal=v; recompute(); }
  function setMode(m){
    if((m==='local'||m==='both') && !ST.input) return;
    if((m==='global'||m==='both') && !globalAvailable()) return;
    ST.mode=m;
    // re-render shell to update scope highlight + columns baseline
    document.getElementById('page').innerHTML=shell(); recompute(); }
  function setF(k,v){ ST[k]= (v===''||v==null)?null:(k==='fProj'?v:parseFloat(v)); if(k==='fProj')ST.fProj=v; renderTable(); }
  function clearFilters(){ ST.fSimMin=ST.fSimMax=ST.fMissMax=null; ST.fProj='all';
    document.getElementById('page').innerHTML=shell(); recompute(); }
  function sortBy(k){ if(k==='rank')return;
    if(ST.sortKey===k){ ST.sortDir*=-1; } else { ST.sortKey=k; ST.sortDir=(k==='id'||k==='name'||k==='run'||k==='bio')?1:-1; }
    renderTable(); }
  function toTree(){ if(!ST.input)return; S.treeInput=ST.input; go('snptree'); }
  function toMatrix(){ if(!ST.input)return; S.matrixInput=ST.input; go('snpmatrix'); }

  function exportCSV(){
    const cs=cols(), rows=viewRows();
    const head=cs.map(c=>c.k==='rank'?'rank':c.t).join(',');
    const lines=rows.map((r,i)=>cs.map(c=>{
      if(c.k==='rank') return i+1;
      let v = c.k==='dsim'?r.dsim : r[c.k];
      if(v==null) return '';
      if(typeof v==='string' && /[",\n]/.test(v)) return '"'+v.replace(/"/g,'""')+'"';
      return v;
    }).join(','));
    const csv=[head].concat(lines).join('\n');
    const blob=new Blob([csv],{type:'text/csv'});
    const u=URL.createObjectURL(blob), a=document.createElement('a');
    a.href=u; a.download=`snpcompare_${ST.mode}_${(ST.focal||'focal').replace(/[^\w.-]/g,'_')}.csv`;
    document.body.appendChild(a); a.click(); a.remove(); setTimeout(()=>URL.revokeObjectURL(u),1500);
  }

  /* ---------------- helpers ---------------- */
  function esc(s){ return String(s==null?'':s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }
  function notice(html){ return `<div style="text-align:center;padding:26px;color:var(--muted);max-width:640px;margin:0 auto">${html}</div>`; }
  function injectCSS(){
    if(document.getElementById('snpcompare-css'))return;
    const s=document.createElement('style'); s.id='snpcompare-css';
    s.textContent=`
      .fl-lbl{font-size:10.5px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.4px;margin-bottom:7px}
      .cmp-f{border:1px solid var(--line);border-radius:9px;padding:8px 10px;font-size:13px;width:110px}
      .mtx-note{font-size:12px;color:var(--muted);background:var(--blue-50,#eef4ff);border:1px solid var(--line);border-radius:8px;padding:9px 12px;margin-bottom:12px}
      .mtx-note b{color:var(--ink)}
      .btn.solid{background:var(--green-50,#e9f4ec);border-color:#bfe5cb;color:var(--green-600,#1f8a4c)}
      .cmp-table{border-collapse:collapse;width:100%;font-size:13px}
      .cmp-table thead th{position:sticky;top:0;background:#eef1f6;border-bottom:1px solid var(--line);
        text-align:left;padding:9px 12px;font-weight:600;color:var(--ink);z-index:1}
      .cmp-table td{padding:8px 12px;border-bottom:1px solid #eef1f5;white-space:nowrap}
      .cmp-table td.num{text-align:right;font-family:var(--mono)}
      .cmp-table tbody tr:nth-child(even){background:#fafbfd}
      .cmp-table tbody tr:hover{background:#eef4ff}
      .cmp-table tr.cmp-focal{background:#fff6cf !important;font-weight:600}`;
    document.head.appendChild(s);
  }

  if(typeof SNPTools!=='undefined') SNPTools.register('snpcompare', { render });

  return { render, setFocalFromInput, setMode, setF, clearFilters, sortBy, toTree, toMatrix, exportCSV,
           // testing / debugging
           dose, localCompute, getGlobal, _CFG:CFG, _ST:ST };
})();
if(typeof window!=='undefined') window.SNPCompare = SNPCompare;
