/* =====================================================================
 *  snpimpact.js — AI-driven prioritization & interpretation (real data).
 *
 *  Ranks, sorts, and filters the variants in a region (regardless of
 *  accessions), using the same VCF SNPVersity produced: DNA + protein
 *  language-model scores (PlantCAD / ESM), predicted consequence, and
 *  real Pfam domain annotation. Receives the region via S.impactInput.
 * ===================================================================== */
(function () {

  const IMP = {
    rows: [], input: null, _sig: null,
    sortKey: 'priority', sortDir: 1,        // 1 = TOP first
    fCons: 'all', fImpact: 'all', fScore: 'all', fDomain: 'all',
    openId: null,
    shortlist: new Set(),
  };
  const CAP = 1500;                          // max rows rendered at once
  const PRIO_RANK = { TOP:0, HIGH:1, MODERATE:2, LOW:3 };

  /* ---------- small helpers ---------- */
  function consPill(r){ return `<span class="cons ${r.consClass}">${r.consequence}</span>`; }
  function prioPill(p){ return `<span class="prio ${(p||'LOW').toLowerCase()}">${p||'LOW'}</span>`; }
  function scoreColor(v){
    const lo=-12, hi=6, t=Math.max(0,Math.min(1,(v-lo)/(hi-lo)));
    const r = t<.5 ? 255 : Math.round(255*(1-(t-.5)*2));
    const g = t<.5 ? Math.round(255*t*2) : 200;
    return `rgb(${r},${Math.max(60,g)},60)`;
  }
  function scorePct(v){ const lo=-12, hi=6; return Math.max(0,Math.min(100,((v-lo)/(hi-lo))*100)); }
  function scoreCell(v){
    if (v==null) return '<span style="color:var(--faint)">—</span>';
    return `<span class="imp-score" style="background:${scoreColor(v)}">${v>0?'+':''}${v.toFixed(1)}</span>`;
  }
  const fmtMaf = m => (m==null ? '—' : (+m).toFixed(3));

  /* ---------- PanEffect jump ---------- */
  /* true when this row is a single amino-acid substitution */
  function isMissense(r){
    return r.consClass==='missense' || /missense/i.test(r.consequence||'');
  }
  /* substitution string PanEffect highlights on, e.g. "A123T".
     SNPVersity rows carry it as r.sub; SNPImpact rows carry the parts. */
  function aaSub(r){
    if (r.sub) return r.sub;
    if (r.aaRef && r.aaAlt && r.resi!=null) return `${r.aaRef}${r.resi}${r.aaAlt}`;
    return '';
  }
  function canPanEffect(r){ return !!(r && r.gene && r.gene!=='—'); }
  /* internal view switch — highlights the substitution when missense,
     otherwise just opens PanEffect on the gene */
  function panEffect(r){
    if (!canPanEffect(r)) return;
    const sub = isMissense(r) ? aaSub(r) : '';
    if (typeof goPanEffect === 'function') goPanEffect(r.gene, sub ? {variant:sub} : undefined);
    else go('paneffect');
  }

  function filtered(){
    let r = IMP.rows.filter(v =>
      (IMP.fCons==='all'   || v.consClass===IMP.fCons) &&
      (IMP.fImpact==='all' || v.priority===IMP.fImpact) &&
      (IMP.fScore==='all'  || (v.combined!=null && v.combined<=-4)) &&
      (IMP.fDomain==='all' || (IMP.fDomain==='dom' ? v.domain!=='—' : v.domain==='—'))
    );
    const k = IMP.sortKey, dir = IMP.sortDir;
    r.sort((a,b)=>{
      if (k==='priority'){
        const d = PRIO_RANK[a.priority]-PRIO_RANK[b.priority];
        if (d) return dir*d;
        const ac=a.combined==null?Infinity:a.combined, bc=b.combined==null?Infinity:b.combined;
        return ac-bc;                                    // tie-break: most deleterious first
      }
      if (k==='gene'||k==='consequence'||k==='domain')
        return dir*String(a[k]).localeCompare(String(b[k]));
      const av=a[k], bv=b[k], na=av==null, nb=bv==null;   // numeric; nulls last
      if (na&&nb) return 0; if (na) return 1; if (nb) return -1;
      return dir*(av-bv);
    });
    return r;
  }

  /* ---------- render ---------- */
  function render(page){
    page = page || document.getElementById('page');
    const input = S.impactInput;
    if (!input || !input.rows || !input.rows.length){ page.innerHTML = emptyState(); return; }

    const sig = `${input.chr}:${input.start}-${input.end}:${input.dataset}:${input.rows.length}`;
    if (IMP._sig !== sig){
      IMP.rows = Data.rankImpact(input.rows);
      IMP._sig = sig; IMP.input = input; IMP.openId = null; IMP.shortlist.clear();
    }
    IMP.sec = Data.hasSecondaryScores(input.dataset);    // show PlantCAD2/ESM2 only for MaizeGDB 2026
    Data.ensureGeneDomains();                             // warm up detail track (non-blocking)
    Data.ensureGeneModels(input.chr);                     // warm up gene-model view for this chromosome

    const region = `${input.chr}:${(+input.start).toLocaleString()}–${(+input.end).toLocaleString()}`;
    const all = filtered();
    const rows = all.slice(0, CAP);
    const top = IMP.rows.filter(r=>r.priority==='TOP'||r.priority==='HIGH').length;

    page.innerHTML = `
      <div class="sec"><div class="bar"></div><div>
        <div class="n">AI-DRIVEN PRIORITIZATION · PlantCAD + ESM</div>
        <h2>Rank candidate variants using AI</h2>
        <p>Combine DNA and protein language-model scores with predicted effects to surface
        the variants most likely to be causal. Sort, filter, and open any variant to see its
        gene-model consequence, domain impact, and AI score summary.</p>
      </div></div>

      <div class="imp-context">
        <span><b>${IMP.rows.length.toLocaleString()}</b> variants in region</span>
        <span class="dot">·</span><span><b>${top.toLocaleString()}</b> high-priority</span>
        <span class="dot">·</span><span>region <span class="mono">${region}</span></span>
        <span class="dot">·</span><span><span class="mono">${esc(input.datasetName||input.dataset)}</span></span>
      </div>

      <div class="filterbar imp-filters">
        ${sel('Category','fCons',[['all','All consequences'],['lof','Loss-of-function'],['splice','Splice'],['missense','Missense'],['indel','In-frame indel'],['syn','Synonymous'],['other','Non-coding / other']])}
        ${sel('Priority','fImpact',[['all','All priorities'],['TOP','TOP'],['HIGH','HIGH'],['MODERATE','MODERATE'],['LOW','LOW']])}
        ${sel('AI score','fScore',[['all','All scores'],['high','AI high-priority (≤ −4)']])}
        ${sel('Domain effect','fDomain',[['all','All'],['dom','In a Pfam domain'],['nodom','No domain hit']])}
        <div class="right">
          <button class="btn" onclick="IMPACT.exportCSV()">${ICONS.download||''} Export CSV</button>
          <button class="btn" onclick="IMPACT.sendCompare()">
            <svg viewBox="0 0 24 24" fill="none"><path d="M4 12h13m0 0l-5-5m5 5l-5 5" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"/></svg>
            Send region to SNPCompare
          </button>
        </div>
      </div>

      <div class="tbl-wrap" style="max-height:none">
        <table class="vcf imp">
          <thead><tr>
            ${th('Gene','gene')}<th>Variant</th>${th('Consequence','consequence')}${th('Domain','domain')}
            ${th('PlantCAD1','plantcad','num')}${IMP.sec?th('PlantCAD2','plantcad2','num'):''}${th('ESM','esm','num')}${IMP.sec?th('ESM2','esm2','num')+th('ESM3','esm3','num'):''}${th('Priority','priority')}<th></th>
          </tr></thead>
          <tbody>${rows.map(rowHTML).join('')}</tbody>
        </table>
      </div>
      <div class="pager"><div class="info">Showing <span class="mono">${rows.length.toLocaleString()}</span> of
        <span class="mono">${all.length.toLocaleString()}</span> variants${all.length>CAP?` (top ${CAP} — refine filters to see more)`:''}</div></div>

      <div id="impDetail">${IMP.openId ? detailHTML(IMP.rows.find(r=>r.id===IMP.openId)) : ''}</div>
    `;
    if (typeof attachTT==='function') attachTT();
  }

  function emptyState(){
    return `<div class="empty-state"><div class="ei">${ICONS.star||''}</div>
      <h3>Send a region from SNPVersity</h3>
      <p>SNPImpact ranks the variants in a queried region by predicted impact — independent of
      which accessions you picked. Run a query in SNPVersity, then use
      <b>Send to SNPImpact</b> (or the “Send selection to…” menu).</p>
      <button class="btn solid" onclick="go('snpversity')">Go to SNPVersity</button></div>`;
  }

  function sel(label, key, opts){
    return `<div class="fld"><label>${label}</label>
      <select onchange="IMPACT.setFilter('${key}', this.value)">
        ${opts.map(o=>`<option value="${o[0]}" ${IMP[key]===o[0]?'selected':''}>${o[1]}</option>`).join('')}
      </select></div>`;
  }
  function th(label, key, cls){
    const active = IMP.sortKey===key;
    const arrow = active ? (IMP.sortDir>0?' ▲':' ▼') : ' ⇅';
    return `<th class="sortable ${cls||''} ${active?'on':''}" onclick="IMPACT.sort('${key}')">${label}<span class="arr">${arrow}</span></th>`;
  }

  function rowHTML(r){
    const star = IMP.shortlist.has(r.id);
    const opened = IMP.openId===r.id;
    const sub = isMissense(r) ? aaSub(r) : '';
    const peJump = (isMissense(r) && canPanEffect(r))
      ? ` <a class="pe-jump" href="#" title="View ${esc(sub||'this substitution')} in PanEffect"
             onclick="event.stopPropagation();IMPACT.panEffect('${r.id}');return false;">effects ↗</a>`
      : '';
    return `<tr class="imp-row ${opened?'open':''}" onclick="IMPACT.open('${r.id}')">
      <td class="gene-link" style="padding-left:11px">${r.gene}</td>
      <td class="c-mono c-alt">${r.variant}</td>
      <td>${consPill(r)}${peJump}</td>
      <td>${domTag(r.domain)}</td>
      <td class="num">${scoreCell(r.plantcad)}</td>
      ${IMP.sec?`<td class="num">${scoreCell(r.plantcad2)}</td>`:''}
      <td class="num">${scoreCell(r.esm)}</td>
      ${IMP.sec?`<td class="num">${scoreCell(r.esm2)}</td><td class="num">${scoreCell(r.esm3)}</td>`:''}
      <td>${prioPill(r.priority)}</td>
      <td style="text-align:center">
        <button class="star-btn ${star?'on':''}" title="Add to shortlist"
          onclick="event.stopPropagation();IMPACT.star('${r.id}')">${star?'★':'☆'}</button>
      </td>
    </tr>`;
  }

  /* ---------- gene model: semantic zoom (affected exon + flanks + UTR + 1kb up/down) ---------- */
  function gmLocate(model, pos){
    const ex = model.exons, N = ex.length, strand = model.strand==='-'?'-':'+';
    const exonNo = i => strand==='+' ? i+1 : N-i;
    for (let i=0;i<N;i++) if (pos>=ex[i][0] && pos<=ex[i][1]) return {kind:'exon', i, N, strand, exonNo:exonNo(i)};
    if (pos < model.start){ return {kind: strand==='+'?'upstream':'downstream', i:0, N, strand}; }
    if (pos > model.end){   return {kind: strand==='+'?'downstream':'upstream', i:N-1, N, strand}; }
    let ia=0; for (let i=0;i<N-1;i++){ if (pos>ex[i][1] && pos<ex[i+1][0]){ ia=i; break; } }
    return {kind:'intron', i:ia, N, strand};
  }
  function gmWindow(loc){
    const N=loc.N; let idx;
    if (loc.kind==='exon')        idx=[loc.i-1, loc.i, loc.i+1];
    else if (loc.kind==='intron') idx=[loc.i, loc.i+1];
    else                          idx=(loc.i===0)?[0,1]:[N-2,N-1];
    idx = [...new Set(idx.filter(i=>i>=0 && i<N))].sort((a,b)=>a-b);
    if (idx.length>3){ idx = idx.sort((a,b)=>Math.abs(a-loc.i)-Math.abs(b-loc.i)).slice(0,3).sort((a,b)=>a-b); }
    return idx;
  }
  function geneModelSVG(r, model){
    if (!model || !model.exons || !model.exons.length) return null;
    const ex=model.exons, cds=model.cds||[], N=ex.length, strand=model.strand==='-'?'-':'+';
    const loc=gmLocate(model, r.pos), idx=gmWindow(loc);
    const cdsSpan=(s,e)=>{ let cs=null,ce=null; for(const [a,b] of cds){ const lo=Math.max(s,a),hi=Math.min(e,b); if(lo<=hi){ cs=cs==null?lo:Math.min(cs,lo); ce=ce==null?hi:Math.max(ce,hi);} } return cs==null?null:[cs,ce]; };
    const exonW=L=>Math.max(36,Math.min(130,Math.round(30+L*0.12)));
    const INTRON=28, FLANK=48, ELL=20, PAD=18, yMid=42, cdsH=20, utrH=8;
    const leftFlank=idx[0]===0, rightFlank=idx[idx.length-1]===N-1, leftEll=idx[0]>0, rightEll=idx[idx.length-1]<N-1;
    let segs=[], x=PAD;
    if (leftFlank){ segs.push({type:'flank',side:'left',x,w:FLANK,label:strand==='+'?'upstream':'downstream'}); x+=FLANK; }
    else if (leftEll){ segs.push({type:'ell',x,w:ELL}); x+=ELL; }
    for (let j=0;j<idx.length;j++){
      const i=idx[j], s=ex[i][0], e=ex[i][1], w=exonW(e-s+1);
      segs.push({type:'exon',i,s,e,x,w,cds:cdsSpan(s,e),exonNo:strand==='+'?i+1:N-i}); x+=w;
      if (j<idx.length-1){ segs.push({type:'intron',x,w:INTRON}); x+=INTRON; }
    }
    if (rightFlank){ segs.push({type:'flank',side:'right',x,w:FLANK,label:strand==='+'?'downstream':'upstream'}); x+=FLANK; }
    else if (rightEll){ segs.push({type:'ell',x,w:ELL}); x+=ELL; }
    const W=x+PAD, H=70;
    const lab=(tx,ty,t,anchor='middle',col='#6b7c98',sz=9.5)=>`<text x="${tx}" y="${ty}" text-anchor="${anchor}" font-size="${sz}" fill="${col}" font-family="Inter,system-ui,sans-serif">${t}</text>`;
    let body=lab(6,yMid+4,strand==='+'?'5′':'3′','start','#33456a',11)+lab(W-6,yMid+4,strand==='+'?'3′':'5′','end','#33456a',11);
    for (const sg of segs){
      if (sg.type==='intron'){
        body+=`<line x1="${sg.x}" y1="${yMid}" x2="${sg.x+sg.w}" y2="${yMid}" stroke="#9fb0c9" stroke-width="1.5"/>`;
        const cx=sg.x+sg.w/2, d=strand==='+'?1:-1;
        body+=`<path d="M${cx-3*d},${yMid-3} L${cx+3*d},${yMid} L${cx-3*d},${yMid+3}" fill="none" stroke="#9fb0c9" stroke-width="1.3"/>`;
      } else if (sg.type==='flank'){
        body+=`<line x1="${sg.x}" y1="${yMid}" x2="${sg.x+sg.w}" y2="${yMid}" stroke="#c7d2e0" stroke-width="1.5" stroke-dasharray="3 2"/>`+lab(sg.x+sg.w/2,yMid+17,'±1kb');
      } else if (sg.type==='ell'){
        body+=lab(sg.x+sg.w/2,yMid+4,'· · ·','middle','#9fb0c9',13);
      } else {
        body+=`<rect x="${sg.x}" y="${yMid-utrH/2}" width="${sg.w}" height="${utrH}" rx="2" fill="#d3deee"/>`;
        if (sg.cds){ const f0=(sg.cds[0]-sg.s)/(sg.e-sg.s), f1=(sg.cds[1]-sg.s)/(sg.e-sg.s);
          body+=`<rect x="${sg.x+f0*sg.w}" y="${yMid-cdsH/2}" width="${Math.max(2,(f1-f0)*sg.w)}" height="${cdsH}" rx="2.5" fill="#2f6ad0"/>`; }
        body+=lab(sg.x+sg.w/2,yMid+17,'Exon '+sg.exonNo);
      }
    }
    // variant marker + caption
    let mx=null, caption='';
    const aff=segs.find(sg=>sg.type==='exon' && r.pos>=sg.s && r.pos<=sg.e);
    if (aff){ mx=aff.x+((r.pos-aff.s)/(aff.e-aff.s))*aff.w;
      if (aff.cds && r.pos>=aff.cds[0] && r.pos<=aff.cds[1]) caption=`Exon ${aff.exonNo} of ${N} (coding)`;
      else { const before=r.pos<(aff.cds?aff.cds[0]:aff.e); caption=((before===(strand==='+'))?'5′':'3′')+' UTR · exon '+aff.exonNo; }
    } else if (loc.kind==='intron'){ const iv=segs.find(sg=>sg.type==='intron'); if(iv){ mx=iv.x+iv.w/2; caption='Intron (between exons)'; } }
    else { const fl=segs.find(sg=>sg.type==='flank' && sg.label===loc.kind); if(fl){ mx=fl.side==='left'?fl.x+7:fl.x+fl.w-7; } caption=(loc.kind==='upstream'?'≤1 kb upstream':'≤1 kb downstream')+' of gene'; }
    if (mx!=null) body+=`<path d="M${mx-4},${yMid-19} l8,0 l-4,7 z" fill="#c0362c"/><line x1="${mx}" y1="${yMid-13}" x2="${mx}" y2="${yMid+cdsH/2}" stroke="#c0362c" stroke-width="1.6"/>`;
    return { svg:`<svg viewBox="0 0 ${W} ${H}" width="100%" style="max-width:${Math.max(W,320)}px;height:auto;overflow:visible">${body}</svg>`, caption };
  }

  /* ---------- variant detail: consequence → real domain track → AI scores ---------- */
  function consequenceBlurb(r){
    if (r.consClass==='lof')     return 'Loss-of-function: '+r.consequence.toLowerCase()+(r.resi!=null?(' at residue '+r.resi):'')+' — likely disrupts the protein.';
    if (r.consClass==='splice')  return 'Splice-site change — may disrupt normal transcript splicing.';
    if (r.consClass==='missense')return 'Missense'+(r.aaRef&&r.aaAlt?(' ('+r.aaRef+r.resi+r.aaAlt+')'):'')+' — single amino-acid substitution.';
    if (r.consClass==='indel')   return r.consequence+' — in-frame change to protein length.';
    if (r.consClass==='syn')     return 'Synonymous — no amino-acid change.';
    return r.consequence+' — non-coding / regulatory region.';
  }
  function geneDomainTrack(r){
    const gd = Data.geneDomains(r.gene);
    if (!gd || !gd.len){
      return `<div class="id-sub">${r.domain==='—'?'No annotated Pfam domain at this position.':'Falls in '+r.domain+'.'}</div>`;
    }
    const L = gd.len;
    const segs = (gd.domains||[]).map(d=>{
      const left = 100*(d.start-1)/L, w = Math.max(1.5, 100*(d.end-d.start+1)/L);
      const hit = r.resi!=null && r.resi>=d.start && r.resi<=d.end;
      return `<div class="dom-seg ${hit?'':'alt'}" style="left:${left.toFixed(1)}%;width:${w.toFixed(1)}%"
                   title="${esc(d.name)} (${d.pfam}) ${d.start}-${d.end}">${esc(d.name)}</div>`;
    }).join('');
    const marker = (r.resi!=null && r.resi<=L) ? `<div class="trunc" style="left:${(100*(r.resi-1)/L).toFixed(1)}%"></div>` : '';
    const inDom = r.resi!=null ? (gd.domains||[]).find(d=>r.resi>=d.start&&r.resi<=d.end) : null;
    return `
      <div class="id-sub">${inDom
          ? `Residue ${r.resi} lies in the <b>${esc(inDom.name)}</b> domain (${inDom.pfam}).`
          : r.resi!=null ? `Residue ${r.resi} is outside annotated domains.` : 'Non-coding variant — not placed on the protein.'}</div>
      <div class="protbar">
        <span class="pp" style="left:0">1</span>${segs}${marker}<span class="pp" style="right:0">${L}</span>
      </div>`;
  }
  function detailHTML(r){
    if (!r) return '';
    const peSub = isMissense(r) ? aaSub(r) : '';
    return `
    <div class="imp-detail card fade">
      <div class="idh">
        <div>
          <div class="idh-t">Variant detail</div>
          <div class="idh-meta">
            <span>Gene: <b class="mono">${r.gene}</b></span>
            <span>Variant: <b class="mono">${r.variant}</b></span>
            <span>${consPill(r)}</span>
            <span>${prioPill(r.priority)}</span>
          </div>
        </div>
        <div class="idh-actions">
          <button class="btn" onclick="goFold('${r.gene}')">${ICONS.fold} View in SNPFold</button>
          <button class="btn" onclick="goFunction('${r.gene}','${IMP.input.dataset}')">${ICONS.leaf||''} Gene in SNPFunction</button>
          <button class="btn ghost" onclick="IMPACT.panEffect('${r.id}')"
            title="${peSub ? 'Highlight '+esc(peSub)+' in PanEffect' : 'View '+esc(r.gene)+' in PanEffect'}"
            ${canPanEffect(r)?'':'disabled'}>${ICONS.effect||''} ${peSub ? `PanEffect · ${esc(peSub)}` : 'PanEffect'} →</button>
          <button class="btn ghost" onclick="IMPACT.close()" aria-label="close">✕</button>
        </div>
      </div>

      <div class="id-steps">
        <div class="id-step">
          <div class="id-step-h"><span class="num-dot">1</span> Predicted consequence <span class="muted">· gene model</span></div>
          <div class="id-sub">${consequenceBlurb(r)}</div>
          ${(()=>{ const gm=geneModelSVG(r, Data.geneModelOf(IMP.input.chr, r.gene));
            return gm
              ? `<div class="gm-wrap" style="margin:8px 0 2px">${gm.svg}</div>
                 <div class="id-foot ${r.consClass==='lof'?'red':''}">${gm.caption} · ${IMP.input.chr}:${r.pos.toLocaleString()} ${r.ref}›${r.alt} · impact ${r.impactLevel||'—'}</div>`
              : `<div class="id-foot ${r.consClass==='lof'?'red':''}">${IMP.input.chr}:${r.pos.toLocaleString()} · ${r.ref}›${r.alt} · impact ${r.impactLevel||'—'}</div>`;
          })()}
        </div>

        <div class="id-step">
          <div class="id-step-h"><span class="num-dot">2</span> Protein / domain view <span class="muted">· Pfam</span></div>
          ${geneDomainTrack(r)}
        </div>

        <div class="id-step">
          <div class="id-step-h"><span class="num-dot">3</span> AI score summary</div>
          <div class="ai-pill ${r.combined!=null&&r.combined<=-4?'hi':''}">AI score: ${r.combined==null?'n/a':r.combined<=-7?'high impact':r.combined<=-4?'elevated':'low impact'}</div>
          ${scoreBar('PlantCAD1', r.plantcad)}
          ${IMP.sec?scoreBar('PlantCAD2', r.plantcad2):''}
          ${scoreBar('ESM', r.esm)}
          ${IMP.sec?scoreBar('ESM2', r.esm2)+scoreBar('ESM3', r.esm3):''}
          <div class="pctl">
            <div class="pctl-l">Region impact percentile</div>
            <div class="pctl-bar"><div class="pctl-fill" style="width:${r.percentile==null?0:r.percentile}%"></div><span class="pctl-v">${r.percentile==null?'n/a':r.percentile+'th'}</span></div>
          </div>
          <div class="muted" style="font-size:11px;margin-top:6px">MAF ${fmtMaf(r.maf)} · imputation r² ${r.r2==null?'—':(+r.r2).toFixed(2)} · MQ ${r.mq==null?'—':r.mq}</div>
        </div>
      </div>

      <div class="id-actions">
        <button class="btn" onclick="IMPACT.close()">${ICONS.compare} Back to list</button>
        <button class="btn" onclick="IMPACT.star('${r.id}')">${IMP.shortlist.has(r.id)?'★ On shortlist':'☆ Add to shortlist'}</button>
      </div>
    </div>`;
  }
  function scoreBar(label, v){
    if (v==null) return `<div class="sbar"><div class="sbar-h"><span>${label} score</span><b class="mono">—</b></div><div class="sbar-track" style="background:linear-gradient(90deg,#e5484d 0%,#f5b545 50%,#2f9e44 100%)"></div></div>`;
    return `<div class="sbar">
      <div class="sbar-h"><span>${label} score</span><b class="mono">${v>0?'+':''}${v.toFixed(1)}</b></div>
      <div class="sbar-track" style="background:linear-gradient(90deg,#e5484d 0%,#f5b545 50%,#2f9e44 100%)"><span class="sbar-mark" style="left:${scorePct(v)}%"></span></div>
    </div>`;
  }

  function esc(s){ return String(s==null?'':s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }

  /* ---------- public handlers ---------- */
  window.IMPACT = {
    setFilter(k,v){ IMP[k]=v; render(); },
    sort(k){ if(IMP.sortKey===k) IMP.sortDir*=-1;
      else { IMP.sortKey=k; IMP.sortDir = (k==='gene'||k==='consequence'||k==='domain')?1:1; } render(); },
    open(id){ IMP.openId = (IMP.openId===id?null:id);
      Promise.all([Data.ensureGeneDomains(), Data.ensureGeneModels(IMP.input.chr)]).then(()=>{ render();
        const d=document.getElementById('impDetail'); if(IMP.openId && d && d.scrollIntoView) d.scrollIntoView({behavior:'smooth',block:'nearest'}); }); },
    close(){ IMP.openId=null; render(); },
    panEffect(id){ panEffect(IMP.rows.find(r=>r.id===id)); },
    star(id){ IMP.shortlist.has(id)?IMP.shortlist.delete(id):IMP.shortlist.add(id); render(); },
    sendCompare(){ if(IMP.input){ S.compareInput=IMP.input; } go('snpcompare'); },
    exportCSV(){
      const cols=['gene','variant','consequence','domain','plantcad']
        .concat(IMP.sec?['plantcad2']:[])
        .concat(['esm'])
        .concat(IMP.sec?['esm2','esm3']:[])
        .concat(['combined','priority','percentile','impactLevel','maf','pos','ref','alt']);
      const rows=filtered();
      const line=r=>cols.map(c=>{ let v=r[c]; if(v==null)return '';
        if(typeof v==='string' && /[",\n]/.test(v)) return '"'+v.replace(/"/g,'""')+'"'; return v; }).join(',');
      const csv=[cols.join(',')].concat(rows.map(line)).join('\n');
      const blob=new Blob([csv],{type:'text/csv'}), u=URL.createObjectURL(blob), a=document.createElement('a');
      a.href=u; a.download=`snpimpact_${(IMP.input.chr||'region')}_${IMP.input.start}_${IMP.input.end}.csv`;
      document.body.appendChild(a); a.click(); a.remove(); setTimeout(()=>URL.revokeObjectURL(u),1500);
    },
  };

  SNPTools.register('snpimpact', { render });
})();
