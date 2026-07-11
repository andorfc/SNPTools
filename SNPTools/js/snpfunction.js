/* =====================================================================
 *  snpfunction.js — Gene function & allele mining.
 *
 *  Gene-scoped view (independent of a variant region): a functional
 *  dossier (Pfam domains, size, links), the gene's variant burden across
 *  the WHOLE panel, and a damaging / knockout allele catalog listing which
 *  accessions carry each damaging allele. Pulls Data.geneFunction(gene).
 * ===================================================================== */
(function () {

  const FN = { gene:null, dataset:null, data:null, loading:false, openId:null };

  function esc(s){ return String(s==null?'':s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }
  function scoreColor(v){ const lo=-12,hi=6,t=Math.max(0,Math.min(1,(v-lo)/(hi-lo))); const r=t<.5?255:Math.round(255*(1-(t-.5)*2)); const g=t<.5?Math.round(255*t*2):200; return `rgb(${r},${Math.max(60,g)},60)`; }
  function scoreCell(v){ return v==null?'<span style="color:var(--faint)">—</span>':`<span class="imp-score" style="background:${scoreColor(v)}">${v>0?'+':''}${v.toFixed(1)}</span>`; }
  function prioPill(p){ return `<span class="prio ${(p||'LOW').toLowerCase()}">${p||'LOW'}</span>`; }
  function consPill(v){ return `<span class="cons ${v.consClass}">${v.consequence}</span>`; }

  /* ---------- render ---------- */
  function render(page){
    page = page || document.getElementById('page');
    injectCSS();
    if (S.functionGene && S.functionGene !== FN.gene){ FN.gene = S.functionGene; FN.dataset = S.functionDataset || FN.dataset; FN.data = null; FN.openId = null; }
    S.functionGene = null;
    if (!FN.gene){ page.innerHTML = searchBar() + emptyState(); return; }

    if (!FN.data && !FN.loading){
      FN.loading = true;
      page.innerHTML = searchBar() + `<div class="loading" style="padding:44px;text-align:center"><div class="spinner"></div><div>Analyzing <b>${esc(FN.gene)}</b> across the panel…</div></div>`;
      Data.geneFunction(FN.gene, FN.dataset).then(d => { FN.data = d; FN.loading = false; render(); })
        .catch(e => { FN.loading = false; FN.data = {gene:FN.gene, error:(e&&e.message)||'failed'}; render(); });
      return;
    }
    const d = FN.data;
    if (d && d.error){ page.innerHTML = searchBar() + notice(`Couldn’t analyze “${esc(FN.gene)}”: ${esc(d.error)}`); return; }

    page.innerHTML = searchBar() + hero(d) + dossier(d) + burden(d) + catalog(d);
    if (typeof attachTT==='function') attachTT();
  }

  function searchBar(){
    return `<div class="card pad" style="margin-bottom:16px;display:flex;gap:10px;align-items:center;flex-wrap:wrap">
      <span style="font-size:10.5px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.4px">Gene model</span>
      <input id="fnGeneInput" value="${esc(FN.gene||'')}" placeholder="e.g. Zm00001eb378140" spellcheck="false"
        style="flex:1;min-width:280px;border:1px solid var(--line);border-radius:9px;padding:9px 11px;font-family:var(--mono);font-size:13px"
        onkeydown="if(event.key==='Enter')FUNCTION.load()">
      <button class="btn" onclick="FUNCTION.load()">Analyze gene</button>
    </div>`;
  }
  function emptyState(){
    return `<div class="empty-state"><div class="ei">${ICONS.leaf||ICONS.star||''}</div>
      <h3>Pick a gene to mine its functional variation</h3>
      <p>SNPFunction summarizes a gene across the whole panel: its Pfam domains, the burden of
      coding variation, and a catalog of damaging / knockout alleles with the accessions that carry them.
      Search a gene model above, or jump here from a gene in SNPVersity or SNPImpact.</p></div>`;
  }
  function notice(html){ return `<div style="text-align:center;padding:26px;color:var(--muted);max-width:640px;margin:0 auto">${html}</div>`; }

  /* ---------- hero + dossier ---------- */
  function hero(d){
    return `<div class="sec"><div class="bar"></div><div>
      <div class="n">GENE FUNCTION &amp; ALLELE MINING</div>
      <h2 style="margin-bottom:2px">${esc(d.gene)}</h2>
      <p>Functional variation for this gene across <b>${d.nAccessions.toLocaleString()}</b> ${esc(d.datasetName)} accessions —
      domains, variant burden, and which lines carry damaging or knocked-out copies.</p>
    </div></div>`;
  }
  function dossier(d){
    const region = `${d.chr}:${(+d.start).toLocaleString()}–${(+d.end).toLocaleString()}${d.strand?` (${d.strand})`:''}`;
    const jb = `https://jbrowse.maizegdb.org/index.html?data=B73&loc=${encodeURIComponent(d.gene)}`;
    const pe = `https://maizegdb.org/effect/maize_v2/index.html?id=${encodeURIComponent(d.gene)}`;
    const mg = `https://maizegdb.org/gene_center/gene/${encodeURIComponent(d.gene)}`;
    const doms = (d.domains||[]).length
      ? d.domains.map(x=>domTag(`${x.name} (${x.pfam})`)).join(' ')
      : '<span class="muted">No Pfam domains annotated</span>';
    return `<div class="card pad" style="margin-bottom:16px">
      <div class="fn-grid">
        <div><div class="fn-k">Location</div><div class="fn-v mono">${region}</div></div>
        <div><div class="fn-k">Protein</div><div class="fn-v">${d.protLen?`${d.protLen} aa`:'—'}${d.protein?` <span class="muted mono">${esc(d.protein)}</span>`:''}</div></div>
        <div><div class="fn-k">Variants in gene</div><div class="fn-v"><b>${d.nVariants.toLocaleString()}</b></div></div>
        <div style="flex:1 1 100%"><div class="fn-k">Pfam domains</div><div class="fn-v">${doms}</div></div>
      </div>
      <div style="display:flex;gap:8px;flex-wrap:wrap;margin-top:12px">
        <button class="btn" onclick="goFold('${esc(d.gene)}')">${ICONS.fold} View structure (SNPFold)</button>
        <a class="btn ghost" href="${pe}" target="_blank" rel="noopener">PanEffect ↗</a>
        <a class="btn ghost" href="${jb}" target="_blank" rel="noopener">JBrowse ↗</a>
        <a class="btn ghost" href="${mg}" target="_blank" rel="noopener">MaizeGDB ↗</a>
      </div>
    </div>`;
  }

  /* ---------- burden ---------- */
  function burden(d){
    const b = d.burden, bc = b.byClass;
    const sec = Data.hasSecondaryScores(d.dataset);
    const total = d.nVariants || 1;
    const barSeg = (n,cls,label)=> n? `<span class="fn-seg ${cls}" style="width:${Math.max(2,100*n/total)}%" title="${label}: ${n}"></span>`:'';
    const nss = b.nonsynSyn==null ? '∞' : b.nonsynSyn;
    const af = b.afSpectrum;
    return `<div class="card pad" style="margin-bottom:16px">
      <div class="fn-h">Variant burden</div>
      <div class="fn-stats">
        ${stat('Non-syn : syn', nss, 'coding constraint (higher = more nonsynonymous)')}
        ${stat('Exon : intron', b.exonIntron==null?'∞':b.exonIntron, 'variants in exons vs introns of the gene model')}
        ${stat('Domain-disrupting', b.domainDisrupting, 'coding variants inside a Pfam domain')}
        ${stat('Candidate KO lines', d.koLines, 'accessions homozygous for a loss-of-function allele')}
        ${stat('Mean PlantCAD1', b.meanPlantcad==null?'—':b.meanPlantcad, 'average PlantCAD1 DNA language-model score')}
        ${sec?stat('Mean PlantCAD2', b.meanPlantcad2==null?'—':b.meanPlantcad2, 'average PlantCAD2 (2026) DNA language-model score'):''}
        ${stat('Mean ESM', b.meanEsm==null?'—':b.meanEsm, 'average ESM protein language-model score')}
        ${sec?stat('Mean ESM2', b.meanEsm2==null?'—':b.meanEsm2, 'average ESM2 (2026) protein language-model score'):''}
        ${sec?stat('Mean ESM3', b.meanEsm3==null?'—':b.meanEsm3, 'average ESM3 (2026) protein language-model score'):''}
      </div>
      <div class="fn-barwrap">
        <div class="fn-bar">
          ${barSeg(bc.lof,'lof','Loss-of-function')}${barSeg(bc.splice,'splice','Splice')}${barSeg(bc.missense,'missense','Missense')}${barSeg(bc.indel,'indel','In-frame indel')}${barSeg(bc.syn,'syn','Synonymous')}${barSeg(bc.other,'other','Non-coding / other')}
        </div>
        <div class="fn-legend">
          ${leg('lof','LOF',bc.lof)} ${leg('splice','Splice',bc.splice)} ${leg('missense','Missense',bc.missense)} ${leg('indel','Indel',bc.indel)} ${leg('syn','Synon.',bc.syn)} ${leg('other','Other',bc.other)}
        </div>
      </div>
      <div class="fn-af">Allele frequency: <b>${af.rare}</b> rare (&lt;1%) · <b>${af.low}</b> low (1–5%) · <b>${af.common}</b> common (≥5%)</div>
    </div>`;
  }
  function stat(k,v,tip){ return `<div class="fn-stat" title="${esc(tip)}"><div class="fn-statv">${v}</div><div class="fn-statk">${k}</div></div>`; }
  function leg(cls,label,n){ return `<span class="fn-lg"><span class="fn-sw ${cls}"></span>${label} ${n}</span>`; }

  /* ---------- damaging / knockout catalog ---------- */
  function catalog(d){
    if (!d.damaging.length)
      return `<div class="card pad"><div class="fn-h">Damaging &amp; knockout alleles</div><div class="muted" style="padding:6px 0">No loss-of-function or high-impact damaging alleles found in this gene across the panel.</div></div>`;
    const sec = Data.hasSecondaryScores(d.dataset);
    const rows = d.damaging.map(v=>{
      const open = FN.openId===v.id;
      return `<tr class="imp-row ${open?'open':''}" onclick="FUNCTION.toggle('${v.id}')">
        <td class="c-mono c-alt" style="padding-left:11px">${v.variant}</td>
        <td>${consPill(v)}</td>
        <td>${domTag(v.domain)}</td>
        <td class="num">${scoreCell(v.plantcad)}</td>
        ${sec?`<td class="num">${scoreCell(v.plantcad2)}</td>`:''}
        <td class="num">${scoreCell(v.esm)}</td>
        ${sec?`<td class="num">${scoreCell(v.esm2)}</td><td class="num">${scoreCell(v.esm3)}</td>`:''}
        <td>${prioPill(v.priority)}</td>
        <td class="num">${v.het}</td>
        <td class="num">${v.hom?`<b>${v.hom}</b>`:'0'}</td>
        <td class="num">${(v.af*100).toFixed(1)}%</td>
      </tr>${open?carrierRow(v,sec):''}`;
    }).join('');
    return `<div class="card pad">
      <div class="fn-h" style="display:flex;align-items:center;gap:10px">Damaging &amp; knockout alleles
        <span class="muted" style="font-weight:400;font-size:12px">${d.damaging.length} alleles · click a row for carrier lines</span>
        <button class="btn" style="margin-left:auto" onclick="FUNCTION.exportCSV()">${ICONS.download||''} Export CSV</button>
      </div>
      <div class="tbl-wrap" style="max-height:none"><table class="vcf imp">
        <thead><tr><th style="padding-left:11px">Allele</th><th>Consequence</th><th>Domain</th><th class="num">PlantCAD1</th>${sec?'<th class="num">PlantCAD2</th>':''}<th class="num">ESM</th>${sec?'<th class="num">ESM2</th><th class="num">ESM3</th>':''}<th>Priority</th><th class="num" title="heterozygous carriers">Het</th><th class="num" title="homozygous carriers">Hom</th><th class="num">AF</th></tr></thead>
        <tbody>${rows}</tbody>
      </table></div>
    </div>`;
  }
  function carrierRow(v, sec){
    const chip = (id,cls)=>`<span class="carrier ${cls}">${esc(id)}</span>`;
    const homs = v.carriersHom.slice(0,60).map(id=>chip(id,'hom')).join('');
    const hets = v.carriersHet.slice(0,60).map(id=>chip(id,'het')).join('');
    return `<tr class="fn-carriers"><td colspan="${sec?12:9}">
      <div class="fn-cwrap">
        <div><div class="fn-k">Homozygous ${v.consClass==='lof'?'(candidate knockouts)':''} · ${v.carriersHom.length}</div>
          <div class="fn-chips">${homs||'<span class="muted">none</span>'}${v.carriersHom.length>60?` <span class="muted">+${v.carriersHom.length-60} more</span>`:''}</div></div>
        <div style="margin-top:8px"><div class="fn-k">Heterozygous · ${v.carriersHet.length}</div>
          <div class="fn-chips">${hets||'<span class="muted">none</span>'}${v.carriersHet.length>60?` <span class="muted">+${v.carriersHet.length-60} more</span>`:''}</div></div>
      </div></td></tr>`;
  }

  /* ---------- public handlers ---------- */
  window.FUNCTION = {
    load(){ const el=document.getElementById('fnGeneInput'); if(!el)return; const g=el.value.trim(); if(!g)return;
      FN.gene=g; FN.data=null; FN.openId=null; render(); },
    toggle(id){ FN.openId = (FN.openId===id?null:id); render(); },
    exportCSV(){
      const d=FN.data; if(!d||!d.damaging) return;
      const sec=Data.hasSecondaryScores(d.dataset);
      const cols=['variant','consequence','domain','plantcad']
        .concat(sec?['plantcad2']:[])
        .concat(['esm'])
        .concat(sec?['esm2','esm3']:[])
        .concat(['combined','priority','het','hom','af','homozygous_carriers','het_carriers']);
      const line=v=>{
        const base=[v.variant,v.consequence,v.domain,v.plantcad];
        if(sec) base.push(v.plantcad2);
        base.push(v.esm);
        if(sec) base.push(v.esm2, v.esm3);
        base.push(v.combined,v.priority,v.het,v.hom,v.af.toFixed(4),
          '"'+v.carriersHom.join(';')+'"','"'+v.carriersHet.join(';')+'"');
        return base.map(x=>x==null?'':x).join(',');
      };
      const csv=[cols.join(',')].concat(d.damaging.map(line)).join('\n');
      const blob=new Blob([csv],{type:'text/csv'}), u=URL.createObjectURL(blob), a=document.createElement('a');
      a.href=u; a.download=`snpfunction_${d.gene}_damaging.csv`; document.body.appendChild(a); a.click(); a.remove(); setTimeout(()=>URL.revokeObjectURL(u),1500);
    },
  };

  function injectCSS(){
    if (document.getElementById('snpfunction-css')) return;
    const s=document.createElement('style'); s.id='snpfunction-css';
    s.textContent=`
      .fn-grid{display:flex;gap:26px;flex-wrap:wrap}
      .fn-k{font-size:10.5px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.4px;margin-bottom:4px}
      .fn-v{font-size:14px;color:var(--ink)} .fn-v .dom-tag{margin:0 4px 4px 0;display:inline-block}
      .fn-h{font-weight:600;font-size:15px;margin-bottom:12px;color:var(--ink)}
      .fn-stats{display:flex;gap:14px;flex-wrap:wrap;margin-bottom:14px}
      .fn-stat{background:#f5f8fc;border:1px solid var(--line);border-radius:10px;padding:10px 14px;min-width:120px}
      .fn-statv{font-size:20px;font-weight:700;color:var(--ink);font-family:var(--mono)}
      .fn-statk{font-size:11px;color:var(--muted);margin-top:2px}
      .fn-bar{display:flex;height:16px;border-radius:8px;overflow:hidden;background:#eef1f5}
      .fn-seg{display:inline-block;height:100%}
      .fn-seg.lof,.fn-sw.lof{background:#c0362c}.fn-seg.splice,.fn-sw.splice{background:#b8862b}
      .fn-seg.missense,.fn-sw.missense{background:#2f6ad0}.fn-seg.indel,.fn-sw.indel{background:#176c3a}
      .fn-seg.syn,.fn-sw.syn{background:#9aa7bd}.fn-seg.other,.fn-sw.other{background:#cdd6e3}
      .fn-legend{display:flex;gap:14px;flex-wrap:wrap;margin-top:8px;font-size:12px;color:var(--muted)}
      .fn-lg{display:inline-flex;align-items:center;gap:5px} .fn-sw{width:10px;height:10px;border-radius:3px;display:inline-block}
      .fn-af{margin-top:12px;font-size:13px;color:var(--muted)} .fn-af b{color:var(--ink)}
      .fn-carriers td{background:#fbfcfe;border-bottom:1px solid var(--line)}
      .fn-cwrap{padding:10px 12px} .fn-chips{display:flex;flex-wrap:wrap;gap:5px;margin-top:5px}
      .carrier{font-family:var(--mono);font-size:11px;padding:2px 7px;border-radius:6px;border:1px solid var(--line)}
      .carrier.hom{background:#fdecea;border-color:#f0c4bd;color:#8f281c}
      .carrier.het{background:#eef4ff;border-color:#cfe0ff;color:#274b8f}`;
    document.head.appendChild(s);
  }

  SNPTools.register('snpfunction', { render });
})();
