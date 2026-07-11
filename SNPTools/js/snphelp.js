/* =====================================================================
 *  snphelp.js — Help & FAQ page for the SNPTools suite.
 *
 *  Registers itself as the 'help' tool so the router can reach it, and
 *  exposes openHelp() for the MaizeGDB masthead link (which doesn't want
 *  to depend on TOOLS[] being populated). Explains every page in the
 *  suite and every dataset defined in data.js, plus a scores/annotations
 *  glossary and a short FAQ.
 *
 *  Depends on core.js (SNPTools, S, go, ICONS, renderNav) and — for the
 *  live dataset table — data.js (Data.datasets). Degrades gracefully if
 *  Data isn't present.
 * ===================================================================== */
const SNPHelp = (function () {

  /* ---- tool descriptions, in user language, kept in step with core.js ---- */
  const LIVE = 'live', SOON = 'in development';
  const PAGES = [
    { id:'snpversity', name:'SNPVersity', icon:'dna', color:'#2563eb', status:LIVE,
      tag:'Build a variant view across accessions',
      what:'The front door of the suite and the starting point for most work. Choose a dataset, type a genomic interval (or a B73 v5 gene model ID), and pick the accessions you want. SNPVersity queries the variant store and returns a color-coded genotype table plus a downloadable VCF — allele states, predicted effects, and DNA/protein language-model scores included.',
      give:'A dataset, a region or gene, and a set of accessions.',
      get:'A genotype table and a VCF. From here, "Send selection to…" hands the same result to any other tool.' },
    { id:'snptrait', name:'SNPTrait', icon:'leaf', color:'#1f8a4c', status:SOON,
      tag:'Connect variation to traits',
      what:'Links genomic variation to phenotype and trait records from the National Germplasm collection. Search, sort, and filter accessions by trait values and metadata, then move a selected set straight into SNPVersity and the rest of the suite.',
      give:'Trait, phenotype, or metadata filters across 20,000+ accessions.',
      get:'A shortlisted set of lines you can hand off to the genomic tools.' },
    { id:'snpimpact', name:'SNPImpact', icon:'star', color:'#7c3aed', status:LIVE,
      tag:'Rank candidate variants',
      what:'Prioritizes the variants in a region regardless of which accessions you picked. It orders them by an AI-based score (PlantCAD DNA-model and ESM protein-model predictions) combined with predicted consequence and Pfam domain annotation, so likely causal changes rise to the top. Filter by consequence, priority, score, or domain, and flag a shortlist.',
      give:'A region sent from SNPVersity (accessions are ignored here).',
      get:'A ranked, filterable variant table and a shortlist of candidate alleles.' },
    { id:'snpfunction', name:'SNPFunction', icon:'func', color:'#2563eb', status:LIVE,
      tag:'Gene function & allele mining',
      what:'A gene-scoped dossier that is independent of any one region. Enter a gene and it summarizes the Pfam domains and protein, computes the gene\u2019s variant burden across the whole panel, and lists which accessions carry each damaging or knockout allele.',
      give:'A gene model ID and a dataset.',
      get:'A domain/protein summary, a burden breakdown, and a damaging-allele catalog with carrier lines.' },
    { id:'snpcompare', name:'SNPCompare', icon:'compare', color:'#0e7490', status:LIVE,
      tag:'Similarity to a focal accession',
      what:'Ranks every accession by identity-by-state similarity to one focal accession. Similarity is the fraction of co-called sites with the same genotype. Works genome-wide (precomputed), for the current region (computed in-browser from a SNPVersity result), or both side-by-side with a delta that surfaces region-specific relatedness such as introgression.',
      give:'A focal accession, and optionally a region handed off from SNPVersity.',
      get:'A ranked similarity table with project, SRA ID, and accession name.' },
    { id:'snptree', name:'SNPTree', icon:'tree', color:'#15803d', status:LIVE,
      tag:'Local phylogeny',
      what:'Builds a local phylogenetic tree from the genotype matrix already in memory, using identity-by-state distances (UPGMA / neighbour-joining). Useful for reading haplotype structure, introgression, and how lines cluster in a region.',
      give:'A SNPVersity result set (up to 250 accessions before a warning).',
      get:'An interactive tree, downloadable as Newick, MEGA, or PHYLIP.' },
    { id:'snpmatrix', name:'SNPMatrix', icon:'grid', color:'#b45309', status:LIVE,
      tag:'Pairwise distance matrix',
      what:'Computes the pairwise identity-by-state distance among your selected accessions and draws it as a heatmap. Reorder by input order or by clustering, switch between IBS distance and % identity, and colour rows by bioproject.',
      give:'A SNPVersity result set (up to 400 accessions before a warning).',
      get:'A heatmap plus downloads: CSV distance matrix, PHYLIP, PNG, and SVG.' },
    { id:'snpimpute', name:'SNPImpute', icon:'impute', color:'#0891b2', status:SOON,
      tag:'Impute sequence & function',
      what:'Pan-genome\u2013guided imputation that fills missing genotypes and carries functional predictions through, across light and deep sequencing depths.',
      give:'A genotype set with missing calls.',
      get:'Imputed genotypes with function predictions.' },
    { id:'snpfold', name:'SNPFold', icon:'fold', color:'#be185d', status:LIVE,
      tag:'Variants on protein structure',
      what:'Maps coding variants onto predicted protein structure. A linear protein browser aligns variants with Pfam domains, secondary structure, and per-residue pLDDT confidence; an on-demand 3D viewer shows the fold coloured by confidence, domain, or impact; and a per-variant readout interprets each change (domain, local confidence, secondary structure, predicted \u0394\u0394G).',
      give:'A gene with an available structure model.',
      get:'A structure-aware, per-variant interpretation of coding changes.' },
    { id:'snpdensity', name:'SNPDensity', icon:'density', color:'#9333ea', status:SOON,
      tag:'Density & burden',
      what:'Measures SNP and INDEL density, burden, and distribution across genes and regions to highlight mutational load, constraint, and diversification.',
      give:'A region or gene set.',
      get:'Density tracks and burden summaries.' },
    { id:'snpgermplasm', name:'SNPGermplasm', icon:'germ', color:'#16a34a', status:SOON,
      tag:'Collection management',
      what:'Applies genotype-driven analytics to germplasm management, identifying redundancy, uniqueness, and priority materials for curation and deployment.',
      give:'A collection of accessions.',
      get:'Redundancy, uniqueness, and priority flags.' },
  ];

  /* ---- glossary ---- */
  const GLOSSARY = [
    ['B73 RefGen v5', 'Every coordinate, gene model, and annotation in the suite is anchored to the B73 version 5 maize reference genome.'],
    ['Genotype dosage', 'Each call is read as 0 (0/0, reference), 1 (heterozygous), or 2 (1/1, alternate homozygous). Missing calls (./.) are left out of a comparison rather than counted as a match.'],
    ['IBS distance / % identity', 'Identity-by-state compares two accessions site by site over the calls they share. Distance is the mean allele difference; % identity is 100 − distance. Used by SNPTree, SNPMatrix, and SNPCompare.'],
    ['Predicted effect & impact', 'Each variant carries a predicted consequence (missense, LOF, splice, indel, synonymous, …) rolled up into an impact tier: HIGH, MODERATE, LOW, or MODIFIER, most severe wins when a site lists several.'],
    ['PlantCAD score', 'A DNA language-model prediction of how disruptive a change is. MaizeGDB 2026 carries a second-generation PlantCAD2; older datasets carry a single DNA score.'],
    ['ESM score', 'A protein language-model prediction of the effect of an amino-acid change. MaizeGDB 2026 carries ESM2 and ESM3 alongside the original; older datasets carry a single AA score.'],
    ['MAF', 'Minor allele frequency, the frequency of the less common allele. SNPVersity can filter a region by a minimum MAF.'],
    ['Pfam domain', 'When a variant falls inside a known protein domain, that domain is shown and linked to InterPro. Domain annotation is still being loaded for some regions, where it reads as \u2014.'],
    ['VCF', 'The Variant Call Format file SNPVersity generates for your query. It is the exact matrix the other tools reuse when you send a selection.'],
  ];

  /* ---- FAQ ---- */
  const FAQ = [
    ['Where does the data come from?',
     'A query sends your region and accession list to the server, which reads the real HDF5 variant store, writes a VCF for exactly that slice, and returns it. The tools parse that VCF into the tables and matrices you see, so everything downstream is one consistent result.'],
    ['How do I move a selection between tools?',
     'Run a query in SNPVersity, then use "Send selection to…" in the top bar (or the buttons on a result). The same genotype matrix is handed to SNPImpact, SNPCompare, SNPTree, SNPMatrix, and SNPTrait without re-querying. SNPMatrix and SNPTree can also pass their set on to each other.'],
    ['Why did a wide region give me a download instead of a table?',
     'Regions larger than one million bases skip the in-browser table and return a downloadable VCF instead, so the page stays responsive. Narrow the interval to get the interactive table back.'],
    ['Why are some cells \u2014 or N/A?',
     'A \u2014 in a domain column means Pfam annotation for that position hasn\u2019t been loaded yet. Blank second-generation scores (PlantCAD2, ESM2, ESM3) mean the dataset is not MaizeGDB 2026 — only that family carries them. MQ and coverage read N/A when the source didn\u2019t record them.'],
    ['Is there a limit on how many accessions I can compare?',
     'The distance tools warn before doing very large computations in the browser: SNPTree above 250 accessions and SNPMatrix above 400. You can build anyway, it just may be slow. SNPImpact renders up to 1,500 variants at a time.'],
    ['What can I download?',
     'A VCF from SNPVersity; a CSV distance matrix, PHYLIP, PNG, and SVG from SNPMatrix; Newick, MEGA, and PHYLIP trees from SNPTree. Comparison and impact tables can be exported from their own pages.'],
    ['Which tools are ready to use now?',
     'SNPVersity, SNPImpact, SNPFunction, SNPCompare, SNPTree, SNPMatrix, and SNPFold are live. SNPTrait, SNPImpute, SNPDensity, and SNPGermplasm are on the roadmap and marked in development in the sidebar.'],
  ];

  function esc(s){ return String(s==null?'':s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }
  function ico(k){ return (typeof ICONS!=='undefined' && ICONS[k]) || ''; }

  /* ---- dataset table (live from data.js when available) ---- */
  function datasets(){
    try { if (typeof Data!=='undefined' && Data.datasets) return Data.datasets(); } catch(e){}
    return [];
  }
  function flag(on){ return on ? '<span class="hp-yes">Yes</span>' : '<span class="hp-no">\u2013</span>'; }
  function datasetTable(){
    const ds = datasets();
    if (!ds.length) return '';
    const rows = ds.map(d=>`<tr>
        <td><b>${esc(d.name)}</b><div class="hp-sub">${esc(d.sub||'')}</div></td>
        <td class="hp-mono">${esc(d.ref)}</td>
        <td class="hp-num">${esc(d.acc)}</td>
        <td class="hp-num">${esc(d.sites)}</td>
        <td>${(d.filters||[]).map(f=>`<span class="hp-chip">${esc(f)}</span>`).join(' ')}</td>
        <td class="hp-c">${flag(d.het)}</td>
        <td class="hp-c">${flag(d.indel)}</td>
        <td class="hp-c">${flag(d.impute)}</td>
      </tr>`).join('');
    return `<div class="hp-tablewrap"><table class="hp-table">
      <thead><tr>
        <th>Dataset</th><th>Reference</th><th>Accessions</th><th>Sites</th>
        <th>Filters</th><th>Het</th><th>INDELs</th><th>Imputed</th>
      </tr></thead><tbody>${rows}</tbody></table></div>
      <p class="hp-fine">Accessions are grouped into families (MaizeGDB 2026, MaizeGDB 2024, Schnable 2023, NAM 2021). Only the <b>MaizeGDB 2026</b> family carries the second-generation language-model scores (PlantCAD2, ESM2, ESM3); the others provide a single DNA score and a single protein score.</p>`;
  }

  function statusPill(s){
    const live = s===LIVE;
    return `<span class="hp-status ${live?'on':'off'}">${live?'Live':'In development'}</span>`;
  }
  function pageCard(p){
    const live = p.status===LIVE;
    const openBtn = live
      ? `<button class="hp-open" onclick="go('${p.id}')">Open ${esc(p.name)} ${ico('caret')}</button>`
      : `<span class="hp-open muted" title="On the roadmap">Coming soon</span>`;
    return `<details class="hp-tool">
      <summary>
        <span class="hp-ti" style="background:${p.color}">${ico(p.icon)}</span>
        <span class="hp-tt"><b>${esc(p.name)}</b><span class="hp-tag">${esc(p.tag)}</span></span>
        ${statusPill(p.status)}
        <span class="hp-chev">${ico('caret')}</span>
      </summary>
      <div class="hp-tbody">
        <p>${esc(p.what)}</p>
        <div class="hp-io">
          <div><span class="hp-lbl">You give it</span>${esc(p.give)}</div>
          <div><span class="hp-lbl">You get back</span>${esc(p.get)}</div>
        </div>
        <div class="hp-actions">${openBtn}</div>
      </div>
    </details>`;
  }

  function render(page){
    injectCSS();
    page = page || document.getElementById('page');
    page.className = 'page fade';
    const crumb=document.getElementById('crumbTool'); if(crumb) crumb.innerHTML='<b>Help &amp; FAQ</b>';

    const toolCards = PAGES.map(pageCard).join('');
    const gloss = GLOSSARY.map(g=>`<div class="hp-gl"><dt>${esc(g[0])}</dt><dd>${esc(g[1])}</dd></div>`).join('');
    const faq = FAQ.map(f=>`<details class="hp-faq"><summary>${esc(f[0])}<span class="hp-chev">${ico('caret')}</span></summary><div>${esc(f[1])}</div></details>`).join('');

    page.innerHTML = `
      <section class="hp-hero">
        <div class="hp-eyebrow">SNPTools · SNPVersity 2.1 · B73 RefGen v5</div>
        <h1>Help &amp; FAQ</h1>
        <p>SNPTools is an integrated suite for exploring maize sequence variation. Everything starts from a genomic
           query and flows between tools without re-running it — this page explains each tool, the datasets behind
           them, and the vocabulary you\u2019ll meet along the way.</p>
        <div class="hp-jump">
          <a href="#hp-flow">How it works</a>
          <a href="#hp-tools">The tools</a>
          <a href="#hp-data">Datasets</a>
          <a href="#hp-gloss">Glossary</a>
          <a href="#hp-faq">FAQ</a>
        </div>
      </section>

      <section id="hp-flow" class="hp-sec">
        <div class="hp-h"><span class="hp-n">01</span><h2>How the suite fits together</h2></div>
        <div class="hp-flow">
          <div class="hp-step"><span class="hp-si" style="background:#2563eb">${ico('search')}</span>
            <b>Query</b><p>In SNPVersity, choose a dataset, a region or gene, and the accessions you care about.</p></div>
          <div class="hp-arrow">${ico('caret')}</div>
          <div class="hp-step"><span class="hp-si" style="background:#1f8a4c">${ico('table')}</span>
            <b>Result</b><p>You get a genotype table and a VCF for exactly that slice of the genome.</p></div>
          <div class="hp-arrow">${ico('caret')}</div>
          <div class="hp-step"><span class="hp-si" style="background:#b45309">${ico('compare')}</span>
            <b>Send onward</b><p>Hand the same matrix to any other tool with "Send selection to…" — no re-query.</p></div>
        </div>
        <p class="hp-fine">A query is run against the real variant store and returned as a VCF; every other tool reuses that
          one result, so a set you build once stays consistent as you rank it, compare it, cluster it, or draw it.</p>
      </section>

      <section id="hp-tools" class="hp-sec">
        <div class="hp-h"><span class="hp-n">02</span><h2>The tools</h2></div>
        <p class="hp-lead">Eleven tools across five stages. Seven are live today; the rest are on the roadmap and share
          the same data and coordinates. Expand any tool for what it does, what it takes, and what it returns.</p>
        <div class="hp-tools">${toolCards}</div>
      </section>

      <section id="hp-data" class="hp-sec">
        <div class="hp-h"><span class="hp-n">03</span><h2>Datasets</h2></div>
        <p class="hp-lead">Each query runs against one dataset. All are called against the B73 v5 reference; they differ
          in how they were filtered, how many accessions and sites they hold, and which score columns they carry.</p>
        ${datasetTable()}
      </section>

      <section id="hp-gloss" class="hp-sec">
        <div class="hp-h"><span class="hp-n">04</span><h2>Scores &amp; annotations</h2></div>
        <dl class="hp-gloss">${gloss}</dl>
      </section>

      <section id="hp-faq" class="hp-sec">
        <div class="hp-h"><span class="hp-n">05</span><h2>Frequently asked</h2></div>
        <div class="hp-faqs">${faq}</div>
      </section>

      <section class="hp-foot">
        <div>
          <b>Still stuck?</b>
          <p>SNPTools is part of MaizeGDB, the Maize Genetics and Genomics Database.</p>
        </div>
        <div class="hp-foot-btns">
          <button class="hp-open" onclick="go('snpversity')">Start in SNPVersity ${ico('caret')}</button>
          <a class="hp-open ghost" href="https://maizegdb.org" target="_blank" rel="noopener">Visit MaizeGDB ${ico('caret')}</a>
        </div>
      </section>`;

    // jump links: smooth scroll within the .scroll container
    page.querySelectorAll('.hp-jump a').forEach(a=>{
      a.addEventListener('click', e=>{
        const t=document.querySelector(a.getAttribute('href')); if(!t) return;
        e.preventDefault(); t.scrollIntoView({behavior:'smooth', block:'start'});
      });
    });
  }

  /* Entry point used by the MaizeGDB masthead link. Mirrors go() but doesn't
     assume 'help' is present in the TOOLS map, so core.js needs no changes. */
  function open(){
    if (typeof S!=='undefined') S.tool='help';
    if (typeof closeMenu==='function') closeMenu();
    if (typeof toggleRail==='function') toggleRail(false);
    const sb=document.getElementById('sendBtn'); if(sb) sb.style.display='none';
    if (typeof renderNav==='function') renderNav();
    document.querySelectorAll('.mg-help').forEach(el=>el.classList.add('on'));
    window.scrollTo(0,0);
    const sc=document.querySelector('.scroll'); if(sc) sc.scrollTop=0;
    render(document.getElementById('page'));
  }

  function injectCSS(){
    if (document.getElementById('snphelp-css')) return;
    const s=document.createElement('style'); s.id='snphelp-css';
    s.textContent = `
      .hp-hero{padding:6px 0 18px;border-bottom:1px solid var(--line,#e6e9ef);margin-bottom:26px}
      .hp-eyebrow{font-family:var(--mono,'IBM Plex Mono',monospace);font-size:11px;letter-spacing:.14em;text-transform:uppercase;color:var(--blue-600,#2563eb);margin-bottom:10px}
      .hp-hero h1{font-family:var(--disp,'Space Grotesk',sans-serif);font-size:40px;line-height:1.05;letter-spacing:-.02em;margin:0 0 12px;color:var(--ink,#141922)}
      .hp-hero p{max-width:70ch;color:var(--muted,#5b6b83);font-size:15px;line-height:1.6;margin:0 0 18px}
      .hp-jump{display:flex;gap:8px;flex-wrap:wrap}
      .hp-jump a{font-size:12.5px;font-weight:600;color:var(--ink,#141922);text-decoration:none;padding:7px 13px;border:1px solid var(--line,#e6e9ef);border-radius:999px;background:#fff;transition:.15s}
      .hp-jump a:hover{border-color:var(--blue-600,#2563eb);color:var(--blue-600,#2563eb)}

      .hp-sec{margin:0 0 42px;scroll-margin-top:16px}
      .hp-h{display:flex;align-items:baseline;gap:12px;margin-bottom:16px}
      .hp-h h2{font-family:var(--disp,'Space Grotesk',sans-serif);font-size:23px;letter-spacing:-.01em;margin:0;color:var(--ink,#141922)}
      .hp-n{font-family:var(--mono,'IBM Plex Mono',monospace);font-size:12px;font-weight:600;color:#c8a24a;border:1px solid #e7d6a8;background:#fdf8ec;border-radius:6px;padding:3px 7px}
      .hp-lead{max-width:74ch;color:var(--muted,#5b6b83);font-size:14px;line-height:1.6;margin:0 0 16px}
      .hp-fine{max-width:78ch;color:var(--muted,#5b6b83);font-size:12.5px;line-height:1.6;margin:12px 0 0}
      .hp-fine b,.hp-lead b{color:var(--ink,#141922)}

      /* flow */
      .hp-flow{display:flex;align-items:stretch;gap:10px;flex-wrap:wrap}
      .hp-step{flex:1 1 210px;border:1px solid var(--line,#e6e9ef);border-radius:12px;background:#fff;padding:16px}
      .hp-step b{display:block;font-family:var(--disp,'Space Grotesk',sans-serif);font-size:15px;color:var(--ink,#141922);margin:10px 0 5px}
      .hp-step p{margin:0;font-size:13px;line-height:1.55;color:var(--muted,#5b6b83)}
      .hp-si{display:inline-flex;width:34px;height:34px;border-radius:9px;align-items:center;justify-content:center;color:#fff}
      .hp-si svg{width:19px;height:19px}
      .hp-arrow{display:flex;align-items:center;color:var(--faint,#aab4c4)}
      .hp-arrow svg{width:22px;height:22px}

      /* tool accordions */
      .hp-tools{display:flex;flex-direction:column;gap:10px}
      .hp-tool{border:1px solid var(--line,#e6e9ef);border-radius:12px;background:#fff;overflow:hidden}
      .hp-tool[open]{border-color:#cdd6e6;box-shadow:0 1px 0 rgba(20,25,34,.03)}
      .hp-tool summary{display:flex;align-items:center;gap:13px;padding:14px 16px;cursor:pointer;list-style:none}
      .hp-tool summary::-webkit-details-marker{display:none}
      .hp-ti{display:inline-flex;width:34px;height:34px;border-radius:9px;align-items:center;justify-content:center;color:#fff;flex:0 0 auto}
      .hp-ti svg{width:19px;height:19px}
      .hp-tt{display:flex;flex-direction:column;gap:2px;min-width:0}
      .hp-tt b{font-family:var(--disp,'Space Grotesk',sans-serif);font-size:15.5px;color:var(--ink,#141922)}
      .hp-tag{font-size:12.5px;color:var(--muted,#5b6b83)}
      .hp-status{margin-left:auto;font-size:10.5px;font-weight:700;letter-spacing:.04em;text-transform:uppercase;padding:3px 9px;border-radius:999px;border:1px solid}
      .hp-status.on{color:#1f8a4c;border-color:#bfe3cc;background:#f0faf3}
      .hp-status.off{color:#b06a12;border-color:#f0dcbb;background:#fdf6ec}
      .hp-chev{color:var(--faint,#aab4c4);transition:transform .18s;flex:0 0 auto}
      .hp-chev svg{width:18px;height:18px;transform:rotate(90deg)}
      .hp-tool[open] .hp-chev svg{transform:rotate(-90deg)}
      .hp-tbody{padding:2px 16px 18px 63px}
      .hp-tbody>p{margin:0 0 14px;font-size:13.5px;line-height:1.62;color:#3a465a;max-width:76ch}
      .hp-io{display:flex;gap:24px;flex-wrap:wrap;padding:13px 15px;background:var(--blue-50,#f2f6ff);border:1px solid #e4ecfb;border-radius:10px}
      .hp-io>div{flex:1 1 240px;font-size:13px;line-height:1.5;color:#3a465a}
      .hp-lbl{display:block;font-size:10px;font-weight:700;letter-spacing:.08em;text-transform:uppercase;color:var(--blue-600,#2563eb);margin-bottom:4px}
      .hp-actions{margin-top:14px}
      .hp-open{display:inline-flex;align-items:center;gap:6px;font-family:inherit;font-size:13px;font-weight:600;color:#fff;background:var(--ink,#141922);border:1px solid var(--ink,#141922);border-radius:9px;padding:9px 15px;cursor:pointer;text-decoration:none;transition:.15s}
      .hp-open:hover{background:#000}
      .hp-open svg{width:15px;height:15px}
      .hp-open.ghost{background:#fff;color:var(--ink,#141922)}
      .hp-open.ghost:hover{border-color:var(--blue-600,#2563eb);color:var(--blue-600,#2563eb)}
      .hp-open.muted{background:#f4f6fa;border-color:var(--line,#e6e9ef);color:var(--muted,#5b6b83);cursor:default}

      /* dataset table */
      .hp-tablewrap{overflow-x:auto;border:1px solid var(--line,#e6e9ef);border-radius:12px;background:#fff}
      .hp-table{border-collapse:collapse;width:100%;font-size:13px;min-width:720px}
      .hp-table th{text-align:left;font-size:10.5px;letter-spacing:.06em;text-transform:uppercase;color:var(--muted,#5b6b83);font-weight:700;padding:11px 14px;border-bottom:1px solid var(--line,#e6e9ef);background:#fafbfd;white-space:nowrap}
      .hp-table td{padding:12px 14px;border-bottom:1px solid #eef1f6;vertical-align:top;color:#3a465a}
      .hp-table tr:last-child td{border-bottom:none}
      .hp-sub{font-size:11.5px;color:var(--muted,#5b6b83);margin-top:2px;font-weight:400}
      .hp-mono{font-family:var(--mono,'IBM Plex Mono',monospace);font-size:12px;white-space:nowrap}
      .hp-num{font-family:var(--mono,'IBM Plex Mono',monospace);font-size:12.5px;white-space:nowrap;color:var(--ink,#141922)}
      .hp-c{text-align:center}
      .hp-chip{display:inline-block;font-size:11px;color:#3a465a;background:#f1f4f9;border:1px solid #e4e9f1;border-radius:6px;padding:2px 7px;margin:2px 2px 0 0;white-space:nowrap}
      .hp-yes{color:#1f8a4c;font-weight:600}
      .hp-no{color:var(--faint,#aab4c4)}

      /* glossary */
      .hp-gloss{margin:0;display:grid;grid-template-columns:repeat(auto-fit,minmax(310px,1fr));gap:2px 26px}
      .hp-gl{padding:13px 0;border-bottom:1px solid #eef1f6}
      .hp-gl dt{font-family:var(--disp,'Space Grotesk',sans-serif);font-size:14px;font-weight:600;color:var(--ink,#141922);margin-bottom:4px}
      .hp-gl dd{margin:0;font-size:13px;line-height:1.55;color:var(--muted,#5b6b83)}

      /* faq */
      .hp-faqs{display:flex;flex-direction:column;gap:8px}
      .hp-faq{border:1px solid var(--line,#e6e9ef);border-radius:10px;background:#fff}
      .hp-faq summary{display:flex;align-items:center;gap:12px;padding:14px 16px;cursor:pointer;font-family:var(--disp,'Space Grotesk',sans-serif);font-size:14.5px;font-weight:600;color:var(--ink,#141922);list-style:none}
      .hp-faq summary::-webkit-details-marker{display:none}
      .hp-faq .hp-chev{margin-left:auto}
      .hp-faq>div{padding:0 16px 16px;font-size:13.5px;line-height:1.62;color:#3a465a;max-width:80ch}

      /* footer */
      .hp-foot{display:flex;align-items:center;gap:18px;flex-wrap:wrap;justify-content:space-between;border:1px solid var(--line,#e6e9ef);border-radius:14px;background:linear-gradient(180deg,#fff,#fbfcfe);padding:20px 22px;margin-bottom:24px}
      .hp-foot b{font-family:var(--disp,'Space Grotesk',sans-serif);font-size:16px;color:var(--ink,#141922)}
      .hp-foot p{margin:4px 0 0;font-size:13px;color:var(--muted,#5b6b83)}
      .hp-foot-btns{display:flex;gap:10px;flex-wrap:wrap}

      @media (max-width:620px){
        .hp-hero h1{font-size:31px}
        .hp-tbody{padding-left:16px}
        .hp-arrow{transform:rotate(90deg);align-self:center}
      }`;
    document.head.appendChild(s);
  }

  if (typeof SNPTools!=='undefined') SNPTools.register('help', { render });
  return { render, open };
})();
if (typeof window!=='undefined'){ window.SNPHelp = SNPHelp; window.openHelp = SNPHelp.open; }
