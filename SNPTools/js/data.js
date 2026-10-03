/* =====================================================================
 *  data.js — the DATA LAYER (LIVE / real HDF5 build).
 *
 *  queryVariants() no longer fabricates demo data: it POSTs the region +
 *  accession list to processForm.php, which runs h5_to_vcf.py against the
 *  real .h5 store, writes a VCF, and returns its path. We then fetch that
 *  VCF and parse it into the exact row shape the tools already expect.
 *
 *  Accession IDs come from the accession catalogue (window.SNP_CATALOG, or the
 *  older flat window.SNP_REAL_ACCESSIONS), which holds the actual column
 *  names inside the .h5 files, so a selection maps to real HDF5 columns.
 *
 *  The "Domain" column reads Pfam blocks by genomic position from
 *  data/domains/ (ensureDomains / domainAt). SNPImpact ranks the queried rows
 *  (rankImpact); SNPFold and SNPFunction read one gene's rows (geneFunction).
 *  Nothing in this layer fabricates data.
 * ===================================================================== */
const Data = (function () {

  /* =============================================================
   *  BACKEND CONFIG — edit these two lines if your paths differ.
   * ============================================================= */
  const CFG = {
    endpoint : 'processForm.php',  // relative to index.html
    vcfDir   : 'vcf/',             // web-served, writable; MUST match processForm.php
    geneEndpoint : 'lookupGeneModel.php',  // gene model -> coordinates
    structDir : 'data/structures/',  // per-gene structure-<gene>.js files (SNPFold)
    domainsDir : 'data/domains/by_chr/',              // per-chromosome Pfam files: <chr>.json (preferred)
    domainsUrl : 'data/domains/domains.by_chr.json',  // combined file (fallback if per-chrom absent)
    domainsGeneUrl : 'data/domains/domains.by_gene.json',  // gene -> canonical protein domains (SNPImpact detail)
    geneModelsDir : 'data/genemodels/by_chr/',        // per-chromosome exon/CDS structure (SNPImpact gene model)
    // Hard span backstop — a region wider than this is download-only regardless
    // of how few accessions are selected. The real gate is the variants ×
    // accessions "table work" budget below (TABLE_WORK_MAX / TABLE_ROWS_MAX);
    // this just catches an absurd region.
    tableMaxSpan : 20_000_000,
    // The largest request the server builds (variants x accessions): h5_to_vcf.py's
    // SNPTOOLS_MAX_CELLS default. Only used to warn in the run bar; the server decides.
    buildCellsMax : 2e9,
  };

  /* ---------------- datasets (UI cards) ----------------
   * `id`     -> sent to processForm.php as dataSet
   * `family` -> which real accession list to show + which .h5 family
   */
  const DATASETS = [
    /* Initial release: the GRIN-linked 2026 set only. The MaizeGDB 2026 sets are kept here,
       commented out, with everything that serves them (processForm.php, ibsCompare.php,
       FIELD_STATUS, FAMILY_META, the accession catalog); uncomment to offer them again.
    {id:'mgdb2026_hq', family:'mgdb2026',     name:'MaizeGDB 2026', sub:'High Quality',   ref:'B73 v5', acc:'2,710', sites:'98M',
     filters:['MQ ≥ 30','Coverage ≥ 50%','LD max R² > 0.5'], het:true,  indel:true,  impute:false},
    {id:'mgdb2026_hc', family:'mgdb2026',     name:'MaizeGDB 2026', sub:'High Coverage',  ref:'B73 v5', acc:'2,710', sites:'290M',
     filters:['MQ ≥ 30','Coverage ≥ 50%'], het:true,  indel:true,  impute:false},
    */
    // GRIN-linked release v1.4: 926 Grzybowski et al. 2023 lines (Beagle-imputed)
    // + 7 lines called separately at the same sites (incl. NAM founder CML103). Sample ids = release
    // sample_id (ZmG_<genotype>); metadata in js/zmgrin.catalog.js.
    {id:'zmgrin2026_imp', family:'zmgrin2026', name:'MaizeGDB GRIN-linked 2026', sub:'Imputed (Grzybowski 2023 sites)',
     ref:'B73 v5', acc:'933', sites:'46M',
     filters:['Grzybowski 2023 GATK filters','Beagle 5 imputation','GRIN-linked'], het:true, indel:true, impute:true},
  ];

  /* SNPVersity annotation columns, in table order. annotationFields(datasetId) says, per
     set, whether a column is filled ('ok'), not available for that set ('na': shown empty,
     with the reason in the header tooltip and the note above the table), waiting for a data
     merge ('pending'), or left out of that set's table altogether ('hidden': a field the
     call set never records, so an always-empty column says nothing). 'auto' (Site QC): shown
     as 'ok' when a row of the result carries a value, left out otherwise, so a store without
     a site-QC sidecar shows the table exactly as before. */
  const ANNOTATION_COLUMNS = [
    {key:'gene',   label:'Gene model'},   {key:'effect', label:'Effect'},
    {key:'impact', label:'SNPEff Impact'},{key:'domain', label:'Domain'},
    {key:'mq',     label:'MQ'},           {key:'comp',   label:'COMP'},
    {key:'r2',     label:'maxR²'},        {key:'maf',    label:'MAF'},
    {key:'qc',     label:'Site QC', status:'auto'},
    {key:'pc1',    label:'PlantCAD1'},    {key:'pc2',    label:'PlantCAD2'},
    {key:'evo2',   label:'Evo2'},
    {key:'esm1',   label:'ESM1'},         {key:'esm2',   label:'ESM2'},
    {key:'esm3',   label:'ESM3'},         {key:'esmc',   label:'ESM-C'},
  ];
  const ZMGRIN_NA_QC = 'Not available for the Grzybowski et al. 2023 call set (source VCFs carry no per-site MQ/coverage).';
  const FIELD_STATUS = {
    mgdb2026_hq: {
      evo2: {status:'na', note:'Not in the MaizeGDB 2026 INFO.'},
      esmc: {status:'na', note:'Not in the MaizeGDB 2026 INFO.'},
    },
    mgdb2026_hc: {
      r2: {status:'na', note:'maxR² is computed only for the High Quality set (its LD filter); the High Coverage set has no MAXR2.'},
      evo2: {status:'na', note:'Not in the MaizeGDB 2026 INFO.'},
      esmc: {status:'na', note:'Not in the MaizeGDB 2026 INFO.'},
    },
    zmgrin2026_imp: {
      mq:   {status:'hidden', note:ZMGRIN_NA_QC},
      comp: {status:'hidden', note:ZMGRIN_NA_QC + ' (COV/CVC/CVP)'},
      r2:   {status:'ok', note:'From the 933 release genotypes: highest PLINK 1.9 r² with any variant 400-5,000 bp away, no MAF/missingness/r² filtering. Blank = no partner variant within 400-5,000 bp, or monomorphic.'},
      pc1:  {status:'ok', note:'PlantCAD1 scores for SNPs (indels have none), rounded to 0.1.'},
      pc2:  {status:'ok', note:'PlantCAD2 scores for SNPs (indels have none), rounded to 0.1.'},
      esm1: {status:'ok', note:'Missense variants only (ESM-1b 650M).'},
      esm2: {status:'ok', note:'Missense variants only (ESM-2 650M, Full_ESM_stack store layer = esm2_store_score).'},
      esm3: {status:'ok', note:'Missense variants only (ESM3 open).'},
      evo2: {status:'ok', note:'Evo2 7B log-likelihood ratio (256-bp left context), SNPs within 1 kb of a gene only, rounded to 0.1.'},
      esmc: {status:'ok', note:'ESM C 600M log-likelihood ratio, missense variants only.'},
      maf:  {status:'ok', note:'Computed from the release genotypes of the 933 lines.'},
      qc:   {status:'auto', note:'Heterozygous and homozygous carriers among all 933 lines (tools/build_site_qc.py sidecars).'},
    },
  };
  /* Precomputed genome-wide IBS files per dataset family (ibsCompare.php?probe=1).
     mgdb2026 keeps its fixed files in ./distance/ and is always offered, as before;
     any other family is offered only when ./distance/<family>/ has the files.
     Resolves {dataset, available, n, sites:['all','snp'?], trees:['nj','upgma'?]}. */
  const _gdProbe = {};
  function globalDistance(family){
    if (family === 'mgdb2026') return Promise.resolve({dataset:family, available:true, sites:['all'], trees:[], legacy:true});
    if (!_gdProbe[family]) _gdProbe[family] = fetch('ibsCompare.php?probe=1&dataset=' + encodeURIComponent(family), {cache:'no-store'})
      .then(r => r.ok ? r.json() : {available:false})
      .catch(() => ({available:false}))
      .then(j => { const v = Object.assign({dataset:family, available:false, sites:[], trees:[]}, j || {});
                   _gdProbe[family]._v = v; return v; });
    return _gdProbe[family];
  }
  function globalDistanceKnown(family){
    if (family === 'mgdb2026') return {dataset:family, available:true, sites:['all'], trees:[], legacy:true};
    return _gdProbe[family] && _gdProbe[family]._v || null;
  }
  function annotationFields(datasetId){
    const over = FIELD_STATUS[datasetId] || {};
    return ANNOTATION_COLUMNS.map(c => Object.assign({key:c.key, label:c.label, status:c.status || 'ok', note:''}, over[c.key] || {}));
  }

  /*const DATASETS = [
    {id:'mgdb2026_hq', family:'mgdb2026',     name:'MaizeGDB 2026', sub:'High Quality',   ref:'B73 v5', acc:'2,710', sites:'98M',
     filters:['MQ ≥ 30','Coverage ≥ 50%','LD max R² > 0.5'], het:true,  indel:true,  impute:false},
    {id:'mgdb2026_hc', family:'mgdb2026',     name:'MaizeGDB 2026', sub:'High Coverage',  ref:'B73 v5', acc:'2,710', sites:'290M',
     filters:['MQ ≥ 30','Coverage ≥ 50%'], het:true,  indel:true,  impute:false},
    {id:'schnable2023',family:'schnable2023', name:'Schnable 2023', sub:'Imputed markers',ref:'B73 v5', acc:'1,515', sites:'12M',
     filters:['Imputed'], het:false, indel:false, impute:true},
    {id:'nam2021',     family:'nam2021',      name:'NAM 2021',      sub:'Founder panel',  ref:'B73 v5', acc:'27',    sites:'78M',
     filters:['MQ ≥ 30','Founder panel'], het:true,  indel:true,  impute:false},
    {id:'mgdb2024_hq', family:'mgdb2024',     name:'MaizeGDB 2024', sub:'High Quality',   ref:'B73 v5', acc:'1,498', sites:'83M',
     filters:['MQ ≥ 30','Coverage ≥ 50%'], het:true,  indel:true,  impute:false},
  ];
  */

  /* colors + a friendly label per family, for the accession picker group */
  const FAMILY_META = {
    mgdb2026:     {name:'MaizeGDB 2026', color:'#2563eb'},
    mgdb2024:     {name:'MaizeGDB 2024', color:'#cf8a12'},
    schnable2023: {name:'Schnable 2023', color:'#1f8a4c'},
    nam2021:      {name:'NAM 2021',      color:'#7c3aed'},
    zmgrin2026:   {name:'MaizeGDB GRIN-linked 2026', color:'#0f766e'},
  };

  /* ---------------- accession catalog (projects -> groups -> accessions) ----------------
   * Source: window.SNP_CATALOG (compiled from data/accessions.tsv + data/projects.tsv).
   * Falls back to the older flat window.SNP_REAL_ACCESSIONS if present.
   */
  const CATALOG = (typeof window !== 'undefined' && window.SNP_CATALOG) || null;
  const REAL    = (typeof window !== 'undefined' && window.SNP_REAL_ACCESSIONS) || {};

  /* An unknown id falls back to the first listed dataset's family (it was a hard-coded
     'mgdb2026', which pointed every tool at a set this release does not offer). */
  function familyOf(datasetId){
    const d = DATASETS.find(x => x.id === datasetId);
    return d ? d.family : DATASETS[0].family;
  }
  // The second-generation language-model scores (PlantCAD2 / ESM2 / ESM3) are carried by
  // the GRIN-linked 2026 release (and by MaizeGDB 2026, when listed). Tools call this to
  // show those columns conditionally; ESM-C is in the GRIN-linked INFO as well.
  function hasSecondaryScores(datasetId){ return ['mgdb2026','zmgrin2026'].indexOf(familyOf(datasetId)) >= 0; }

  /* Language-model score columns for a dataset, in display order. `key` indexes
     the row object parseVcf() builds (pc1/pc2 = DNA, esm1..esm3/esmc = protein).
     Used by SNPGeo's variant table; any tool can render row[model.key]. Families
     without second-generation scores get the single DNA/protein pair their
     INFO carries (DNA_SCORE -> pc1, AA_SCORE -> esm1). */
  function scoreModels(datasetId){
    if (!hasSecondaryScores(datasetId)) return [
      {key:'pc1',  kind:'dna',     label:'PlantCAD', tip:'PlantCaduceus DNA language-model score (INFO DNA_SCORE / plantcad1_score).'},
      {key:'esm1', kind:'protein', label:'ESM1b',    tip:'ESM1b protein language-model score (INFO AA_SCORE / ESM1_score).'},
    ];
    return [
      {key:'pc1',  kind:'dna',     label:'PlantCAD1', tip:'PlantCAD1 DNA language-model score (INFO plantcad1_score, or DNA_SCORE; 0.1 steps in the GRIN-linked 2026 set).'},
      {key:'pc2',  kind:'dna',     label:'PlantCAD2', tip:'PlantCAD2 DNA language-model score (INFO plantcad2_score; 0.1 steps in the GRIN-linked 2026 set).'},
      {key:'evo2', kind:'dna',     label:'Evo2',      tip:'Evo2 7B DNA language-model log-likelihood ratio (INFO evo2_score; SNPs within 1 kb of a gene; 0.1 steps).'},
      {key:'esm1', kind:'protein', label:'ESM1',      tip:'ESM1b protein language-model score (INFO ESM1_score, or AA_SCORE).'},
      {key:'esm2', kind:'protein', label:'ESM2',      tip:'ESM2 protein language-model score (INFO ESM2_score).'},
      {key:'esm3', kind:'protein', label:'ESM3',      tip:'ESM3 protein language-model score (INFO ESM3_score).'},
      {key:'esmc', kind:'protein', label:'ESM-C',     tip:'ESM-C protein language-model score (INFO ESMC_score).'},
    ];
  }
  function famNode(datasetId){
    return CATALOG ? CATALOG.families[familyOf(datasetId)] : null;
  }

  // Projects (bioproject sections) with metadata + groups, for the picker.
  function projectsFor(datasetId){
    const fam = famNode(datasetId);
    if (fam) return fam.projects;
    // legacy fallback: single synthetic project from the flat list
    const list = REAL[familyOf(datasetId)] || [];
    return [{id:familyOf(datasetId), title:familyOf(datasetId), bioprojects:[], color:'#2563eb',
             count:list.length, namFounders:[],
             groups:[{name:null, accessions:list}]}];
  }

  // Flat accession list (one object per accession) for search, chips, table headers.
  const _accCache = {};
  function accessionsFor(datasetId){
    const fam = familyOf(datasetId);
    if (_accCache[fam]) return _accCache[fam];
    const node = famNode(datasetId);
    let out = [];
    if (node){
      node.projects.forEach(p => p.groups.forEach(g => g.accessions.forEach(a => {
        // Spread first so every catalogue metadata key (panel, subpop, country,
        // grin, ...) reaches SNPTrait/SNPGeo; the computed keys below still win.
        out.push({...a, id:a.id, run:a.run, founder:a.founder, rep:a.rep, reps:a.reps,
                  label:a.label, group:g.name, namFounder:a.namFounder,
                  proj:p.id, projColor:p.color, projTitle:p.title});
      })));
    } else {
      out = (REAL[fam] || []).map(a => ({...a}));
    }
    /* Line QC (js/<prefix>.lineqc.js, when loaded): the share of clean sites where the line is
       heterozygous, and its class. The generated catalogue file is left as it is. */
    const lq = (typeof window !== 'undefined' && window.SNP_LINE_QC) ? window.SNP_LINE_QC[fam] : null;
    if (lq && lq.lines) out.forEach(a => { const x = lq.lines[a.id];
      if (x){ a.hetShare = x[0]; a.sampleQC = LINE_QC_LABELS[x[1]] || null; } });
    _accCache[fam] = out;
    return out;
  }
  const LINE_QC_LABELS = {I:'Inbred', E:'Elevated heterozygosity', H:'Heterozygous sample'};

  function namFoundersFor(datasetId){
    const fam = famNode(datasetId);
    return fam ? (fam.namFounders || []) : [];
  }

  // Default selection: for a dataset with tagged NAM founders (2026), preselect
  // one accession per NAM founder; otherwise one per the first 12 founders.
  function defaultSelectionFor(datasetId){
    const list = accessionsFor(datasetId);
    const nam = namFoundersFor(datasetId);
    if (nam.length){
      const pick = {}, ids = [];
      for (const a of list){
        if (a.namFounder && !pick[a.namFounder]){ pick[a.namFounder] = 1; ids.push(a.id); }
      }
      return ids;
    }
    const ids = [], seen = new Set();
    for (const a of list){
      if (seen.has(a.founder)) continue;
      seen.add(a.founder); ids.push(a.id);
      if (ids.length >= 12) break;
    }
    return ids;
  }

  // union index for accessionById()
  let _byId = null;
  function idIndex(){
    if (_byId) return _byId;
    _byId = new Map();
    ['mgdb2026','mgdb2024','schnable2023','nam2021','zmgrin2026'].forEach(fam => {
      const dsid = (DATASETS.find(d => d.family === fam) || {}).id;
      if (dsid) accessionsFor(dsid).forEach(a => { if (!_byId.has(a.id)) _byId.set(a.id, a); });
    });
    return _byId;
  }
  function accessionById(id){ return idIndex().get(id) || null; }

  /* ---------------- genome geometry ---------------- */
  /* Example gene intervals, B73 v5 (Zm00001eb.1), as returned by
     lookupGeneModel.php from gff/genes_data.serialized. The previous values were
     placeholders that did not match the annotation (e.g. Zm00001eb374090 is on
     chr9, not chr8:163.45 Mb). Live lookups always go through lookupGene(). */
  const GENE_MODELS = {
    'Zm00001eb374090':{chr:'chr9', start:12838008,  end:12843999},
    'Zm00001eb067740':{chr:'chr2', start:4493424,   end:4497434},
    'Zm00001eb404760':{chr:'chr10',start:218406,    end:220251},
    'Zm00001eb404740':{chr:'chr10',start:129631,    end:131683},
    'Zm00001eb233650':{chr:'chr5', start:90721578,  end:90727950},
    'Zm00001eb313510':{chr:'chr7', start:123685735, end:123691964},
  };
  const CHR_LEN = {chr1:308452471,chr2:243675191,chr3:238017767,chr4:250330460,chr5:226353449,
    chr6:181357234,chr7:185808916,chr8:182411202,chr9:163004744,chr10:152435371};
  const CENTRO = {chr10:.34};

  /* ---------------- gene model -> coordinates (live) ----------------
   * Resolves a B73 v5 gene model ID to its interval via lookupGeneModel.php,
   * which reads the serialized GFF store on the server. Returns
   *   {id, chr, start, end}  or  null when the ID isn't found.
   */
  const EXAMPLE_GENES = [
    'Zm00001eb374090','Zm00001eb067740','Zm00001eb374230',
    'Zm00001eb056510','Zm00001eb233650','Zm00001eb313510',
  ];
  /* One lookup per gene id for the session: SNPFold alone asked three or four times per gene,
     and every ask unserializes the whole gene store on the server. A failed lookup is not
     kept, so it is retried next time. Callers get their own copy of the answer. */
  const _geneLookups = new Map();
  async function lookupGene(id){
    id = (id || '').trim();
    if (!id) return null;
    if (!_geneLookups.has(id)){
      const p = _lookupGeneUncached(id);
      _geneLookups.set(id, p);
      p.catch(() => { _geneLookups.delete(id); });
    }
    const g = await _geneLookups.get(id);
    return g ? Object.assign({}, g) : null;
  }
  async function _lookupGeneUncached(id){
    const url = `${CFG.geneEndpoint}?geneModelId=${encodeURIComponent(id)}`;
    const resp = await fetch(url, {cache:'no-store'});
    if (!resp.ok) throw new Error('Gene lookup failed (HTTP ' + resp.status + ')');
    const raw = await resp.text();
    let d;
    try { d = JSON.parse(raw); }
    catch (e) { throw new Error('lookupGeneModel.php did not return JSON:\n' + raw.slice(0, 600)); }
    // The PHP returns id:'empty' (and 0/0) when the gene isn't in the store.
    if (!d || d.id === 'empty' || d.chromosome == null) return null;
    const start = parseInt(d.start, 10), end = parseInt(d.end, 10);
    if (!Number.isFinite(start) || !Number.isFinite(end) || (start === 0 && end === 0)) return null;
    return {id, chr: d.chromosome, start, end};
  }

  /* =============================================================
   *  SNPVERSITY — LIVE query.
   *  region + accession ids -> VCF (via processForm.php) -> rows
   * ============================================================= */
  const SEVERITY = {HIGH:3, MODERATE:2, LOW:1, MODIFIER:0};

  function uniqueOutName(lo, hi){
    const ts = Date.now();
    const rnd = Math.random().toString().slice(2, 11);
    return `${CFG.vcfDir}snpv_${ts}_${rnd}_${lo}_${hi}.vcf.gz`;
  }

  /* processForm.php returns a gzip-compressed VCF (.vcf.gz). Inflate it in
     the browser; fall back to treating the bytes as plain text if they are
     not gzip-framed (a proxy that already decompressed, or a plain .vcf
     from some other caller). */
  async function readVcfResponse(resp){
    const buf = new Uint8Array(await resp.arrayBuffer());
    const gzipped = buf.length > 2 && buf[0] === 0x1f && buf[1] === 0x8b;
    if (gzipped && typeof DecompressionStream === 'function'){
      const stream = new Blob([buf]).stream().pipeThrough(new DecompressionStream('gzip'));
      return await new Response(stream).text();
    }
    return new TextDecoder().decode(buf);
  }

  /* Parse the whole INFO column (';'-separated key=value tags) once per row
     into a plain object. This replaced ~20 calls/row to a helper that
     compiled a fresh RegExp each time — millions of compilations for a wide
     region. A missing key reads back as undefined, which numOrNull /
     firstNum / cleanTok already treat the same as null. */
  function parseInfo(infoStr){
    const o = {};
    if (!infoStr) return o;
    const parts = infoStr.split(';');
    for (let i = 0; i < parts.length; i++){
      const p = parts[i];
      const eq = p.indexOf('=');
      if (eq === -1) o[p] = ''; else o[p.slice(0, eq)] = p.slice(eq + 1);
    }
    return o;
  }
  const cleanTok = s => {
    if (s == null) return '';
    const tok = String(s).split(/[;,]+/)[0].replace(/_/g, ' ').trim();
    return tok === '.' ? '' : tok;
  };
  function numOrNull(v){
    if (v == null || v === '.' || v === '') return null;
    const f = parseFloat(v);
    return Number.isNaN(f) ? null : f;
  }
  function intOrNull(v){
    if (v == null || v === '' || v === '.') return null;
    const n = parseInt(v, 10);
    return Number.isNaN(n) ? null : n;
  }
  function firstNum(v){                       // first of a possibly comma-listed value
    if (v == null) return null;
    return numOrNull(String(v).split(/[;,]+/)[0]);
  }

  /* ---- Site QC ----
     A store with a site-QC sidecar gives every variant NHET and NHOM, its heterozygous and
     homozygous-alternate carriers among ALL the release's lines, and SITEQC (h5_to_vcf.py). The
     class comes from the two counts by one rule, first match wins; the same rule is
     site_qc_codes in tools/build_site_qc.py and h5_to_vcf.py and site_qc in
     tools/annotate_release_info.py. Usable = PASS or HET_ELEVATED (or no class at all);
     flagged = HET_ONLY or HET_EXCESS. */
  function siteQc(nHet, nHom){
    if (nHet == null || nHom == null) return null;
    const c = nHet + nHom;
    if (c === 0) return 'NO_CARRIER';
    if (nHom === 0) return 'HET_ONLY';
    if (nHet > nHom) return 'HET_EXCESS';
    if (nHet * 4 >= c) return 'HET_ELEVATED';
    return 'PASS';
  }
  const SITE_QC = {
    PASS:         {label:'Pass',         group:'pass',      tip:'Fewer than a quarter of the carriers are heterozygous.'},
    HET_ELEVATED: {label:'Het elevated', group:'caution',   tip:'A quarter or more of the carriers are heterozygous. Read with caution.'},
    HET_EXCESS:   {label:'Het excess',   group:'flagged',   tip:'Heterozygous carriers outnumber homozygous ones. In inbred lines this usually means reads from another copy of the sequence map here.'},
    HET_ONLY:     {label:'Het only',     group:'flagged',   tip:'No line is homozygous for this allele.'},
    NO_CARRIER:   {label:'No carrier',   group:'nocarrier', tip:'No line in this release carries the allele. The site comes from the larger panel the site list was built on.'},
  };
  function qcGroup(code){ return (code && SITE_QC[code]) ? SITE_QC[code].group : null; }
  function qcUsable(code){ return code == null || code === 'PASS' || code === 'HET_ELEVATED'; }
  function qcFlagged(code){ return code === 'HET_ONLY' || code === 'HET_EXCESS'; }
  /* The views' filter choices: 'all' (All sites), 'usable' (Hide flagged and no-carrier),
     'pass' (Passing only). A row with no class passes every choice. */
  function qcKeep(code, mode){
    if (code == null || !mode || mode === 'all') return true;
    return mode === 'pass' ? code === 'PASS' : qcUsable(code);
  }
  /* Each view's default choice, in one place. SNPFunction ('func') always applies it; the
     others only to results that carry the class. */
  const SITE_QC_DEFAULTS = {versity:'all', impact:'usable', func:'usable', fold:'usable', geo:'all', distance:'usable'};
  const SITE_QC_CHOICES = [['all', 'All sites'], ['usable', 'Hide flagged and no-carrier'], ['pass', 'Passing only']];
  function releaseLines(){ return (accessionsFor(DATASETS[0].id) || []).length || 933; }
  function qcTip(code, nHet, nHom){
    const q = SITE_QC[code]; if (!q) return '';
    return q.tip + (nHet != null && nHom != null
      ? ` ${nHet.toLocaleString()} heterozygous, ${nHom.toLocaleString()} homozygous among ${releaseLines().toLocaleString()} lines.` : '');
  }
  /* rows whose class is usable, or every row when none carries a class (SNPTree, SNPMatrix,
     SNPCompare region scope) */
  function hasSiteQc(rows){ return !!(rows && rows.some(r => r.qc != null)); }
  function usableRows(rows){
    if (!hasSiteQc(rows)) return rows || [];
    return rows.filter(r => qcUsable(r.qc));
  }

  /* ---- SNPCurate (window.SNP_CURATE, js/snpcurate.data.js) ----
     Curated alleles by site, by gene and by interval, for the marks the other tools draw; nothing
     when the file is not loaded. The index follows the loaded entries, so it is rebuilt when an
     entry is added at run time. */
  let _curIdx = null;
  function curateIndex(){
    const C = (typeof window !== 'undefined' && window.SNP_CURATE) ? window.SNP_CURATE : null;
    const E = (C && C.entries) || null;
    if (!E) return null;
    if (_curIdx && _curIdx.E === E && _curIdx.n === E.length) return _curIdx;
    const site = new Map(), gene = new Map();
    E.forEach(e => {
      if (e.status === 'site') site.set(`${e.chr}:${e.pos}:${e.ref}:${e.alt}`, e);
      [e.gene, e.site_gene].filter(Boolean).forEach(g => {
        if (!gene.has(g)) gene.set(g, []);
        if (gene.get(g).indexOf(e) < 0) gene.get(g).push(e);
      });
    });
    _curIdx = {E, n:E.length, site, gene};
    return _curIdx;
  }
  /* the site entry at exactly this chr, pos, ref and alt, or null */
  function curatedAt(chr, pos, ref, alt){
    const I = curateIndex();
    return I ? (I.site.get(`${chr}:${pos}:${ref}:${alt}`) || null) : null;
  }
  /* every entry whose gene or site_gene is this model, sites and not */
  function curatedForGene(gene){ const I = curateIndex(); return I ? (I.gene.get(gene) || []).slice() : []; }
  /* entries whose site, position, codon or interval falls in chr:lo-hi */
  function curatedInInterval(chr, lo, hi){
    const I = curateIndex(); if (!I) return [];
    return I.E.filter(e => e.chr === chr && (
      (e.pos != null && e.pos >= lo && e.pos <= hi) ||
      (e.codon_positions || []).some(p => p >= lo && p <= hi) ||
      (e.start != null && e.end != null && e.start <= hi && e.end >= lo)));
  }
  const CURATE_STATUS = {not_a_site:'Not a site in this release.', not_locatable:'Cannot be placed on this release\u2019s gene model.',
    structural:'A structural variant this release cannot genotype.', unreliable:'Calls in this gene are unreliable in this release.',
    not_annotated:'The gene model is not in this release\u2019s annotation.'};
  /* one sentence on why the release cannot show a grey entry ('' for a site) */
  function curateStatusText(e){ return (e && CURATE_STATUS[e.status]) || ''; }

  /* A region's genotype matrix is variants × accessions — by far the biggest
     thing in memory. Store each call as a 1-byte code in an Int8Array rather
     than a string:  0 = ref/ref, 1 = het, 2 = alt/alt, 3 = missing (null).
     ~8x smaller than an Array of strings, and the downstream dosage() calls
     become the identity. GT_CODE is the fast path — a direct lookup on the
     raw VCF cell, no allocation. */
  const GT_CODE = Object.create(null);
  GT_CODE['0/0'] = 0; GT_CODE['0|0'] = 0;
  GT_CODE['1/0'] = 1; GT_CODE['0/1'] = 1; GT_CODE['1|0'] = 1; GT_CODE['0|1'] = 1;
  GT_CODE['1/1'] = 2; GT_CODE['1|1'] = 2;
  GT_CODE['./.'] = 3; GT_CODE['.|.'] = 3; GT_CODE['.'] = 3;
  function gtCode(cell){
    if (cell == null) return 3;
    const q = GT_CODE[cell];
    if (q !== undefined) return q;
    // rare: FORMAT subfields (0/1:...) or non-0/1 allele indices
    let s = cell;
    const c = s.indexOf(':'); if (c !== -1) s = s.slice(0, c);
    const q2 = GT_CODE[s]; if (q2 !== undefined) return q2;
    const m = s.split(/[/|]/);
    if (m.length < 2 || m[0] === '.' || m[1] === '.' || m[0] === '' || m[1] === '') return 3;
    return (m[0] !== '0' ? 1 : 0) + (m[1] !== '0' ? 1 : 0);
  }

  function parseVcf(text, ids){
    const rows = [];
    if (!text) return {rows, sampleCols:{}, header:[]};
    const sampleCols = {};      // sampleName -> column index in the row
    const nIds = ids.length;
    const N = text.length;
    /* Walk lines by index rather than text.split('\n') — the split would
       hold a second full copy of the payload as an array of line strings. */
    let p = 0;
    while (p < N){
      let nl = text.indexOf('\n', p);
      if (nl === -1) nl = N;
      let end = nl;
      if (end > p && text.charCodeAt(end - 1) === 13) end--;   // strip trailing \r
      if (end === p){ p = nl + 1; continue; }                  // blank line
      if (text.charCodeAt(p) === 35){                          // '#'
        if (text.charCodeAt(p + 1) !== 35){                    // '#CHROM' header row, not '##meta'
          const cols = text.slice(p + 1, end).split('\t');
          for (let i = 9; i < cols.length; i++) sampleCols[cols[i]] = i;
        }
        p = nl + 1; continue;
      }
      const t = text.slice(p, end).split('\t');
      p = nl + 1;
      if (t.length < 8) continue;
      const II = parseInfo(t[7] || '');

      // impact: take the most severe when a site lists several
      let impact = 'MODIFIER';
      const effTokens = (II.EFFECT || '').split(/[;,]+/);
      let best = -1;
      for (let k = 0; k < effTokens.length; k++){
        const key = effTokens[k].trim().toUpperCase();
        if (key in SEVERITY && SEVERITY[key] > best){ best = SEVERITY[key]; impact = key; }
      }

      // genotypes in the exact order the caller selected, as 1-byte codes
      const gts = new Int8Array(nIds);
      for (let k = 0; k < nIds; k++){
        const ci = sampleCols[ids[k]];
        gts[k] = (ci == null) ? 3 : gtCode(t[ci]);
      }

      const mq = firstNum(II.MQ), cvp = firstNum(II.CVP), pos = parseInt(t[1], 10);
      const nHet = intOrNull(II.NHET), nHom = intOrNull(II.NHOM), qv = II.SITEQC;
      const qc = (qv && SITE_QC[qv]) ? qv : siteQc(nHet, nHom);
      rows.push({
        pos,
        ref:    t[3],
        alt:    t[4],
        gene:   cleanTok(II.GENEMODEL) || '—',
        effect: cleanTok(II.TYPE) || 'intergenic',
        impact,
        sub:    cleanTok(II.SUB),
        // a site with several consequences (several genes): the raw parallel lists, read
        // per gene by rowForGene(); null for the usual single consequence
        anns:   (II.GENEMODEL && II.GENEMODEL.indexOf(',') !== -1) ? [II.GENEMODEL, II.TYPE, II.EFFECT, II.SUB] : null,
        domain: domainAt(t[0], pos),     // Pfam domain covering this position (or '—')
        mq:     (mq != null) ? Math.round(mq) : 'N/A',
        comp:   (cvp != null) ? cvp : 'N/A',
        r2:     firstNum(II.MAXR2),
        maf:    firstNum(II.MAF),         // null when the store has no MAF (shown empty, not 0)
        // site QC (h5_to_vcf.py, from the store's sidecar): carriers among ALL the release's lines
        nHet, nHom, qc,
        // 2026 uses plantcad1/2 + ESM1/2/3; older projects (2024/Schnable/NAM)
        // use a single DNA_SCORE (PlantCaduceus) and AA_SCORE (ESM1b) -> map to col 1.
        pc1:    numOrNull(II.plantcad1_score != null ? II.plantcad1_score : II.DNA_SCORE),
        pc2:    numOrNull(II.plantcad2_score),
        evo2:   numOrNull(II.evo2_score != null ? II.evo2_score : II.EVO2_score),   // genic +/-1 kb SNPs only
        esm1:   numOrNull(II.ESM1_score != null ? II.ESM1_score : II.AA_SCORE),
        esm2:   numOrNull(II.ESM2_score),
        esm3:   numOrNull(II.ESM3_score),
        esmc:   numOrNull(II.ESMC_score != null ? II.ESMC_score : (II.ESMC_SCORE != null ? II.ESMC_SCORE : II.esmc_score)),
        gts,
      });
    }
    return {rows, sampleCols};
  }

  /* ---- one site, several consequences ----
     SnpEff writes one consequence per affected gene (and transcript): GENEMODEL, TYPE,
     EFFECT and SUB are parallel comma lists, most severe first. In the GRIN-linked 2026 test
     windows 1,464 of 3,744 sites list more than one (20 of 414 in a MaizeGDB 2026 window): a site intronic in one
     gene and downstream of its neighbour, or missense in two overlapping genes. A row's own
     gene/effect/sub are the FIRST entry, which is right for a region view (SNPVersity,
     SNPImpact). A view of one gene (SNPFunction, SNPFold, SNPGeo by gene) must read that
     gene's entry instead, or it drops the site, or shows the neighbour's substitution:
     rowForGene(). The ESM scores belong to one entry only, the first missense one of the
     list (tools/annotate_release_info.py, pick_esm), so they go with that entry alone;
     PlantCAD/Evo2 score the allele, not a protein, and stay on every entry. */
  function annotationsOf(r){
    if (!r) return [];
    const one = [{gene:r.gene, effect:r.effect, impact:r.impact, sub:r.sub, esm:true}];
    if (!r.anns) return one;
    const [G, T, E, U] = r.anns.map(x => String(x == null ? '' : x).split(','));
    if (T.length !== G.length) return one;            // not parallel: keep the first-entry reading
    const esmAt = T.findIndex(t => /missense/i.test(t));
    return G.map((g, i) => {
      const imp = String(E[i] || '').trim().toUpperCase();
      return {gene: cleanTok(g) || '—', effect: cleanTok(T[i]) || 'intergenic',
              impact: (imp in SEVERITY) ? imp : 'MODIFIER', sub: cleanTok(U[i]), esm: i === esmAt};
    });
  }
  /* The row as seen from one gene: that gene's most severe entry (the first on a tie),
     its ESM scores blanked when they were computed for another entry. null when the site
     has no consequence in the gene. Genotypes are shared, not copied. */
  function rowForGene(r, gene){
    if (!r || !gene) return null;
    if (!r.anns) return r.gene === gene ? r : null;
    let best = null;
    for (const a of annotationsOf(r)){
      if (a.gene === gene && (!best || SEVERITY[a.impact] > SEVERITY[best.impact])) best = a;
    }
    if (!best) return null;
    const v = Object.assign({}, r, {gene:best.gene, effect:best.effect, impact:best.impact, sub:best.sub});
    if (!best.esm){ v.esm1 = v.esm2 = v.esm3 = v.esmc = null; }
    return v;
  }
  /* Every single gene model a row has a consequence in (intergenic "A_B" spans dropped). */
  function genesOf(r){
    return [...new Set(annotationsOf(r).map(a => a.gene))].filter(g => g && g !== '—' && !/\s/.test(g));
  }

  // Build the accession objects the table header needs, in selection order.
  function accsFor(datasetId, ids){
    const map = new Map(accessionsFor(datasetId).map(a => [a.id, a]));
    return ids.map(id => map.get(id) || {id, run:id, founder:id, proj:familyOf(datasetId), projColor:'#8a94a6'});
  }

  // Order the submitted accessions by project (catalog order), then by accession
  // name — so each project's color renders as one contiguous block in the table.
  function sortIds(datasetId, ids){
    const order = {};
    projectsFor(datasetId).forEach((p, i) => { order[p.id] = i; });
    const map = new Map(accessionsFor(datasetId).map(a => [a.id, a]));
    const keyOf = id => {
      const a = map.get(id);
      const pi = a && (order[a.proj] != null) ? order[a.proj] : 9999;
      return {pi, name: a ? a.id : id};
    };
    return ids.slice().sort((x, y) => {
      const kx = keyOf(x), ky = keyOf(y);
      if (kx.pi !== ky.pi) return kx.pi - ky.pi;
      return kx.name.localeCompare(ky.name, undefined, {numeric:true, sensitivity:'base'});
    });
  }

  /* Beyond this much table work (variants × accessions, or raw variant count)
     the in-browser parse + render would freeze the tab for many seconds, so the
     result is download-only — no table. Chosen so a typical selection (tens of
     accessions, a few Mb) always renders, but a full panel or a genome-arm-wide
     region does not. */
  const TABLE_WORK_MAX = 4e7;
  const TABLE_ROWS_MAX  = 400000;

  /* genome-wide variant density per bp for a dataset, from its manifest `sites`
     count ('98M' etc.) over the B73 v5 assembly length — used only to *predict*
     in the run bar whether a query will come back as a table or a VCF download
     (the real decision uses the exact server count). Regional density varies,
     so this is a rough guide. */
  const GENOME_LEN = Object.keys(CHR_LEN).reduce((s, c) => s + CHR_LEN[c], 0);
  function siteCount(s){
    const m = String(s || '').match(/([\d.]+)\s*([KMG]?)/i);
    if (!m) return 0;
    return parseFloat(m[1]) * ({K:1e3, M:1e6, G:1e9}[m[2].toUpperCase()] || 1);
  }
  /* estimateResult(datasetId, spanBp, nAccessions) -> {estVariants, willDownload, overLimit} */
  function estimateResult(datasetId, spanBp, nAccessions){
    const d = DATASETS.find(x => x.id === datasetId);
    const density = d ? siteCount(d.sites) / GENOME_LEN : 0;
    const estVariants = Math.round(Math.max(0, spanBp) * density);
    const willDownload = spanBp > CFG.tableMaxSpan
      || estVariants > TABLE_ROWS_MAX
      || estVariants * Math.max(1, nAccessions) > TABLE_WORK_MAX;
    const overLimit = estVariants * Math.max(1, nAccessions) > CFG.buildCellsMax;
    return {estVariants, willDownload, overLimit};
  }

  /**
   * queryVariants(dataset, chr, lo, hi, ids, opts) -> Promise<{rows, accs, chr, vcfUrl, span, wide, empty, variants}>
   * `wide` is true when the result is too big for the table — the caller should
   * offer the VCF download only. `opts.forceTable` skips the size gate (used by
   * the internal gene-scoped callers, whose regions are always tiny).
   */
  async function queryVariants(dataset, chr, lo, hi, ids, opts){
    opts = opts || {};
    const forceTable = !!opts.forceTable;
    ids = sortIds(dataset, ids);      // group by project, then accession name
    const accs = accsFor(dataset, ids);
    const span = Math.max(hi - lo, 0);
    const outName = uniqueOutName(lo, hi);

    const body = new URLSearchParams({
      start: String(lo), end: String(hi), chr: chr, dataSet: dataset,
      genotypes: JSON.stringify(ids), outName: outName,
    });
    const resp = await fetch(CFG.endpoint, {
      method: 'POST',
      headers: {'Content-Type': 'application/x-www-form-urlencoded'},
      body,
    });
    const raw = await resp.text();
    let json;
    try { json = JSON.parse(raw); }
    catch (e) { throw new Error('processForm.php did not return JSON:\n' + raw.slice(0, 1000)); }

    // A genuinely empty interval is a valid result, not a failure.
    if (json.status === 'empty'){
      return {rows:[], accs, chr, vcfUrl:null, span, wide:false, empty:true};
    }
    if (!resp.ok || json.status !== 'success'){
      if (json.output)  console.error('[processForm.php] h5_to_vcf.py output:\n' + json.output);
      if (json.command) console.error('[processForm.php] command:\n' + json.command);
      const err = new Error(json.message || ('Request failed (HTTP ' + resp.status + ')'));
      err.detail = {command: json.command, output: json.output};
      throw err;
    }

    const vcfUrl = json.outFile || outName;
    const variants = (json.variants != null && isFinite(+json.variants)) ? +json.variants : null;

    // Too big for an in-browser table — download only. The server tells us the
    // exact variant count, so this decision happens *before* the fetch + parse.
    const tooBig = span > CFG.tableMaxSpan
      || (variants != null && (variants > TABLE_ROWS_MAX || variants * ids.length > TABLE_WORK_MAX));
    if (tooBig && !forceTable){
      return {rows:[], accs, chr, vcfUrl, span, variants, wide:true, empty:false};
    }

    // Fetch the VCF, retrying briefly — right after the server writes it, a
    // same-machine dev server can 404 it for a moment (Windows file handle) or
    // lose the race against a user-triggered download of the same URL.
    let vcfResp = null;
    for (let attempt = 0; attempt < 4; attempt++){
      if (attempt) await new Promise(r => setTimeout(r, 250 * attempt));
      try { vcfResp = await fetch(vcfUrl, {cache:'no-store'}); } catch (e) { vcfResp = null; }
      if (vcfResp && vcfResp.ok) break;
    }
    if (!vcfResp || !vcfResp.ok){
      const err = new Error('The VCF was built but could not be retrieved (' + (vcfResp ? 'HTTP ' + vcfResp.status : 'network error') + '). It is still available to download.');
      err.detail = {vcfUrl};
      throw err;
    }
    const vcfText = await readVcfResponse(vcfResp);
    await ensureDomains(chr);                    // load just this chromosome's Pfam file (cached)
    const {rows} = parseVcf(vcfText, ids);
    return {rows, accs, chr, vcfUrl, span, variants, wide:false, empty: rows.length === 0};
  }

  /**
   * queryVariantsByGene(dataset, geneId, ids, opts) -> Promise<{rows, accs, chr, start, end, gene, vcfUrl, span, wide, empty, variants}>
   * Resolves the gene through lookupGeneModel.php (gff/genes_data.serialized);
   * when that endpoint is unreachable (e.g. a static file server without PHP)
   * the built-in GENE_MODELS table covers the example genes. `ids` defaults to the
   * dataset's default selection — pass every accession id when per-population
   * statistics are needed (SNPGeo does).
   */
  async function queryVariantsByGene(dataset, geneId, ids, opts){
    dataset = dataset || DATASETS[0].id;
    geneId = String(geneId || '').trim();
    let gene = null, lookupErr = null;
    try { gene = await lookupGene(geneId); }
    catch (e) { lookupErr = e; }
    if (!gene && GENE_MODELS[geneId]) gene = Object.assign({id: geneId}, GENE_MODELS[geneId]);
    if (!gene){
      if (lookupErr) throw new Error('Gene lookup unavailable (' + lookupErr.message + ') and ' + geneId + ' is not a built-in example gene.');
      throw new Error('Gene not found: ' + geneId);
    }
    ids = (ids && ids.length) ? ids : defaultSelectionFor(dataset);
    const r = await queryVariants(dataset, gene.chr, gene.start, gene.end, ids, opts);
    return Object.assign({}, r, {chr: gene.chr, start: gene.start, end: gene.end, gene: geneId});
  }

  /* =============================================================
   *  SNPFold — protein structure files + coding variants.
   * ============================================================= */
  function structureFor(gene){ return (window.SNPFOLD_STRUCT||{})[gene] || null; }
  function pdbFor(gene){ return (window.SNPFOLD_PDB||{})[gene] || null; }
  /* Lazily load js/structures/structure-<gene>.js on demand. Resolves when the
     gene's model is available; rejects if there's no file for it. */
  function ensureStructure(gene){
    gene = (gene||'').trim();
    if (!gene) return Promise.reject(new Error('no gene'));
    if ((window.SNPFOLD_STRUCT||{})[gene] || (window.SNPFOLD_PDB||{})[gene]) return Promise.resolve(true);
    return new Promise((resolve, reject)=>{
      const s = document.createElement('script');
      s.src = `${CFG.structDir}structure-${encodeURIComponent(gene)}.js`;
      s.onload  = ()=> ((window.SNPFOLD_STRUCT||{})[gene] || (window.SNPFOLD_PDB||{})[gene])
        ? resolve(true) : reject(new Error('no model data for '+gene));
      s.onerror = ()=> reject(new Error('no structure file for '+gene));
      document.head.appendChild(s);
    });
  }

  /* ---- Pfam domains by genomic position (SNPVersity Domain column / SNPImpact) ----
     Loads ONE chromosome's file on demand (data/domains/by_chr/<chr>.json), cached
     per chromosome. Falls back to a combined domains.by_chr.json if per-chrom is absent,
     and degrades to '—' if neither exists. */
  const _domByChr = {};            // chr -> sorted intervals
  const _domChrProm = {};          // chr -> in-flight promise
  let _domCombined = null, _domCombinedProm = null;
  function _loadCombinedDomains(){
    if (_domCombined) return Promise.resolve(_domCombined);
    if (_domCombinedProm) return _domCombinedProm;
    _domCombinedProm = fetch(CFG.domainsUrl, {cache:'force-cache'})
      .then(r => r.ok ? r.json() : {}).catch(() => ({}))
      .then(x => { _domCombined = x || {}; return _domCombined; });
    return _domCombinedProm;
  }
  function ensureDomains(chr){
    if (chr in _domByChr) return Promise.resolve(_domByChr[chr]);
    if (_domChrProm[chr]) return _domChrProm[chr];
    _domChrProm[chr] = fetch(CFG.domainsDir + encodeURIComponent(chr) + '.json', {cache:'force-cache'})
      .then(r => { if (!r.ok) throw 0; return r.json(); })
      .then(arr => { _domByChr[chr] = arr || []; return _domByChr[chr]; })
      .catch(() => _loadCombinedDomains().then(idx => { _domByChr[chr] = (idx && idx[chr]) || []; return _domByChr[chr]; }));
    return _domChrProm[chr];
  }
  // rows: [g_start, g_end, name, pfam, type], sorted by g_start.
  // returns "Name (PFxxxxx)" for the most specific domain covering pos, else '—'.
  function domainAt(chr, pos){
    const a = _domByChr[chr];
    if (!a || !a.length) return '—';
    let lo = 0, hi = a.length;                 // first index with g_start > pos
    while (lo < hi){ const m = (lo + hi) >> 1; if (a[m][0] <= pos) lo = m + 1; else hi = m; }
    let best = null;
    for (let i = lo - 1; i >= 0; i--){
      const iv = a[i];
      if (iv[0] < pos - 100000) break;         // domain blocks are exon-sized
      if (iv[0] <= pos && iv[1] >= pos){
        if (!best || (iv[1] - iv[0]) < (best[1] - best[0])) best = iv;   // smallest = most specific
      }
    }
    return best ? `${best[2]}${best[3] ? ` (${best[3]})` : ''}` : '—';
  }

  /* Coding-consequence classification for the structural view.
       missense -> full residue-level treatment (AA change + LM scores)
       lof      -> truncation / position marker (stop gained, frameshift, start/stop lost)
       indel    -> in-frame insertion / deletion (position marker)
       null     -> not shown (synonymous, stop_retained, splice, UTR, intron, intergenic ...) */
  function classifyConsequence(effect){
    const e = String(effect || '').toLowerCase().replace(/[\s]+/g, '_');
    if (/missense|protein_altering|non[_-]?synonymous/.test(e)) return {klass:'missense', label:'Missense',            structural:true};
    if (/stop_gained|nonsense/.test(e))                         return {klass:'lof',      label:'Stop gained',         structural:false};
    if (/frameshift/.test(e))                                   return {klass:'lof',      label:'Frameshift',          structural:false};
    if (/start_lost|initiator_codon/.test(e))                   return {klass:'lof',      label:'Start lost',          structural:false};
    if (/stop_lost/.test(e))                                    return {klass:'lof',      label:'Stop lost',           structural:false};
    if (/inframe_insertion/.test(e))                            return {klass:'indel',    label:'In-frame insertion',  structural:false};
    if (/inframe_deletion/.test(e))                             return {klass:'indel',    label:'In-frame deletion',   structural:false};
    return null;
  }

  const AA3 = {ALA:'A',ARG:'R',ASN:'N',ASP:'D',CYS:'C',GLN:'Q',GLU:'E',GLY:'G',HIS:'H',ILE:'I',
    LEU:'L',LYS:'K',MET:'M',PHE:'F',PRO:'P',SER:'S',THR:'T',TRP:'W',TYR:'Y',VAL:'V',TER:'*',SEC:'U'};
  /* Parse INFO SUB into {ref, resi, alt}. Handles 1-letter (A1V, R441*),
     3-letter (Ala1Val, Trp441Ter) and frameshift forms. null if no residue. */
  function parseSub(sub){
    if (sub == null) return null;
    const s = String(sub).trim().replace(/^p\./i, '');
    let m = s.match(/^([A-Z*])(\d+)(fs\*?|[A-Z*]|=)?$/);
    if (m) return {ref:m[1], resi:+m[2], alt:(m[3]||'').replace('fs*','fs') || null};
    m = s.match(/^([A-Za-z]{3})(\d+)([A-Za-z]{3}|Ter|\*|fs)?/);
    if (m){
      const ref = AA3[m[1].toUpperCase()] || m[1][0].toUpperCase();
      let alt = null;
      if (m[3]) alt = AA3[m[3].toUpperCase()] || (m[3] === '*' ? '*' : (/fs/i.test(m[3]) ? 'fs' : m[3][0].toUpperCase()));
      return {ref, resi:+m[2], alt};
    }
    m = s.match(/(\d+)/);
    return m ? {ref:'', resi:+m[1], alt:null} : null;
  }
  function hgvsProtein(p, cls){
    if (!p) return null;
    switch (cls.label){
      case 'Frameshift':  return 'p.' + (p.ref || '') + p.resi + 'fs';
      case 'Stop gained': return 'p.' + (p.ref || '') + p.resi + '*';
      case 'Start lost':  return 'p.' + (p.ref || 'M') + p.resi + '?';
      case 'Stop lost':   return 'p.*' + p.resi + (p.alt && p.alt !== '*' ? p.alt : 'ext');
    }
    if (cls.klass === 'indel') return 'p.' + (p.ref || '') + p.resi + (cls.label.indexOf('insertion') >= 0 ? 'ins' : 'del');
    return 'p.' + (p.ref || '') + p.resi + (p.alt || '');   // missense
  }

  /* LIVE: coding variants for one gene, from the same HDF5 -> VCF pipeline as
     SNPVersity. Residue + AA change come from INFO SUB; consequence from TYPE;
     scores from PlantCAD/ESM (pc1/esm1). Nothing fabricated. Async.
     `ids` defaults to a representative selection (site-level INFO is panel-wide,
     so this stays small); pass a broader set for an exhaustive allele catalog. */
  async function queryFoldVariants(gene, dataset, ids){
    dataset = dataset || DATASETS[0].id;
    const g = await lookupGene(gene);
    if (!g || !g.chr) return [];
    ids = ids || defaultSelectionFor(dataset);
    let res;
    try { res = await queryVariants(dataset, g.chr, g.start, g.end, ids, {forceTable:true}); }
    catch (e){ console.warn('queryFoldVariants:', e && e.message); return []; }
    return foldVariantsFromRows(res.rows, gene);
  }
  /* The coding variants SNPFold draws, from the rows of a gene's interval. They are site-level
     (INFO only), so any accession set gives the same list: geneFunction() builds it from its
     full-panel rows, which spares SNPFold a second query of the same interval. */
  function foldVariantsFromRows(rows, gene){
    const out = [];
    for (const r0 of (rows || [])){
      // this gene's consequence at the site, even when SnpEff lists another gene first
      const r = gene ? (rowForGene(r0, gene) || ((!r0.gene || r0.gene === '—') ? r0 : null)) : r0;
      if (!r) continue;
      const cls = classifyConsequence(r.effect);
      if (!cls) continue;
      const p = parseSub(r.sub);
      if (!p || p.resi == null) continue;                 // need a residue to place it
      const pc  = (r.pc1  != null ? r.pc1  : null);
      const esm = (r.esm1 != null ? r.esm1 : null);
      const combined = (pc != null && esm != null) ? +(((pc + esm) / 2)).toFixed(2)
                     : (pc != null ? pc : (esm != null ? esm : null));
      out.push({
        id:'f' + out.length, gene, resi:p.resi, ref:p.ref, alt:p.alt,
        variant: hgvsProtein(p, cls), consequence: cls.label, consClass: cls.klass,
        structural: cls.structural, pos:r.pos, refNt:r.ref, altNt:r.alt,
        impact: r.impact || null, maf: (r.maf != null ? r.maf : null),
        plantcad: pc, esm: esm,
        plantcad2: (r.pc2 != null ? r.pc2 : null), esm2: (r.esm2 != null ? r.esm2 : null),
        esm3: (r.esm3 != null ? r.esm3 : null),
        evo2: (r.evo2 != null ? r.evo2 : null), esmc: (r.esmc != null ? r.esmc : null),
        combined, qc: r.qc, nHet: r.nHet, nHom: r.nHom,
        priority: combined == null ? null
                : (combined <= -7 ? 'TOP' : combined <= -4 ? 'HIGH' : combined <= -1 ? 'MODERATE' : 'LOW'),
      });
    }
    out.sort((a, b) => a.resi - b.resi);
    return out;
  }

  /* ---------------- SNPImpact: rank a region's variants ---------------- */
  // broader than the fold classifier: also covers splice / UTR / intron / intergenic
  function impactClass(effect){
    const e = String(effect || '').toLowerCase().replace(/[\s]+/g, '_');
    if (/missense|protein_altering|non[_-]?synonymous/.test(e)) return {klass:'missense', label:'Missense',           severe:false};
    if (/stop_gained|nonsense/.test(e))         return {klass:'lof',    label:'Stop gained',        severe:true};
    if (/frameshift/.test(e))                    return {klass:'lof',    label:'Frameshift',         severe:true};
    if (/start_lost|initiator/.test(e))          return {klass:'lof',    label:'Start lost',         severe:true};
    if (/stop_lost/.test(e))                     return {klass:'lof',    label:'Stop lost',          severe:true};
    if (/splice_(acceptor|donor)/.test(e))       return {klass:'splice', label:'Splice site',        severe:true};
    if (/splice/.test(e))                        return {klass:'splice', label:'Splice region',      severe:false};
    if (/inframe_insertion/.test(e))             return {klass:'indel',  label:'In-frame insertion', severe:false};
    if (/inframe_deletion/.test(e))              return {klass:'indel',  label:'In-frame deletion',  severe:false};
    if (/synonymous|stop_retained/.test(e))      return {klass:'syn',    label:'Synonymous',         severe:false};
    if (/intron/.test(e))                        return {klass:'other',  label:'Intron',             severe:false};
    if (/5_prime_utr|five_prime/.test(e))        return {klass:'other',  label:'5\u2032 UTR',        severe:false};
    if (/3_prime_utr|three_prime/.test(e))       return {klass:'other',  label:'3\u2032 UTR',        severe:false};
    if (/upstream/.test(e))                      return {klass:'other',  label:'Upstream',           severe:false};
    if (/downstream/.test(e))                    return {klass:'other',  label:'Downstream',         severe:false};
    if (/intergenic/.test(e))                    return {klass:'other',  label:'Intergenic',         severe:false};
    return {klass:'other', label: (cleanTok(effect) || 'Other'), severe:false};
  }
  function impactPriority(cls, combined, level){
    if (cls.severe) return 'TOP';                          // LOF, splice donor/acceptor
    if (cls.klass === 'missense' || cls.klass === 'indel'){
      if (combined != null){
        if (combined <= -7) return 'TOP';
        if (combined <= -4) return 'HIGH';
        if (combined <= -1) return 'MODERATE';
        return 'LOW';
      }
      return level === 'HIGH' ? 'HIGH' : level === 'MODERATE' ? 'MODERATE' : 'LOW';
    }
    return level === 'HIGH' ? 'HIGH' : 'LOW';               // syn / non-coding
  }
  // rank every variant in a region (from parseVcf rows). Accessions are irrelevant here.
  function rankImpact(rows){
    const out = [];
    for (const r of (rows || [])){
      const cls = impactClass(r.effect);
      const pc = r.pc1 != null ? r.pc1 : null, esm = r.esm1 != null ? r.esm1 : null;
      const combined = (pc != null && esm != null) ? +(((pc + esm) / 2)).toFixed(2)
                     : (pc != null ? pc : (esm != null ? esm : null));
      const p = parseSub(r.sub);
      const coding = (cls.klass === 'missense' || cls.klass === 'lof' || cls.klass === 'indel');
      const variant = (p && p.resi != null && coding)
        ? hgvsProtein(p, {label: cls.label, klass: cls.klass})
        : `${r.pos} ${r.ref}>${r.alt}`;
      out.push({
        id: 'i' + out.length, gene: r.gene, pos: r.pos, ref: r.ref, alt: r.alt,
        variant, consequence: cls.label, consClass: cls.klass,
        resi: p ? p.resi : null, aaRef: p ? p.ref : null, aaAlt: p ? p.alt : null,
        domain: r.domain || '\u2014', impactLevel: r.impact || null,
        plantcad: pc, esm: esm, combined,
        plantcad2: (r.pc2 != null ? r.pc2 : null),
        pc1: r.pc1, pc2: r.pc2, esm1: r.esm1, esm2: r.esm2, esm3: r.esm3,
        evo2: (r.evo2 != null ? r.evo2 : null), esmc: (r.esmc != null ? r.esmc : null),
        maf: r.maf, r2: r.r2, mq: r.mq, qc: r.qc, nHet: r.nHet, nHom: r.nHom,
        priority: impactPriority(cls, combined, r.impact), percentile: null,
      });
    }
    // real region percentile: most deleterious (most negative combined) -> highest
    const scored = out.filter(v => v.combined != null).slice().sort((a, b) => a.combined - b.combined);
    const N = scored.length;
    scored.forEach((v, i) => { v.percentile = N ? Math.round(100 * (N - i) / N) : null; });
    return out;
  }

  /* gene -> canonical protein domains (protein coords) for the SNPImpact detail track */
  let _geneDom = null, _geneDomProm = null;
  function ensureGeneDomains(){
    if (_geneDom) return Promise.resolve(_geneDom);
    if (_geneDomProm) return _geneDomProm;
    _geneDomProm = fetch(CFG.domainsGeneUrl, {cache:'force-cache'})
      .then(r => r.ok ? r.json() : {}).catch(() => ({}))
      .then(x => { _geneDom = x || {}; return _geneDom; });
    return _geneDomProm;
  }
  function geneDomains(gene){ return (_geneDom || {})[gene] || null; }

  /* per-chromosome exon/CDS structure (canonical transcript) for SNPImpact's gene-model view */
  const _gmByChr = {}, _gmProm = {};
  function ensureGeneModels(chr){
    if (chr in _gmByChr) return Promise.resolve(_gmByChr[chr]);
    if (_gmProm[chr]) return _gmProm[chr];
    _gmProm[chr] = fetch(CFG.geneModelsDir + encodeURIComponent(chr) + '.json', {cache:'force-cache'})
      .then(r => r.ok ? r.json() : {}).catch(() => ({}))
      .then(x => { _gmByChr[chr] = x || {}; return _gmByChr[chr]; });
    return _gmProm[chr];
  }
  function geneModelOf(chr, gene){ return (_gmByChr[chr] || {})[gene] || null; }

  /* ---------- SNPFunction: gene-scoped functional dossier + allele mining ---------- */
  // gts entries are now 1-byte codes from parseVcf: 0/1/2 = dosage, 3 = missing.
  function _dose(g){ return (g == null || g === 3) ? null : g; }
  function _avg(a){ return a.length ? +(a.reduce((s,x)=>s+x,0)/a.length).toFixed(2) : null; }

  // Aggregate a gene across the WHOLE panel: which accessions carry damaging/LOF alleles.
  async function geneFunction(gene, dataset){
    dataset = dataset || DATASETS[0].id;
    const dsName = (DATASETS.find(d=>d.id===dataset)||{}).name || dataset;
    const g = await lookupGene(gene);
    if (!g || !g.chr) return {gene, dataset, datasetName:dsName, error:'No gene-model coordinates found for '+gene+'.'};
    const ids = (accessionsFor(dataset)||[]).map(a=>a.id);            // FULL panel
    await Promise.all([ensureDomains(g.chr), ensureGeneDomains(), ensureGeneModels(g.chr)]);
    let res;
    try { res = await queryVariants(dataset, g.chr, g.start, g.end, ids, {forceTable:true}); }
    catch (e){ return {gene, chr:g.chr, start:g.start, end:g.end, dataset, datasetName:dsName, error:'Variant query failed: '+(e&&e.message)}; }

    const accs = res.accs || [];
    const qcDiffer = [];
    // every site with a consequence in this gene, read as this gene's consequence
    // (rowForGene): a first-listed-gene filter dropped 31 of Zm00001eb374230's 51 coding sites
    const rows = (res.rows || []).map(r => rowForGene(r, gene)).filter(Boolean);
    const gd = geneDomains(gene), gm = geneModelOf(g.chr, gene);

    const variants = rows.map((r, vi) => {
      const cls = impactClass(r.effect);
      const pc = r.pc1 != null ? r.pc1 : null, esm = r.esm1 != null ? r.esm1 : null;
      const pc2 = r.pc2 != null ? r.pc2 : null, esm2 = r.esm2 != null ? r.esm2 : null;
      const esm3 = r.esm3 != null ? r.esm3 : null;
      const evo2 = r.evo2 != null ? r.evo2 : null, esmc = r.esmc != null ? r.esmc : null;
      const combined = (pc != null && esm != null) ? +(((pc + esm) / 2)).toFixed(2) : (pc != null ? pc : (esm != null ? esm : null));
      let het=0, hom=0, called=0; const homIds=[], hetIds=[], missIds=[];
      for (let k=0;k<accs.length;k++){ const d=_dose(r.gts[k]);
        if (d==null){ missIds.push(accs[k].id); continue; } called++;
        if (d===2){ hom++; homIds.push(accs[k].id); } else if (d===1){ het++; hetIds.push(accs[k].id); } }
      const an = 2*called, ac = het + 2*hom, af = an ? ac/an : 0;
      const p = parseSub(r.sub);
      const coding = (cls.klass==='missense'||cls.klass==='lof'||cls.klass==='indel');
      const variant = (p && p.resi!=null && coding) ? hgvsProtein(p, {label:cls.label, klass:cls.klass}) : `${r.pos} ${r.ref}>${r.alt}`;
      // site QC from this variant's own full-panel counts; the store's class should agree
      const qc = siteQc(het, hom);
      if (r.qc != null && r.qc !== qc) qcDiffer.push(r.pos);
      return {id:'fx'+vi, pos:r.pos, ref:r.ref, alt:r.alt, variant, consequence:cls.label, consClass:cls.klass, severe:cls.severe,
        domain:r.domain||'\u2014', resi:p?p.resi:null, aaRef:p?p.ref:null, aaAlt:p?p.alt:null, plantcad:pc, esm, plantcad2:pc2, esm2, esm3, evo2, esmc, combined,
        priority: impactPriority(cls, combined, r.impact),
        het, hom, af, carriersHom:homIds, carriersHet:hetIds, qc, nHet:het, nHom:hom,
        nRef: called - het - hom, carriersMiss: missIds};       // reference lines: Data.referenceLines(v)
    });
    if (qcDiffer.length) console.warn(`[geneFunction] ${gene}: the store's SITEQC differs from the panel counts at ${qcDiffer.length} site(s), e.g. ${qcDiffer[0]}`);

    /* Site QC: the allele list, knockouts and burden use usable variants only (PASS or
       HET_ELEVATED); flagged and no-carrier alleles are kept apart for the page to show on demand. */
    const qcCounts = {PASS:0, HET_ELEVATED:0, HET_EXCESS:0, HET_ONLY:0, NO_CARRIER:0};
    variants.forEach(v => { qcCounts[v.qc]++; });
    const nVariable = variants.length - qcCounts.NO_CARRIER;
    const used = variants.filter(v => qcUsable(v.qc));
    const proteinChanging = v => v.consClass==='missense' || v.consClass==='lof' || v.consClass==='indel' || (v.consClass==='splice' && v.severe);
    const pcCarrier = variants.filter(v => proteinChanging(v) && v.qc !== 'NO_CARRIER');
    const pcFlagged = pcCarrier.filter(v => qcFlagged(v.qc)).length;
    const codingFlaggedShare = {share: pcCarrier.length ? +(pcFlagged / pcCarrier.length).toFixed(4) : 0,
                                flagged: pcFlagged, withCarrier: pcCarrier.length};

    const byClass = {missense:0, lof:0, splice:0, indel:0, syn:0, other:0};
    used.forEach(v => { byClass[v.consClass] = (byClass[v.consClass]||0) + 1; });
    const nonsyn = byClass.missense + byClass.lof + byClass.indel + byClass.splice, syn = byClass.syn;
    const domainDisrupting = used.filter(v => v.domain!=='\u2014' && (v.consClass==='missense'||v.consClass==='lof'||v.consClass==='indel')).length;
    const afSpectrum = {
      rare:   used.filter(v => v.af>0 && v.af<0.01).length,
      low:    used.filter(v => v.af>=0.01 && v.af<0.05).length,
      common: used.filter(v => v.af>=0.05).length,
    };
    // exon vs intron: prefer the real gene model; else fall back to consequence
    let exonic=0, intronic=0;
    for (const v of used){
      let inExon;
      if (gm && gm.exons) inExon = gm.exons.some(e => v.pos>=e[0] && v.pos<=e[1]);
      else inExon = (v.consClass!=='other') || /utr/i.test(v.consequence);
      if (/intron/i.test(v.consequence) || (gm && gm.exons && !inExon)) intronic++;
      else if (inExon) exonic++;
    }
    const PR = ['TOP','HIGH','MODERATE','LOW'];
    const allDamaging = variants
      .filter(v => v.consClass==='lof' || v.severe || (v.consClass==='missense' && v.combined!=null && v.combined<=-4))
      .sort((a,b) => (PR.indexOf(a.priority)-PR.indexOf(b.priority)) || ((a.combined==null?0:a.combined)-(b.combined==null?0:b.combined)));
    const damaging = allDamaging.filter(v => qcUsable(v.qc));
    const damagingFlagged = allDamaging.filter(v => qcFlagged(v.qc));
    const damagingNoCarrier = allDamaging.filter(v => v.qc === 'NO_CARRIER');
    const koGenotypes = damaging.filter(v=>v.consClass==='lof').reduce((n,v)=>n+v.hom, 0);
    const koLines = new Set(); damaging.filter(v=>v.consClass==='lof').forEach(v=>v.carriersHom.forEach(id=>koLines.add(id)));
    // lines homozygous only for a flagged loss-of-function allele
    const koFlagged = new Set();
    damagingFlagged.filter(v=>v.consClass==='lof').forEach(v=>v.carriersHom.forEach(id=>{ if (!koLines.has(id)) koFlagged.add(id); }));
    // SNPFold's list carries the same class, from the same counts
    const qcAt = new Map(variants.map(v => [v.pos + '|' + v.ref + '|' + v.alt, v]));
    const foldVariants = foldVariantsFromRows(res.rows, gene);
    foldVariants.forEach(f => { const v = qcAt.get(f.pos + '|' + f.refNt + '|' + f.altNt);
      if (v){ f.qc = v.qc; f.nHet = v.het; f.nHom = v.hom; } });

    return {
      gene, chr:g.chr, start:g.start, end:g.end, strand: gm?gm.strand:null, dataset, datasetName:dsName,
      nAccessions: accs.length, nVariants: variants.length, nVariable, qcCounts, codingFlaggedShare,
      protLen: gd?gd.len:null, protein: gd?gd.protein:null, domains: gd?(gd.domains||[]):[],
      burden: { byClass, nonsyn, syn, nonsynSyn: syn ? +(nonsyn/syn).toFixed(2) : (nonsyn?null:0),
                exonic, intronic, exonIntron: intronic ? +(exonic/intronic).toFixed(2) : (exonic?null:0),
                domainDisrupting, meanPlantcad:_avg(used.map(v=>v.plantcad).filter(x=>x!=null)),
                meanEsm:_avg(used.map(v=>v.esm).filter(x=>x!=null)),
                meanPlantcad2:_avg(used.map(v=>v.plantcad2).filter(x=>x!=null)),
                meanEsm2:_avg(used.map(v=>v.esm2).filter(x=>x!=null)),
                meanEsm3:_avg(used.map(v=>v.esm3).filter(x=>x!=null)),
                meanEvo2:_avg(used.map(v=>v.evo2).filter(x=>x!=null)),
                meanEsmc:_avg(used.map(v=>v.esmc).filter(x=>x!=null)), afSpectrum,
                nUsed: used.length, nSites: variants.length },
      damaging, damagingFlagged, damagingNoCarrier, koGenotypes, koLines: koLines.size, koLinesFlagged: koFlagged.size, variants,
      foldVariants,                                          // SNPFold's list, same rows
    };
  }


  /* The lines homozygous for the reference allele at a geneFunction() variant: the panel minus its
     homozygous and heterozygous carriers and its missing calls (kept as a count, nRef, and derived
     on demand rather than listed: most of the panel at most sites). */
  function referenceLines(v, dataset){
    const skip = new Set([].concat(v.carriersHom || [], v.carriersHet || [], v.carriersMiss || []));
    return (accessionsFor(dataset || DATASETS[0].id) || []).map(a => a.id).filter(id => !skip.has(id));
  }

  const DEFAULT_DS = DATASETS[0].id;
  return {
    datasets:    () => DATASETS,
    // dataset-aware accession accessors
    projectsFor, accessionsFor, defaultSelectionFor, familyOf, hasSecondaryScores, scoreModels, namFoundersFor,
    familyMeta:  () => FAMILY_META,
    // parsing helpers shared with SNPGeo (and the Node smoke test)
    parseVcf, parseSub, classifyConsequence,
    annotationsOf, rowForGene, genesOf,   // per-gene reading of multi-consequence sites
    annotationFields,       // per-set status of the SNPVersity annotation columns
    // site QC: the class rule, its labels and tips, and each view's default filter
    siteQc, SITE_QC, SITE_QC_DEFAULTS, SITE_QC_CHOICES, qcGroup, qcUsable, qcFlagged, qcKeep, qcTip, hasSiteQc, usableRows,
    // SNPCurate: curated alleles by site, gene and interval (window.SNP_CURATE)
    curatedAt, curatedForGene, curatedInInterval, curateStatusText,
    globalDistance, globalDistanceKnown,   // precomputed genome-wide IBS / trees per family
    // backwards-compatible defaults (first dataset)
    projects:    () => projectsFor(DEFAULT_DS),
    accessions:  () => accessionsFor(DEFAULT_DS),
    defaultSelection: () => defaultSelectionFor(DEFAULT_DS),
    geneModels:  () => GENE_MODELS,
    exampleGenes:() => EXAMPLE_GENES,
    lookupGene,             // async: gene model id -> {id, chr, start, end} | null
    chromLengths:() => CHR_LEN,
    centromeres: () => CENTRO,
    estimateResult,                        // run-bar prediction: table vs VCF download
    accessionById,
    queryVariants,          // now async (returns a Promise)
    queryVariantsByGene,    // async: gene id -> queryVariants() over the gene interval
    structureFor,
    pdbFor,
    ensureStructure,
    queryFoldVariants,
    rankImpact,
    ensureGeneDomains,
    geneDomains,
    ensureGeneModels,
    geneModelOf,
    geneFunction,
    referenceLines,         // lines homozygous for the reference allele at a geneFunction() variant
  };
})();

/* Global helper: the Site QC mark of a variant (Data.SITE_QC): a pill with the class label for
   the four non-pass classes, a faint "pass" for PASS, nothing without a class. The label is
   always text, so color is never the only cue; the tooltip ends with the carrier counts.
   Used by SNPVersity, SNPImpact, SNPFunction, SNPFold and SNPGeo. */
function siteQcPill(code, nHet, nHom){
  if (code == null || typeof Data === 'undefined' || !Data.SITE_QC[code]) return '';
  var q = Data.SITE_QC[code];
  var tt = String(Data.qcTip(code, nHet, nHom)).replace(/&/g,'&amp;').replace(/"/g,'&quot;').replace(/</g,'&lt;');
  return '<span class="qcp ' + q.group + '" data-tt="' + tt + '">' + (code === 'PASS' ? 'pass' : q.label) + '</span>';
}

/* Global helper: the "het" marker of a line whose sample is heterozygous as sequenced (line QC
   class H: heterozygous at more than 5% of clean sites). Takes an accession object or id;
   nothing for other lines or without js/zmgrin.lineqc.js. Used by SNPVersity's header cells and
   selected-line chips and SNPFunction's carrier chips. */
function sampleHetMark(acc){
  if (typeof Data === 'undefined') return '';
  var a = (acc && typeof acc === 'object') ? acc : Data.accessionById(acc);
  if (!a || a.sampleQC !== 'Heterozygous sample' || a.hetShare == null) return '';
  return '<span class="het-mark" data-tt="This sample is heterozygous at ' + (a.hetShare * 100).toFixed(1) +
    '% of clean sites. Inbred lines are near 0.8%.">het</span>';
}

/* Global helper: a score column's header in a table. "PlantCAD1" / "PlantCAD2" break into two
   lines (Plant / CAD1), so those columns are no wider than the other score columns (Evo2, ESM1...),
   whose width the numbers set. Text headers (CSV, tooltips) keep the one-word name. */
function scoreHeadHTML(label){ return String(label).replace(/^PlantCAD(\d?)$/, 'Plant<br>CAD$1'); }

/* Global helpers: the SNPCurate mark of a curated allele (window.SNP_CURATE, js/snpcurate.data.js)
   as a small link that opens its record. Each mark has its own glyph, so color is not the only
   cue: gold (validated causal change, genotyped here), outline (published marker, associated
   change or tag, genotyped here), grey (known allele this release cannot show). */
var CURATE_GLYPH = {gold:'\u2605', outline:'\u2606', grey:'\u25CB'};
function curateBadge(e){
  if (!e || !e.mark) return '';
  var tips = (window.SNP_CURATE && window.SNP_CURATE.marks) || {};
  var tt = String((tips[e.mark] || '') + ' ' + (e.symbol || '') + ' ' + (e.label || '') + ' (' + e.id + ')')
    .replace(/&/g,'&amp;').replace(/"/g,'&quot;').replace(/</g,'&lt;');
  return '<a href="#" class="cur-badge cur-' + e.mark + '" data-curate="' + e.id + '" data-tt="' + tt + '" aria-label="' + tt +
    '" onclick="event.stopPropagation();openCurate(\'' + e.id + '\');return false;">' + (CURATE_GLYPH[e.mark] || CURATE_GLYPH.grey) + '</a>';
}
function openCurate(id){
  if (typeof S !== 'undefined') S.curateId = id;
  if (typeof go === 'function') go('snpcurate');
}

/* Global helper: render a "Name (PFxxxxx)" domain string as a chip with the
   Pfam accession linked to InterPro. Used by SNPVersity / SNPImpact / SNPFunction. */
function pfamHref(pf){ return 'https://www.ebi.ac.uk/interpro/entry/pfam/' + pf + '/'; }
function domTag(dom){
  if (dom == null || dom === '\u2014' || dom === 'N/A' || dom === '')
    return '<span style="color:var(--faint)">\u2014</span>';
  var esc = function(s){ return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;'); };
  var m = String(dom).match(/(PF\d{4,6})/);
  if (m){
    var pf = m[1];
    var name = String(dom).replace(/\s*\(?\s*PF\d{4,6}\s*\)?\s*/, ' ').trim();
    return '<span class="dom-tag">' + esc(name) +
      ' <a href="' + pfamHref(pf) + '" target="_blank" rel="noopener" title="Pfam ' + pf +
      '" style="color:inherit;text-decoration:underline;text-decoration-style:dotted">' + pf + '</a></span>';
  }
  return '<span class="dom-tag">' + esc(dom) + '</span>';
}
