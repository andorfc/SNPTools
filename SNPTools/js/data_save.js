/* =====================================================================
 *  data.js — the DATA LAYER (LIVE / real HDF5 build).
 *
 *  queryVariants() no longer fabricates demo data: it POSTs the region +
 *  accession list to processForm.php, which runs h5_to_vcf.py against the
 *  real .h5 store, writes a VCF, and returns its path. We then fetch that
 *  VCF and parse it into the exact row shape the tools already expect.
 *
 *  Accession IDs come from accessions.real.js (window.SNP_REAL_ACCESSIONS),
 *  which holds the actual column names inside the .h5 files, so a selection
 *  maps to real HDF5 columns.
 *
 *  PFAM / protein-domain data is not in the HDF5 yet — the "Domain" column
 *  is filled with "N/A" until that separate data structure is built.
 *
 *  SNPImpact / SNPFold still use demo generators (they are not backed by
 *  these .h5 files); those bodies are left untouched below.
 *
 *  Depends on rnd() and pick() from core.js.
 * ===================================================================== */
const Data = (function () {

  /* =============================================================
   *  BACKEND CONFIG — edit these two lines if your paths differ.
   * ============================================================= */
  const CFG = {
    endpoint : 'processForm.php',  // relative to index.html
    vcfDir   : 'vcf/',             // web-served, writable; MUST match processForm.php
    geneEndpoint : 'lookupGeneModel.php',  // gene model -> coordinates
    // Regions wider than this fall back to a download link instead of a
    // full in-browser table (mirrors the original app's behavior).
    tableMaxSpan : 1_000_000,
  };

  /* ---------------- datasets (UI cards) ----------------
   * `id`     -> sent to processForm.php as dataSet
   * `family` -> which real accession list to show + which .h5 family
   */
  const DATASETS = [
    {id:'mgdb2026_hq', family:'mgdb2026',     name:'MaizeGDB 2026', sub:'High Quality',   ref:'B73 v5', acc:'2,710', sites:'98M',
     filters:['MQ ≥ 30','Coverage ≥ 50%','LD max R² > 0.5'], het:true,  indel:true,  impute:false},
    {id:'mgdb2026_hc', family:'mgdb2026',     name:'MaizeGDB 2026', sub:'High Coverage',  ref:'B73 v5', acc:'2,710', sites:'290M',
     filters:['MQ ≥ 30','Coverage ≥ 50%'], het:true,  indel:true,  impute:false},
    {id:'schnable2023',family:'schnable2023', name:'Schnable 2023', sub:'Imputed markers',ref:'B73 v5', acc:'1,515', sites:'12M',
     filters:['Imputed'], het:true, indel:true, impute:true},
    {id:'nam2021',     family:'nam2021',      name:'NAM 2021',      sub:'Founder panel',  ref:'B73 v5', acc:'27',    sites:'78M',
     filters:['MQ ≥ 30','Founder panel'], het:true,  indel:false,  impute:false},
    {id:'mgdb2024_hq', family:'mgdb2024',     name:'MaizeGDB 2024', sub:'High Quality',   ref:'B73 v5', acc:'1,498', sites:'83M',
     filters:['MQ ≥ 30','Coverage ≥ 50%'], het:true,  indel:true,  impute:false},
    {id:'mgdb2024_hc', family:'mgdb2024',     name:'MaizeGDB 2024', sub:'High Coverage',   ref:'B73 v5', acc:'1,498', sites:'228M',
      filters:['MQ ≥ 30','Coverage ≥ 50%','LD max R² > 0.5'], het:true,  indel:true,  impute:false},
  ];

  /* colors + a friendly label per family, for the accession picker group */
  const FAMILY_META = {
    mgdb2026:     {name:'MaizeGDB 2026', color:'#2563eb'},
    mgdb2024:     {name:'MaizeGDB 2024', color:'#cf8a12'},
    schnable2023: {name:'Schnable 2023', color:'#1f8a4c'},
    nam2021:      {name:'NAM 2021',      color:'#7c3aed'},
  };

  /* ---------------- accession catalog (projects -> groups -> accessions) ----------------
   * Source: window.SNP_CATALOG (compiled from data/accessions.tsv + data/projects.tsv).
   * Falls back to the older flat window.SNP_REAL_ACCESSIONS if present.
   */
  const CATALOG = (typeof window !== 'undefined' && window.SNP_CATALOG) || null;
  const REAL    = (typeof window !== 'undefined' && window.SNP_REAL_ACCESSIONS) || {};

  function familyOf(datasetId){
    const d = DATASETS.find(x => x.id === datasetId);
    return d ? d.family : 'mgdb2026';
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
        out.push({id:a.id, run:a.run, founder:a.founder, rep:a.rep, reps:a.reps,
                  label:a.label, group:g.name, namFounder:a.namFounder,
                  proj:p.id, projColor:p.color, projTitle:p.title});
      })));
    } else {
      out = (REAL[fam] || []).map(a => ({...a}));
    }
    _accCache[fam] = out;
    return out;
  }

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
    ['mgdb2026','mgdb2024','schnable2023','nam2021'].forEach(fam => {
      const dsid = (DATASETS.find(d => d.family === fam) || {}).id;
      if (dsid) accessionsFor(dsid).forEach(a => { if (!_byId.has(a.id)) _byId.set(a.id, a); });
    });
    return _byId;
  }
  function accessionById(id){ return idIndex().get(id) || null; }

  /* ---------------- genome geometry ---------------- */
  const GENE_MODELS = {
    'Zm00001eb374090':{chr:'chr8', start:163450112, end:163454880},
    'Zm00001eb067740':{chr:'chr2', start:21008440,  end:21013990},
    'Zm00001eb404760':{chr:'chr10',start:9821400,   end:9826110},
    'Zm00001eb404740':{chr:'chr10',start:9788220,   end:9794010},
    'Zm00001eb233650':{chr:'chr5', start:8841220,   end:8849510},
    'Zm00001eb313510':{chr:'chr7', start:174221000, end:174229800},
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
  async function lookupGene(id){
    id = (id || '').trim();
    if (!id) return null;
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
    return `${CFG.vcfDir}snpv_${ts}_${rnd}_${lo}_${hi}.vcf`;
  }

  // Pull one INFO tag out of the INFO column (values are ';'-separated).
  function info(infoStr, key){
    const m = infoStr.match(new RegExp('(?:^|;)' + key + '=([^;\\t]*)'));
    return m ? m[1] : null;
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
  function firstNum(v){                       // first of a possibly comma-listed value
    if (v == null) return null;
    return numOrNull(String(v).split(/[;,]+/)[0]);
  }

  function parseVcf(text, ids){
    const rows = [];
    if (!text) return {rows, sampleCols:{}, header:[]};
    const lines = text.split(/\r?\n/);
    let sampleCols = {};        // sampleName -> column index in the row
    for (const line of lines){
      if (!line) continue;
      if (line.startsWith('##')) continue;
      if (line.startsWith('#CHROM') || line.startsWith('#')){
        const cols = line.replace(/^#/, '').split('\t');
        for (let i = 9; i < cols.length; i++) sampleCols[cols[i]] = i;
        continue;
      }
      const t = line.split('\t');
      if (t.length < 8) continue;
      const infoStr = t[7] || '';

      // impact: take the most severe when a site lists several
      let impact = 'MODIFIER';
      const effTokens = (info(infoStr, 'EFFECT') || '').split(/[;,]+/);
      let best = -1;
      effTokens.forEach(e => {
        const key = e.trim().toUpperCase();
        if (key in SEVERITY && SEVERITY[key] > best){ best = SEVERITY[key]; impact = key; }
      });

      // genotypes in the exact order the caller selected
      const gts = ids.map(id => {
        const ci = sampleCols[id];
        if (ci == null || t[ci] == null) return './.';
        return (t[ci].split(':')[0] || './.').trim();
      });

      rows.push({
        pos:    parseInt(t[1], 10),
        ref:    t[3],
        alt:    t[4],
        gene:   cleanTok(info(infoStr, 'GENEMODEL')) || '—',
        effect: cleanTok(info(infoStr, 'TYPE')) || 'intergenic',
        impact,
        sub:    cleanTok(info(infoStr, 'SUB')),
        domain: 'N/A',                                  // PFAM data not in HDF5 yet
        mq:     (firstNum(info(infoStr, 'MQ')) != null) ? Math.round(firstNum(info(infoStr, 'MQ'))) : 'N/A',
        comp:   (firstNum(info(infoStr, 'CVP')) != null) ? firstNum(info(infoStr, 'CVP')) : 'N/A',
        r2:     firstNum(info(infoStr, 'MAXR2')),
        maf:    (firstNum(info(infoStr, 'MAF')) != null) ? firstNum(info(infoStr, 'MAF')) : 0,
        pc1:    numOrNull(info(infoStr, 'plantcad1_score')),
        pc2:    numOrNull(info(infoStr, 'plantcad2_score')),
        esm1:   numOrNull(info(infoStr, 'ESM1_score')),
        esm2:   numOrNull(info(infoStr, 'ESM2_score')),
        esm3:   numOrNull(info(infoStr, 'ESM3_score')),
        gts,
      });
    }
    return {rows, sampleCols};
  }

  // Build the accession objects the table header needs, in selection order.
  function accsFor(datasetId, ids){
    const map = new Map(accessionsFor(datasetId).map(a => [a.id, a]));
    return ids.map(id => map.get(id) || {id, run:id, founder:id, proj:familyOf(datasetId), projColor:'#8a94a6'});
  }

  /**
   * queryVariants(dataset, chr, lo, hi, ids) -> Promise<{rows, accs, chr, vcfUrl, span, wide, empty}>
   * `wide` is true when the interval exceeds tableMaxSpan (offer download instead of table).
   */
  async function queryVariants(dataset, chr, lo, hi, ids){
    const accs = accsFor(dataset, ids);
    const span = Math.max(hi - lo, 0);
    const outName = uniqueOutName(lo, hi);

    const body = new URLSearchParams({
      start: String(lo),
      end: String(hi),
      chr: chr,                         // e.g. "chr10" — matches the .h5 filename token
      dataSet: dataset,
      genotypes: JSON.stringify(ids),
      outName: outName,
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
      // Print the raw script output so the true cause is visible in the console.
      if (json.output)  console.error('[processForm.php] h5_to_vcf.py output:\n' + json.output);
      if (json.command) console.error('[processForm.php] command:\n' + json.command);
      const err = new Error(json.message || ('Request failed (HTTP ' + resp.status + ')'));
      err.detail = {command: json.command, output: json.output};
      throw err;
    }

    const vcfUrl = json.outFile || outName;

    // Very wide interval: don't try to render a giant table.
    if (span > CFG.tableMaxSpan){
      return {rows:[], accs, chr, vcfUrl, span, wide:true, empty:false};
    }

    const vcfResp = await fetch(vcfUrl, {cache:'no-store'});
    if (!vcfResp.ok){
      // No VCF written usually means the range returned no variants.
      return {rows:[], accs, chr, vcfUrl, span, wide:false, empty:true};
    }
    const {rows} = parseVcf(await vcfResp.text(), ids);
    return {rows, accs, chr, vcfUrl, span, wide:false, empty: rows.length === 0};
  }

  /* =============================================================
   *  SNPImpact query  (DEMO — not backed by these .h5 files)
   * ============================================================= */
  const BASES = ['A','C','G','T'];
  const CONSEQ = [
    {t:'Loss-of-function', cls:'lof',      base:-8.0},
    {t:'Loss-of-domain',  cls:'lod',      base:-6.0},
    {t:'Splice',          cls:'splice',   base:-3.0},
    {t:'Missense',        cls:'missense', base:-1.4},
    {t:'In-frame deletion',cls:'indel',   base:-0.6},
    {t:'Synonymous',      cls:'syn',      base: 0.3},
  ];
  const DOM_NAMES = ['Kinase domain','NB-ARC','bZIP','NAC domain','WRKY','DNA-binding domain','PPR repeat','F-box'];
  function priorityFromScore(s){ return s<=-7 ? 'TOP' : s<=-4 ? 'HIGH' : s<=-1 ? 'MODERATE' : 'LOW'; }

  function queryImpact(opts){
    opts = opts || {};
    const out = [];
    for (let i=0; i<46; i++){
      const c = pick(CONSEQ);
      const plantcad = +(c.base + rnd(-2,2)).toFixed(1);
      const esm      = +(c.base*0.72 + rnd(-1.5,1.5)).toFixed(1);
      const combined = +((plantcad + esm)/2).toFixed(2);
      const hasDom   = (c.cls==='lod' || c.cls==='missense' || Math.random()<.35);
      const aa = 90 + Math.floor(Math.random()*520);
      const exons = 4 + Math.floor(Math.random()*4);
      const affectedExon = 1 + Math.floor(Math.random()*exons);
      out.push({
        id:'v'+i,
        gene:'Zm00001eb'+(100000+Math.floor(Math.random()*899999)),
        variant: c.cls==='lof'      ? 'p.'+pick(['W','Q','R','E','K'])+aa+'*'
               : c.cls==='missense' ? 'p.'+pick(['A','G','R','D','V'])+aa+pick(['R','K','L','P','S'])
               : c.cls==='lod'      ? 'Δ Exon '+affectedExon
               : c.cls==='splice'   ? 'splice-site'
               : c.cls==='indel'    ? 'deletion'
               :                      'c.'+aa+pick(BASES)+'>'+pick(BASES),
        consequence:c.t, consClass:c.cls,
        domain: hasDom ? pick(DOM_NAMES) : '—',
        plantcad, esm, combined,
        priority: priorityFromScore(combined),
        percentile: Math.max(1, Math.min(99, Math.round(50 - combined*5 + rnd(-4,4)))),
        protLen: 280 + Math.floor(Math.random()*520),
        exons, affectedExon, aa,
      });
    }
    out.sort((a,b) => a.combined - b.combined);
    return out;
  }

  /* =============================================================
   *  SNPFold — protein structure + coding variants (DEMO curated).
   * ============================================================= */
  function structureFor(gene){ return (window.SNPFOLD_STRUCT||{})[gene] || null; }
  function pdbFor(gene){ return (window.SNPFOLD_PDB||{})[gene] || null; }

  const FOLD_VARIANTS = {
    'Zm00001eb378140':[
      {resi:31,  ref:'S', alt:'L', cons:'Missense',        cls:'missense', plantcad:-0.6, esm:-0.4},
      {resi:88,  ref:'P', alt:'S', cons:'Missense',        cls:'missense', plantcad:-1.1, esm:-0.9},
      {resi:145, ref:'A', alt:'V', cons:'Missense',        cls:'missense', plantcad:-2.8, esm:-2.1},
      {resi:206, ref:'G', alt:'R', cons:'Missense',        cls:'missense', plantcad:-7.9, esm:-6.2},
      {resi:268, ref:'D', alt:'N', cons:'Missense',        cls:'missense', plantcad:-6.4, esm:-5.1},
      {resi:300, ref:'R', alt:'C', cons:'Missense',        cls:'missense', plantcad:-3.6, esm:-3.0},
      {resi:350, ref:'T', alt:'M', cons:'Missense',        cls:'missense', plantcad:-1.4, esm:-1.0},
      {resi:441, ref:'W', alt:'*', cons:'Stop gained',     cls:'lof',      plantcad:-8.7, esm:-6.1},
      {resi:486, ref:'V', alt:'F', cons:'Missense',        cls:'missense', plantcad:-4.2, esm:-3.8},
      {resi:520, ref:'K', alt:'*', cons:'Frameshift',      cls:'lof',      plantcad:-7.2, esm:-5.4},
      {resi:555, ref:'E', alt:'K', cons:'Missense',        cls:'missense', plantcad:-2.2, esm:-1.7},
    ],
  };
  function queryFoldVariants(gene){
    const list = FOLD_VARIANTS[gene] || [];
    return list.map((v,i)=>{
      const combined = +((v.plantcad+v.esm)/2).toFixed(2);
      const variant = v.cls==='lof'
        ? (v.cons==='Frameshift' ? 'p.'+v.ref+v.resi+'fs' : 'p.'+v.ref+v.resi+'*')
        : 'p.'+v.ref+v.resi+v.alt;
      const ddg = v.cls==='lof' ? null : +(Math.abs(v.plantcad)*0.45 + rnd(-0.4,0.6)).toFixed(1);
      return {
        id:'f'+i, gene, resi:v.resi, ref:v.ref, alt:v.alt,
        variant, consequence:v.cons, consClass:v.cls,
        plantcad:v.plantcad, esm:v.esm, combined,
        priority: combined<=-7?'TOP':combined<=-4?'HIGH':combined<=-1?'MODERATE':'LOW',
        ddg,
        percentile: Math.max(1,Math.min(99, Math.round(50 - combined*5)))
      };
    });
  }

  /* ---------------- public API ---------------- */
  const DEFAULT_DS = DATASETS[0].id;
  return {
    datasets:    () => DATASETS,
    // dataset-aware accession accessors
    projectsFor, accessionsFor, defaultSelectionFor, familyOf, namFoundersFor,
    // backwards-compatible defaults (first dataset)
    projects:    () => projectsFor(DEFAULT_DS),
    accessions:  () => accessionsFor(DEFAULT_DS),
    defaultSelection: () => defaultSelectionFor(DEFAULT_DS),
    geneModels:  () => GENE_MODELS,
    exampleGenes:() => EXAMPLE_GENES,
    lookupGene,             // async: gene model id -> {id, chr, start, end} | null
    chromLengths:() => CHR_LEN,
    centromeres: () => CENTRO,
    accessionById,
    queryVariants,          // now async (returns a Promise)
    queryImpact,
    structureFor,
    pdbFor,
    queryFoldVariants,
  };
})();
