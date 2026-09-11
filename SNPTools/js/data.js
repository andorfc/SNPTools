/* =====================================================================
 *  data.js — the DATA LAYER (Fusarium, LIVE HDF5 build).
 *
 *  queryVariants() POSTs the region + accession list to processForm.php,
 *  which runs h5_to_vcf.py against the real Fusarium .h5 store (haploid),
 *  writes a VCF, and returns its path. We fetch + parse it into the row
 *  shape the tools expect.
 *
 *  Accession IDs are the real HDF5 column names:
 *    - F. graminearum : numeric isolate IDs, grouped by population (strains.js)
 *    - F. verticillioides 7600 / MRC826 : column names read from the HDF5
 *  All of this is compiled into window.SNP_CATALOG by build_strains_catalog.py.
 *
 *  Fusarium INFO carries MQ, CVC, CVP, TYPE, EFFECT, GENEMODEL, SUB, MAXR2,
 *  DNA_SCORE, AA_SCORE, MAF.  DNA_SCORE -> pc1, AA_SCORE -> esm1 (no secondary
 *  language-model scores yet).  Pfam domains / structures are not wired yet, so
 *  the Domain column degrades to '—'.
 *
 *  Depends on rnd() and pick() from core.js, and the global S (current query).
 * ===================================================================== */
const Data = (function () {

  const CFG = {
    endpoint     : 'processForm.php',
    vcfDir       : 'vcf/',
    geneEndpoint : 'lookupGeneModel.php',
    structDir    : 'data/structures/',
    // Pfam domains are per-reference (three separate assemblies with colliding
    // contig names — chr1/chr2/… mean different things in each genome), so the
    // domain stores are keyed by family, not shared.
    domainDirByFamily : {
      graminearum : 'data/domains_fusarium/Fgram_ph1/',
      vert7600    : 'data/domains_fusarium/Fvert_7600/',
      vertMRC826  : 'data/domains_fusarium/Fvert_mrc/',
    },
    // graminearum's by_chr uses source-GFF seqids: chr1-3 are named, but chr4
    // (and scaffolds) keep their RefSeq accession. Map the app's chr name to the
    // seqid the domain file is keyed under.
    domainChrAlias : {
      graminearum : { chr4: 'NC_026477.1' },
    },
    geneModelsDir  : 'data/genemodels/by_chr/',
    tableMaxSpan : 1_000_000,
  };

  /* ---------------- datasets (UI cards) ----------------
   * id     -> sent to processForm.php as dataSet
   * family -> which accession catalog + which .h5 family + which gene store
   */
  // Isolate + variant counts are the COMPLETE per-genome totals (all chromosomes),
  // from the 2026 site-summary tables. `sites` = SNP + INDEL total for that filter tier;
  // `snps`/`indels` break it down; `sitesShort` is the compact card label.
  const DATASETS = [
    {id:'gram_hq',     family:'graminearum', db:'graminearum', name:'F. graminearum 2026', sub:'High Quality',  ref:'PH-1 (FGSG)', acc:'512',
     sites:'120,390', sitesShort:'120K', snps:'111,754', indels:'8,636', chr:4, asm:'36.4 Mb',
     filters:['MQ ≥ 30','Coverage ≥ 50%','LD max R² > 0.5'], het:false, indel:true,  impute:false},
    {id:'gram_hc',     family:'graminearum', db:'graminearum', name:'F. graminearum 2026', sub:'High Coverage', ref:'PH-1 (FGSG)', acc:'512',
     sites:'2,485,599', sitesShort:'2.49M', snps:'2,339,005', indels:'146,594', chr:4, asm:'36.4 Mb',
     filters:['MQ ≥ 30','Coverage ≥ 50%'], het:false, indel:true,  impute:false},
    {id:'vert7600_hq', family:'vert7600',    db:'vert7600',   name:'F. verticillioides 7600', sub:'High Quality',  ref:'7600 (FVEG)', acc:'113',
     sites:'937,497', sitesShort:'937K', snps:'872,140', indels:'65,357', chr:11, asm:'41.1 Mb',
     filters:['MQ ≥ 30','Coverage ≥ 50%','LD max R² > 0.5'], het:false, indel:true,  impute:false},
    {id:'vert7600_hc', family:'vert7600',    db:'vert7600',   name:'F. verticillioides 7600', sub:'High Coverage', ref:'7600 (FVEG)', acc:'113',
     sites:'1,227,771', sitesShort:'1.23M', snps:'1,139,829', indels:'87,942', chr:11, asm:'41.1 Mb',
     filters:['MQ ≥ 30','Coverage ≥ 50%'], het:false, indel:true,  impute:false},
    {id:'vertmrc_hq',  family:'vertMRC826',  db:'vertMRC826', name:'F. verticillioides MRC826', sub:'High Quality',  ref:'MRC826 (FVERT4)', acc:'113',
     sites:'1,064,463', sitesShort:'1.06M', snps:'990,183', indels:'74,280', chr:12, asm:'42.9 Mb',
     filters:['MQ ≥ 30','Coverage ≥ 50%','LD max R² > 0.5'], het:false, indel:true,  impute:false},
    {id:'vertmrc_hc',  family:'vertMRC826',  db:'vertMRC826', name:'F. verticillioides MRC826', sub:'High Coverage', ref:'MRC826 (FVERT4)', acc:'113',
     sites:'1,402,131', sitesShort:'1.40M', snps:'1,304,907', indels:'97,224', chr:12, asm:'42.9 Mb',
     filters:['MQ ≥ 30','Coverage ≥ 50%'], het:false, indel:true,  impute:false},
  ];

  const FAMILY_META = {
    graminearum: {name:'F. graminearum 2026',      color:'#7c3aed'},
    vert7600:    {name:'F. verticillioides 7600',  color:'#2563eb'},
    vertMRC826:  {name:'F. verticillioides MRC826', color:'#0e7490'},
  };

  /* chromosome geometry per reference — EXACT lengths (bp) from the reference GFFs
     (region features / ##sequence-region). Used for the region ribbon; the query
     uses the entered start/end directly. */
  const CHR_LEN_BY_FAMILY = {
    // F. graminearum PH-1 — 4 chromosomes
    graminearum: {chr1:11697295, chr2:8914601, chr3:7713129, chr4:8033942},
    // F. verticillioides 7600 — 11 chromosomes
    vert7600: {chr1:6218892, chr2:4689494, chr3:4638610, chr4:4230616, chr5:4245961,
               chr6:3897820, chr7:3246105, chr8:2858895, chr9:2756060, chr10:2257929, chr11:2039996},
    // F. verticillioides MRC826 — 12 chromosomes
    vertMRC826: {chr1:6256590, chr2:4714944, chr3:4895763, chr4:4254065, chr5:4357113,
               chr6:4064809, chr7:3327862, chr8:2926371, chr9:2825676, chr10:2460917, chr11:2064835, chr12:693693},
  };

  /* example gene models per reference (verticillioides IDs resolve to matching
     chr tokens; graminearum FGSG_ IDs currently resolve to RefSeq scaffolds — see
     MISSING_DATA.md — so the graminearum autofill needs a chr remap to be useful). */
  // FungiDB-style (unpadded) IDs on purpose: the gene lookup normalizes zero-padding
  // (e.g. FVEG_03144 -> stored FVEG_003144), so the example buttons exercise the same
  // tolerance as a gene ID pasted straight from FungiDB. All resolve to placed (chr) genes.
  // Canonical (reference-GFF) display forms: FGSG 5-digit; FVEG / FVERT4 6-digit.
  const EXAMPLE_GENES_BY_FAMILY = {
    graminearum: ['FGSG_00777','FGSG_00778','FGSG_03537','FGSG_10375'],
    vert7600:    ['FVEG_003144','FVEG_000765','FVEG_014423','FVEG_005000'],
    vertMRC826:  ['FVERT4_000001','FVERT4_000765','FVERT4_005000','FVERT4_010000'],
  };
  const CATALOG = (typeof window !== 'undefined' && window.SNP_CATALOG) || null;
  const REAL    = (typeof window !== 'undefined' && window.SNP_REAL_ACCESSIONS) || {};

  function familyOf(datasetId){
    const d = DATASETS.find(x => x.id === datasetId);
    return d ? d.family : 'graminearum';
  }
  function dbOf(datasetId){
    const d = DATASETS.find(x => x.id === datasetId);
    return d ? d.db : 'graminearum';
  }
  /* Language-model score columns to show for a dataset, in order. As of the 2026
     rebuild ALL three references carry the same six models — two DNA (FunDLM, EVO2)
     and four protein (ESM1, ESM2, ESM3, ESM-C) — in the INFO
     (FUNDLM/EVO2/ESM1/ESM2/ESM3/ESMC_SCORE). `key` indexes the parsed row
     (pc1/pc2 = DNA, esm1..esm3/esmc = protein). */
  function scoreModels(datasetId){
    return [
      {key:'pc1',  kind:'dna',     label:'FunDLM', tip:'FunDLM DNA language-model score (INFO FUNDLM_SCORE); more extreme = more disruptive.'},
      {key:'pc2',  kind:'dna',     label:'EVO2',   tip:'EVO2 DNA language-model score (INFO EVO2_SCORE).'},
      {key:'esm1', kind:'protein', label:'ESM1',   tip:'ESM1 protein language-model score (INFO ESM1_SCORE).'},
      {key:'esm2', kind:'protein', label:'ESM2',   tip:'ESM2 protein language-model score (INFO ESM2_SCORE).'},
      {key:'esm3', kind:'protein', label:'ESM3',   tip:'ESM3 protein language-model score (INFO ESM3_SCORE).'},
      {key:'esmc', kind:'protein', label:'ESM-C',  tip:'ESM-C protein language-model score (INFO ESMC_SCORE).'},
    ];
  }
  // All three references now carry the protein (secondary) models.
  function hasSecondaryScores(datasetId){
    return ['graminearum','vert7600','vertMRC826'].indexOf(familyOf(datasetId)) >= 0;
  }
  function famNode(datasetId){ return CATALOG ? CATALOG.families[familyOf(datasetId)] : null; }

  function projectsFor(datasetId){
    const fam = famNode(datasetId);
    if (fam) return fam.projects;
    const list = REAL[familyOf(datasetId)] || [];
    return [{id:familyOf(datasetId), title:familyOf(datasetId), bioprojects:[], color:'#2563eb',
             count:list.length, namFounders:[], groups:[{name:null, accessions:list}]}];
  }

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
                  // graminearum metadata
                  population:a.population, species:a.species, host:a.host,
                  chemotype:a.chemotype, country:a.country,
                  // verticillioides metadata (from Fvert_MetaData xlsx)
                  strain:a.strain, region:a.region, substrate:a.substrate,
                  substrateDetail:a.substrateDetail, geo:a.geo,
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

  function defaultSelectionFor(datasetId){
    const list = accessionsFor(datasetId);
    const ids = [], seen = new Set();
    for (const a of list){
      if (seen.has(a.founder)) continue;
      seen.add(a.founder); ids.push(a.id);
      if (ids.length >= 12) break;
    }
    return ids;
  }

  let _byId = null;
  function idIndex(){
    if (_byId) return _byId;
    _byId = new Map();
    ['graminearum','vert7600','vertMRC826'].forEach(fam => {
      const dsid = (DATASETS.find(d => d.family === fam) || {}).id;
      if (dsid) accessionsFor(dsid).forEach(a => { if (!_byId.has(a.id)) _byId.set(a.id, a); });
    });
    return _byId;
  }
  function accessionById(id){ return idIndex().get(id) || null; }

  /* ---------------- chromosome geometry (family-aware) ---------------- */
  function chromLengths(datasetId){
    return CHR_LEN_BY_FAMILY[familyOf(datasetId)] || CHR_LEN_BY_FAMILY.graminearum;
  }
  const CENTRO = {};   // Fusarium centromere positions not wired yet

  /* ---------------- gene model -> coordinates (live) ----------------
   * Resolves a gene-model ID to its interval via lookupGeneModel.php, choosing
   * the serialized GFF store from the current dataset's reference (database=).
   */
  function exampleGenes(datasetId){
    return EXAMPLE_GENES_BY_FAMILY[familyOf(datasetId)] || EXAMPLE_GENES_BY_FAMILY.graminearum;
  }
  async function lookupGene(id, database){
    id = (id || '').trim();
    if (!id) return null;
    // default the reference store from the active dataset when not passed explicitly
    if (!database){
      try { database = dbOf((typeof S !== 'undefined' && S.dataset) ? S.dataset : DATASETS[0].id); }
      catch(e){ database = 'graminearum'; }
    }
    const url = `${CFG.geneEndpoint}?geneModelId=${encodeURIComponent(id)}&database=${encodeURIComponent(database)}`;
    const resp = await fetch(url, {cache:'no-store'});
    if (!resp.ok) throw new Error('Gene lookup failed (HTTP ' + resp.status + ')');
    const raw = await resp.text();
    let d;
    try { d = JSON.parse(raw); }
    catch (e) { throw new Error('lookupGeneModel.php did not return JSON:\n' + raw.slice(0, 600)); }
    if (!d || d.id === 'empty' || d.chromosome == null) return null;
    const start = parseInt(d.start, 10), end = parseInt(d.end, 10);
    const isChr  = /^chr[0-9]+$/.test(String(d.chromosome));
    const placed = (d.placed != null) ? !!d.placed : isChr;
    if (!placed){
      // The gene exists but sits on an unplaced scaffold / mitochondrion, not on a
      // chr1..chrN chromosome — so it has no region in the variant store. Return a
      // flagged object (chr:null) so callers can show a clear message instead of
      // trying to query an invalid chromosome.
      return {id, chr:null, placed:false, scaffold:(d.scaffold || d.chromosome),
              start:Number.isFinite(start)?start:null, end:Number.isFinite(end)?end:null};
    }
    if (!Number.isFinite(start) || !Number.isFinite(end) || (start === 0 && end === 0)) return null;
    return {id, chr: d.chromosome, start, end, placed:true};
  }

  /* =============================================================
   *  SNPVERSITY — LIVE query.  region + accession ids -> VCF -> rows
   * ============================================================= */
  const SEVERITY = {HIGH:3, MODERATE:2, LOW:1, MODIFIER:0};

  function uniqueOutName(lo, hi){
    const ts = Date.now();
    const rnd = Math.random().toString().slice(2, 11);
    return `${CFG.vcfDir}snpv_${ts}_${rnd}_${lo}_${hi}.vcf`;
  }

  function info(infoStr, key){
    const m = infoStr.match(new RegExp('(?:^|;)' + key + '=([^;\\t]*)'));
    return m ? m[1] : null;
  }
  const cleanTok = s => {
    if (s == null) return '';
    const tok = String(s).split(/[;,]+/)[0].replace(/_/g, ' ').trim();
    return tok === '.' ? '' : tok;
  };
  // Gene-model IDs must keep their underscores (FVEG_000765, FGSG_00777, FVERT4_010000)
  // — cleanTok's _->space rewrite would mangle them and break the gene->annotation join.
  // Take the first listed gene; leave boundary spans ("A..B") intact for isSingleGeneModel().
  const cleanGene = s => {
    if (s == null) return '';
    const tok = String(s).split(/[;,]+/)[0].trim();
    return tok === '.' ? '' : tok;
  };
  /* Canonical gene id: prefix + zero-padding-stripped number, uppercased
     (FVEG_003144 ⇔ FVEG_03144 ⇔ FVEG_3144). Lets a gene id compare equal to the
     VCF's GENEMODEL regardless of how each zero-pads the numeric part. */
  // NOTE: the prefix can itself contain digits (FVERT4_), so the prefix group must be
  // ".*_" (up to the final underscore), NOT [A-Za-z]+ — otherwise FVERT4 ids fail to
  // match and silently skip padding normalization.
  function geneCanon(g){
    const m = String(g == null ? '' : g).trim().match(/^(.*_)0*(\d+)$/);
    return m ? (m[1].toUpperCase() + m[2]) : String(g == null ? '' : g).trim().toUpperCase();
  }
  /* Canonical DISPLAY id — the FungiDB/reference-GFF form. The reference pads the
     numeric part to a fixed width per genome, EXCEPT a high-numbered tail block that
     is left at its natural (unpadded) length. Partition taken from the reference GFFs
     (2026-08, no numeric value straddles two widths):
       graminearum (FGSG_)   : 5-digit across the whole range (matches the VCF already)
       vert7600    (FVEG_)    : 6-digit up to 15512, natural for 15513..15545
       vertMRC826  (FVERT4_)  : 6-digit up to 15520, natural for 15521..15581
     The variant store (SNPVersity) writes these 5-digit, so applying this at parse time
     makes SNPVersity + every hand-off show the same 6-digit id as SNPFold / SNPFunction /
     SNPImpact. Matching stays padding-tolerant (geneCanon), so a differently-padded id
     still resolves. */
  const CANON_ID = {
    graminearum: { width: 5, natFrom: Infinity },
    vert7600:    { width: 6, natFrom: 15513 },
    vertMRC826:  { width: 6, natFrom: 15521 },
  };
  function familyFromGeneId(g){
    const s = String(g == null ? '' : g).toUpperCase();
    if (s.indexOf('FGSG')   === 0) return 'graminearum';
    if (s.indexOf('FVERT4') === 0) return 'vertMRC826';   // test before FVEG (FVERT4 ≠ FVEG)
    if (s.indexOf('FVEG')   === 0) return 'vert7600';
    return null;
  }
  function canonicalGeneId(gene){
    const s = String(gene == null ? '' : gene).trim();
    const m = s.match(/^(.*_)0*(\d+)$/);              // prefix incl. any digits (FVERT4_) up to final "_"
    if (!m) return s;                                 // compound/boundary/non-numeric: leave as-is
    const n = parseInt(m[2], 10);
    const cfg = CANON_ID[familyFromGeneId(s)];
    if (!cfg) return m[1] + String(n);                // unknown prefix: just strip leading zeros
    const num = (n >= cfg.natFrom) ? String(n) : String(n).padStart(cfg.width, '0');
    return m[1] + num;
  }
  /* A variant belongs to a gene when the gene is its PRIMARY annotation — the
     first GENEMODEL token, already stored as r.gene — compared padding-tolerantly.
     Deliberately NOT matching secondary comma-members ("FGSG_00778,FGSG_00777"):
     the residue we read (r.sub's first value) and the effect belong to that
     primary gene, so pulling the variant into an overlapping neighbor's list would
     attach the WRONG residue (appearing beyond the neighbor's protein/structure)
     and surface variants under a gene that SNPVersity assigns elsewhere. */
  function geneRowMatches(r, gene){
    return geneCanon(r && r.gene) === geneCanon(gene);
  }
  function numOrNull(v){
    if (v == null || v === '.' || v === '') return null;
    const f = parseFloat(v);
    return Number.isNaN(f) ? null : f;
  }
  function firstNum(v){
    if (v == null) return null;
    return numOrNull(String(v).split(/[;,]+/)[0]);
  }

  function parseVcf(text, ids, family){
    const rows = [];
    if (!text) return {rows, sampleCols:{}, header:[]};
    const lines = text.split(/\r?\n/);
    let sampleCols = {};
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

      let impact = 'MODIFIER';
      const effTokens = (info(infoStr, 'EFFECT') || '').split(/[;,]+/);
      let best = -1;
      effTokens.forEach(e => {
        const key = e.trim().toUpperCase();
        if (key in SEVERITY && SEVERITY[key] > best){ best = SEVERITY[key]; impact = key; }
      });

      // Fusarium is haploid: a sample column is a single allele index ("0","1",".").
      const gts = ids.map(id => {
        const ci = sampleCols[id];
        if (ci == null || t[ci] == null) return '.';
        return (t[ci].split(':')[0] || '.').trim();
      });

      // INDEL detection: REF/ALT differ in length (insertion/deletion). The DNA
      // substitution language models (FunDLM / DNABERT, and Evo2) are point-mutation
      // scorers — their values for indels are artifacts, so treat them as missing,
      // the same way the protein models already are for indels. Protein models
      // (ESM*) are '.' for indels in the store, so no extra guard is needed there.
      const isIndel = String(t[3] || '').length !== String(t[4] || '').length;
      rows.push({
        pos:    parseInt(t[1], 10),
        ref:    t[3],
        alt:    t[4],
        gene:   canonicalGeneId(cleanGene(info(infoStr, 'GENEMODEL'))) || '—',
        geneRaw: info(infoStr, 'GENEMODEL') || '',
        effect: cleanTok(info(infoStr, 'TYPE')) || 'intergenic',
        impact,
        sub:    cleanTok(info(infoStr, 'SUB')),
        domain: domainAt(family, t[0], parseInt(t[1], 10)),
        mq:     (firstNum(info(infoStr, 'MQ')) != null) ? Math.round(firstNum(info(infoStr, 'MQ'))) : 'N/A',
        comp:   (firstNum(info(infoStr, 'CVP')) != null) ? firstNum(info(infoStr, 'CVP')) : 'N/A',
        r2:     firstNum(info(infoStr, 'MAXR2')),
        maf:    (firstNum(info(infoStr, 'MAF')) != null) ? firstNum(info(infoStr, 'MAF')) : 0,
        // Language-model scores. 2026 format (graminearum): two DNA models
        // FUNDLM_SCORE/EVO2_SCORE + four protein models ESM1/2/3/ESMC_SCORE.
        // 2025 format (verticillioides, until its rebuild): DNA_SCORE (DNABERT) +
        // AA_SCORE (ESM1). Read the new names first, fall back to the old ones so
        // both formats parse. Column keys: pc1/pc2 = DNA, esm1..esm3/esmc = protein.
        pc1:    isIndel ? null : (numOrNull(info(infoStr, 'FUNDLM_SCORE')) != null ? numOrNull(info(infoStr, 'FUNDLM_SCORE')) : numOrNull(info(infoStr, 'DNA_SCORE'))),
        pc2:    isIndel ? null : numOrNull(info(infoStr, 'EVO2_SCORE')),
        esm1:   numOrNull(info(infoStr, 'ESM1_SCORE')) != null ? numOrNull(info(infoStr, 'ESM1_SCORE')) : numOrNull(info(infoStr, 'AA_SCORE')),
        esm2:   numOrNull(info(infoStr, 'ESM2_SCORE')),
        esm3:   numOrNull(info(infoStr, 'ESM3_SCORE')),
        esmc:   numOrNull(info(infoStr, 'ESMC_SCORE')),
        gts,
      });
    }
    return {rows, sampleCols};
  }

  function accsFor(datasetId, ids){
    const map = new Map(accessionsFor(datasetId).map(a => [a.id, a]));
    return ids.map(id => map.get(id) || {id, run:id, founder:id, proj:familyOf(datasetId), projColor:'#8a94a6'});
  }

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

  async function queryVariants(dataset, chr, lo, hi, ids){
    ids = sortIds(dataset, ids);
    const accs = accsFor(dataset, ids);
    const span = Math.max(hi - lo, 0);
    const outName = uniqueOutName(lo, hi);

    const body = new URLSearchParams({
      start: String(lo), end: String(hi), chr: chr,
      dataSet: dataset, genotypes: JSON.stringify(ids), outName: outName,
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
    if (span > CFG.tableMaxSpan){
      return {rows:[], accs, chr, vcfUrl, span, wide:true, empty:false};
    }
    const vcfResp = await fetch(vcfUrl, {cache:'no-store'});
    if (!vcfResp.ok){
      return {rows:[], accs, chr, vcfUrl, span, wide:false, empty:true};
    }
    const vcfText = await vcfResp.text();
    const fam = familyOf(dataset);
    await ensureDomains(fam);
    const {rows} = parseVcf(vcfText, ids, fam);
    return {rows, accs, chr, vcfUrl, span, wide:false, empty: rows.length === 0};
  }

  /* SNPGeo gene lookup: given a gene name, look up its coordinates and query variants in that region */
  async function queryVariantsByGene(dataset, geneId, ids){
    dataset = dataset || DATASETS[0].id;
    const db = dbOf(dataset);
    
    // Look up gene coordinates
    const gene = await lookupGene(geneId, db);
    if (!gene){
      throw new Error('Gene not found: ' + geneId);
    }
    if (gene.chr == null){
      // Gene is on unplaced scaffold, not in variant store
      throw new Error('Gene ' + geneId + ' is on scaffold "' + gene.scaffold + '", not in variant database');
    }
    
    // Use default accession selection if none provided
    ids = ids || defaultSelectionFor(dataset);
    
    // Query variants in the gene region
    const result = await queryVariants(dataset, gene.chr, gene.start, gene.end, ids);
    
    // Return with gene context
    return {
      rows: result.rows,
      accs: result.accs,
      chr: gene.chr,
      start: gene.start,
      end: gene.end,
      gene: geneId,
      vcfUrl: result.vcfUrl,
      span: result.span,
      wide: result.wide,
      empty: result.empty
    };
  }

  /* SNPImpact ranks REAL variants: the tool queries a region via queryVariants()
     (or receives a SNPVersity hand-off) and calls rankImpact(rows) below. The old
     synthetic queryImpact() generator was removed once the real path landed. */

  /* =============================================================
   *  SNPFold — protein structure + coding variants
   * ============================================================= */
  function structureFor(gene){ return (window.SNPFOLD_STRUCT||{})[gene] || null; }
  function pdbFor(gene){ return (window.SNPFOLD_PDB||{})[gene] || null; }
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

  /* ---- Pfam domains by genomic position — one store per reference family ----
   * Each family's domains.by_chr.json is {contig: [[gstart,gend,name,pfam,type,
   * tx,instance], …]} sorted by gstart (domains crossing introns span several
   * rows). Loaded whole, once per family, and cached. domainAt() takes the
   * family so the three colliding contig namespaces stay separate. */
  const _domByFam = {};        // family -> {contig:[rows...]}
  const _domFamProm = {};
  function domainFolder(family){
    return CFG.domainDirByFamily[family] || CFG.domainDirByFamily.graminearum;
  }
  function domainChrKey(family, chr){
    const al = (CFG.domainChrAlias || {})[family];
    return (al && al[chr]) || chr;
  }
  function ensureDomains(family){
    if (family in _domByFam) return Promise.resolve(_domByFam[family]);
    if (_domFamProm[family]) return _domFamProm[family];
    _domFamProm[family] = fetch(domainFolder(family) + 'domains.by_chr.json', {cache:'force-cache'})
      .then(r => r.ok ? r.json() : {}).catch(() => ({}))
      .then(idx => { _domByFam[family] = idx || {}; return _domByFam[family]; });
    return _domFamProm[family];
  }
  function domainAt(family, chr, pos){
    const idx = _domByFam[family];
    if (!idx) return '—';
    const a = idx[domainChrKey(family, chr)];
    if (!a || !a.length) return '—';
    let lo = 0, hi = a.length;
    while (lo < hi){ const m = (lo + hi) >> 1; if (a[m][0] <= pos) lo = m + 1; else hi = m; }
    let best = null;
    for (let i = lo - 1; i >= 0; i--){
      const iv = a[i];
      if (iv[0] < pos - 100000) break;          // domain blocks are short; safe cutoff
      if (iv[0] <= pos && iv[1] >= pos){
        if (!best || (iv[1] - iv[0]) < (best[1] - best[0])) best = iv;   // tightest wins
      }
    }
    return best ? `${best[2]}${best[3] ? ` (${best[3]})` : ''}` : '—';
  }

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
    return 'p.' + (p.ref || '') + p.resi + (p.alt || '');
  }

  async function queryFoldVariants(gene, dataset, ids){
    dataset = dataset || DATASETS[0].id;
    const g = await lookupGene(gene, dbOf(dataset));
    if (!g || !g.chr) return [];
    ids = ids || defaultSelectionFor(dataset);
    let res;
    try { res = await queryVariants(dataset, g.chr, g.start, g.end, ids); }
    catch (e){ console.warn('queryFoldVariants:', e && e.message); return []; }
    const out = [];
    for (const r of (res.rows || [])){
      if (gene && !geneRowMatches(r, gene) && r.gene !== '—') continue;   // padding-/compound-tolerant
      const cls = classifyConsequence(r.effect);
      if (!cls) continue;
      const p = parseSub(r.sub);
      if (!p || p.resi == null) continue;
      const combined = compositeScore(r);
      out.push({
        id:'f' + out.length, gene, resi:p.resi, ref:p.ref, alt:p.alt,
        variant: hgvsProtein(p, cls), consequence: cls.label, consClass: cls.klass,
        structural: cls.structural, pos:r.pos, refNt:r.ref, altNt:r.alt,
        impact: r.impact || null, maf: (r.maf != null ? r.maf : null),
        ...lmScores(r), combined,
        priority: combined == null ? null
                : (combined <= -7 ? 'TOP' : combined <= -4 ? 'HIGH' : combined <= -1 ? 'MODERATE' : 'LOW'),
      });
    }
    out.sort((a, b) => a.resi - b.resi);
    return out;
  }

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
    if (/5_prime_utr|five_prime/.test(e))        return {klass:'other',  label:'5′ UTR',        severe:false};
    if (/3_prime_utr|three_prime/.test(e))       return {klass:'other',  label:'3′ UTR',        severe:false};
    if (/upstream/.test(e))                      return {klass:'other',  label:'Upstream',           severe:false};
    if (/downstream/.test(e))                    return {klass:'other',  label:'Downstream',         severe:false};
    if (/intergenic/.test(e))                    return {klass:'other',  label:'Intergenic',         severe:false};
    return {klass:'other', label: (cleanTok(effect) || 'Other'), severe:false};
  }
  function impactPriority(cls, combined, level){
    if (cls.severe) return 'TOP';
    if (cls.klass === 'missense' || cls.klass === 'indel'){
      if (combined != null){
        if (combined <= -7) return 'TOP';
        if (combined <= -4) return 'HIGH';
        if (combined <= -1) return 'MODERATE';
        return 'LOW';
      }
      return level === 'HIGH' ? 'HIGH' : level === 'MODERATE' ? 'MODERATE' : 'LOW';
    }
    return level === 'HIGH' ? 'HIGH' : 'LOW';
  }
  function rankImpact(rows){
    const out = [];
    for (const r of (rows || [])){
      const cls = impactClass(r.effect);
      const combined = compositeScore(r);
      const p = parseSub(r.sub);
      const coding = (cls.klass === 'missense' || cls.klass === 'lof' || cls.klass === 'indel');
      const variant = (p && p.resi != null && coding)
        ? hgvsProtein(p, {label: cls.label, klass: cls.klass})
        : `${r.pos} ${r.ref}>${r.alt}`;
      out.push({
        id: 'i' + out.length, gene: r.gene, pos: r.pos, ref: r.ref, alt: r.alt,
        variant, consequence: cls.label, consClass: cls.klass,
        resi: p ? p.resi : null, aaRef: p ? p.ref : null, aaAlt: p ? p.alt : null,
        domain: r.domain || '—', impactLevel: r.impact || null,
        ...lmScores(r), combined,
        maf: r.maf, r2: r.r2, mq: r.mq,
        priority: impactPriority(cls, combined, r.impact), percentile: null,
      });
    }
    const scored = out.filter(v => v.combined != null).slice().sort((a, b) => a.combined - b.combined);
    const N = scored.length;
    scored.forEach((v, i) => { v.percentile = N ? Math.round(100 * (N - i) / N) : null; });
    return out;
  }

  /* Per-gene protein-domain index (SNPFold / SNPImpact tracks). The three
   * references live in separate files but gene IDs are prefix-unique
   * (FGSG_/FVEG_/FVERT4_), so all three merge into one lookup keyed by gene. */
  let _geneDom = null, _geneDomProm = null;
  function ensureGeneDomains(){
    if (_geneDom) return Promise.resolve(_geneDom);
    if (_geneDomProm) return _geneDomProm;
    const fams = Object.keys(CFG.domainDirByFamily);
    _geneDomProm = Promise.all(fams.map(f =>
      fetch(domainFolder(f) + 'domains.by_gene.json', {cache:'force-cache'})
        .then(r => r.ok ? r.json() : {}).catch(() => ({}))
    )).then(parts => {
      const merged = {};
      parts.forEach(p => { for (const k in p) if (!(k in merged)) merged[k] = p[k]; });
      _geneDom = merged; return _geneDom;
    });
    return _geneDomProm;
  }
  function _geneDomCandidates(gene){
    const g = String(gene || '').trim();
    if (!g) return [];
    const out = [g];
    const m = g.match(/^(.*_)(\d+)$/);   // pad/unpad tolerance (FVEG_765 vs FVEG_000765; prefix may contain digits, e.g. FVERT4_)
    if (m){
      const pre = m[1], num = m[2].replace(/^0+/, '') || '0';
      for (const w of [5, 6, 7]) out.push(pre + num.padStart(w, '0'));
      out.push(pre + num);
    }
    return out;
  }
  function geneDomains(gene){
    const map = _geneDom || {};
    for (const c of _geneDomCandidates(gene)){ if (map[c]) return map[c]; }
    return null;
  }

  const _gmByChr = {}, _gmProm = {};
  function ensureGeneModels(chr){
    if (chr in _gmByChr) return Promise.resolve(_gmByChr[chr]);
    if (_gmProm[chr]) return _gmProm[chr];
    _gmProm[chr] = fetch(CFG.geneModelsDir + encodeURIComponent(chr) + '.json', {cache:'force-cache'})
      .then(r => r.ok ? r.json() : {}).catch(() => ({}))
      .then(x => { _gmByChr[chr] = x || {}; return _gmByChr[chr]; });
    return _gmProm[chr];
  }
  /* A variant's GENEMODEL can be a single gene ("FGSG_00780"), a comma list of
   * overlapping genes ("FGSG_00778,FGSG_00777"), or a boundary span for an
   * intergenic/up-downstream variant ("FGSG_00780..FGSG_00781"). Split all of
   * those into individual gene ids so the gene-model view can resolve one. */
  function _geneTokens(geneStr){
    const s = String(geneStr == null ? '' : geneStr).trim();
    if (!s || s === '.') return [];
    const out = [];
    s.split(/[;,]+/).forEach(part => {
      part = part.trim();
      if (!part || part === '.') return;
      if (part.indexOf('..') !== -1) part.split('..').forEach(p => { p = p.trim(); if (p && p !== '.') out.push(p); });
      else out.push(part);
    });
    return out;
  }
  function geneModelOf(chr, gene, pos){
    const map = _gmByChr[chr] || {};
    const hits = []; const seen = new Set();
    for (const tok of _geneTokens(gene)){
      for (const c of _geneDomCandidates(tok)){               // padding-tolerant
        if (map[c] && !seen.has(c)){ seen.add(c); hits.push(map[c]); break; }
      }
    }
    if (!hits.length) return null;
    if (pos == null || hits.length === 1) return hits[0];
    // multiple candidates (intergenic between genes): prefer the gene containing
    // the variant, else the nearest — so the marker lands in the right flank.
    let best = hits[0], bestD = Infinity;
    for (const m of hits){
      const lo = Math.min.apply(null, m.exons.map(e => e[0]));
      const hi = Math.max.apply(null, m.exons.map(e => e[1]));
      if (pos >= lo && pos <= hi) return m;
      const d = pos < lo ? lo - pos : pos - hi;
      if (d < bestD){ bestD = d; best = m; }
    }
    return best;
  }

  function _dose(g){
    if (g == null) return null;
    const s = String(g); if (s==='./.'||s==='.'||s==='') return null;
    // Fusarium is HAPLOID: a single-allele genotype is the isolate's only copy, so
    // an alt call is a FULL-dosage carrier (dose 2 = "homozygous"/knockout-eligible),
    // never heterozygous. (Matches the dosage convention in SNPTree/Matrix/Compare.)
    if (!/[\/|]/.test(s)) return s === '0' ? 0 : 2;
    const a = s.split(/[\/|]/); if (a.length < 2 || a[0]==='.' || a[1]==='.') return null;
    return (a[0]!=='0'?1:0) + (a[1]!=='0'?1:0);
  }
  function _avg(a){ return a.length ? +(a.reduce((s,x)=>s+x,0)/a.length).toFixed(2) : null; }

  /* ---- language-model scores on a variant ----
     Carry all six model scores (pc1/pc2 = DNA: FunDLM/EVO2; esm1/esm2/esm3/esmc =
     protein: ESM1/ESM2/ESM3/ESM-C) from a parsed row onto a variant object, keyed
     exactly like Data.scoreModels() so tools can render `v[model.key]` directly.
     `plantcad`/`esm` are kept as legacy aliases (= primary DNA / primary protein). */
  function lmScores(r){
    return { pc1:r.pc1, pc2:r.pc2, esm1:r.esm1, esm2:r.esm2, esm3:r.esm3, esmc:r.esmc,
             plantcad:r.pc1, esm:r.esm1, plantcad2:r.pc2 };
  }
  function _meanOf(s, keys){
    const a = keys.map(k => s[k]).filter(x => x != null);
    return a.length ? +(a.reduce((x,y)=>x+y,0)/a.length).toFixed(2) : null;
  }
  /* Composite = consensus of the PROTEIN models (ESM1/2/3/ESM-C) — a comparable
     scale and the direct amino-acid-impact signal — falling back to the DNA models
     (FunDLM/EVO2) for non-coding sites. Drives priority tiers, percentile, and the
     severity labels. (Old builds averaged one DNA + one protein score.) */
  function compositeScore(s){
    const prot = _meanOf(s, ['esm1','esm2','esm3','esmc']);
    return prot != null ? prot : _meanOf(s, ['pc1','pc2']);
  }
  function meanScoresOf(variants){
    const keys = ['pc1','pc2','esm1','esm2','esm3','esmc'];
    const out = {};
    keys.forEach(k => { out[k] = _avg(variants.map(v => v[k]).filter(x => x != null)); });
    return out;
  }

  async function geneFunction(gene, dataset){
    dataset = dataset || DATASETS[0].id;
    const dsName = (DATASETS.find(d=>d.id===dataset)||{}).name || dataset;
    const g = await lookupGene(gene, dbOf(dataset));
    if (!g || !g.chr) return {gene, dataset, datasetName:dsName, error:'No gene-model coordinates found for '+gene+'.'};
    const ids = (accessionsFor(dataset)||[]).map(a=>a.id);
    await Promise.all([ensureDomains(familyOf(dataset)), ensureGeneDomains(), ensureGeneModels(g.chr)]);
    let res;
    try { res = await queryVariants(dataset, g.chr, g.start, g.end, ids); }
    catch (e){ return {gene, chr:g.chr, start:g.start, end:g.end, dataset, datasetName:dsName, error:'Variant query failed: '+(e&&e.message)}; }

    const accs = res.accs || [];
    const rows = (res.rows || []).filter(r => geneRowMatches(r, gene));   // padding-/compound-tolerant
    const gd = geneDomains(gene), gm = geneModelOf(g.chr, gene);

    const variants = rows.map((r, vi) => {
      const cls = impactClass(r.effect);
      const combined = compositeScore(r);
      let het=0, hom=0, called=0; const homIds=[], hetIds=[];
      for (let k=0;k<accs.length;k++){ const d=_dose(r.gts[k]);
        if (d==null) continue; called++;
        if (d===2){ hom++; homIds.push(accs[k].id); } else if (d===1){ het++; hetIds.push(accs[k].id); } }
      const an = called, ac = het + hom, af = an ? ac/an : 0;   // haploid: allele count == carrier count
      const p = parseSub(r.sub);
      const coding = (cls.klass==='missense'||cls.klass==='lof'||cls.klass==='indel');
      const variant = (p && p.resi!=null && coding) ? hgvsProtein(p, {label:cls.label, klass:cls.klass}) : `${r.pos} ${r.ref}>${r.alt}`;
      return {id:'fx'+vi, pos:r.pos, ref:r.ref, alt:r.alt, variant, consequence:cls.label, consClass:cls.klass, severe:cls.severe,
        domain:r.domain||'—', resi:p?p.resi:null, aaRef:p?p.ref:null, aaAlt:p?p.alt:null, ...lmScores(r), combined,
        priority: impactPriority(cls, combined, r.impact),
        het, hom, af, carriersHom:homIds, carriersHet:hetIds};
    });

    const byClass = {missense:0, lof:0, splice:0, indel:0, syn:0, other:0};
    variants.forEach(v => { byClass[v.consClass] = (byClass[v.consClass]||0) + 1; });
    const nonsyn = byClass.missense + byClass.lof + byClass.indel + byClass.splice, syn = byClass.syn;
    const domainDisrupting = variants.filter(v => v.domain!=='—' && (v.consClass==='missense'||v.consClass==='lof'||v.consClass==='indel')).length;
    const afSpectrum = {
      rare:   variants.filter(v => v.af>0 && v.af<0.01).length,
      low:    variants.filter(v => v.af>=0.01 && v.af<0.05).length,
      common: variants.filter(v => v.af>=0.05).length,
    };
    let exonic=0, intronic=0;
    for (const v of variants){
      let inExon;
      if (gm && gm.exons) inExon = gm.exons.some(e => v.pos>=e[0] && v.pos<=e[1]);
      else inExon = (v.consClass!=='other') || /utr/i.test(v.consequence);
      if (/intron/i.test(v.consequence) || (gm && gm.exons && !inExon)) intronic++;
      else if (inExon) exonic++;
    }
    const PR = ['TOP','HIGH','MODERATE','LOW'];
    const damaging = variants
      .filter(v => v.consClass==='lof' || v.severe || (v.consClass==='missense' && v.combined!=null && v.combined<=-4))
      .sort((a,b) => (PR.indexOf(a.priority)-PR.indexOf(b.priority)) || ((a.combined==null?0:a.combined)-(b.combined==null?0:b.combined)));
    const koGenotypes = damaging.filter(v=>v.consClass==='lof').reduce((n,v)=>n+v.hom+v.het, 0);
    const koLines = new Set(); damaging.filter(v=>v.consClass==='lof').forEach(v=>{v.carriersHom.forEach(id=>koLines.add(id));v.carriersHet.forEach(id=>koLines.add(id));});

    return {
      gene, chr:g.chr, start:g.start, end:g.end, strand: gm?gm.strand:null, dataset, datasetName:dsName,
      nAccessions: accs.length, nVariants: variants.length,
      protLen: gd?gd.len:null, protein: gd?gd.protein:null, domains: gd?(gd.domains||[]):[],
      burden: { byClass, nonsyn, syn, nonsynSyn: syn ? +(nonsyn/syn).toFixed(2) : (nonsyn?null:0),
                exonic, intronic, exonIntron: intronic ? +(exonic/intronic).toFixed(2) : (exonic?null:0),
                domainDisrupting, meanByModel:meanScoresOf(variants),
                meanPlantcad:_avg(variants.map(v=>v.pc1).filter(x=>x!=null)),
                meanEsm:_avg(variants.map(v=>v.esm1).filter(x=>x!=null)), afSpectrum },
      damaging, koGenotypes, koLines: koLines.size, variants,
    };
  }

  const DEFAULT_DS = DATASETS[0].id;
  return {
    datasets:    () => DATASETS,
    projectsFor, accessionsFor, defaultSelectionFor, familyOf, dbOf, hasSecondaryScores, scoreModels, namFoundersFor,
    canonicalGeneId, geneCanon, parseSub, classifyConsequence,
    projects:    () => projectsFor(DEFAULT_DS),
    accessions:  () => accessionsFor(DEFAULT_DS),
    defaultSelection: () => defaultSelectionFor(DEFAULT_DS),
    geneModels:  () => ({}),
    exampleGenes,
    lookupGene,
    chromLengths: chromLengths,
    centromeres: () => CENTRO,
    accessionById,
    queryVariants, queryVariantsByGene,
    structureFor, pdbFor, ensureStructure, queryFoldVariants, rankImpact,
    ensureDomains, domainAt,
    ensureGeneDomains, geneDomains, ensureGeneModels, geneModelOf, geneFunction,
  };
})();

function pfamHref(pf){ return 'https://www.ebi.ac.uk/interpro/entry/pfam/' + pf + '/'; }
function domTag(dom){
  if (dom == null || dom === '—' || dom === 'N/A' || dom === '')
    return '<span style="color:var(--faint)">—</span>';
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
