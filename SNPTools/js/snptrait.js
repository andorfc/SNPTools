/* =====================================================================
 *  snptrait.js — SNPTrait = metadata selector ("Line Selector" on maize).
 *  Registers 'snptrait'. Lets the user filter and select samples by
 *  metadata, then hand the selected set to SNPVersity to build a VCF.
 *
 *  Maize (main): the schema for a family comes from
 *  window.SNPTRAIT_SCHEMA[family], written by build_maize_trait_catalog.py
 *  into js/zmgrin.catalog.js (family zmgrin2026: panel, subpopulation,
 *  GRIN origin, kernel type and other coded GRIN descriptors; numeric GRIN
 *  traits are filtered by range on their per-accession means). Families
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
  if (fam && gen[fam]) return withLineQcFacet(gen[fam], dataset);
  if (fam && TRAIT_SCHEMA[fam]) return withLineQcFacet(TRAIT_SCHEMA[fam], dataset);
  return TRAIT_NEUTRAL_SCHEMA;
}
/* Line QC (js/zmgrin.lineqc.js): when the lines carry a sample heterozygosity class, the schema
   gains a "Sample heterozygosity" facet. Added here, on a copy, not in the generated schema. */
function withLineQcFacet(sc, dataset){
  const list = (typeof Data!=='undefined' && Data.accessionsFor) ? (Data.accessionsFor(dataset) || []) : [];
  if (!list.some(a => a.sampleQC) || sc.facets.some(f => f[0] === 'sampleQC')) return sc;
  return Object.assign({}, sc, {facets: sc.facets.concat([['sampleQC', 'Sample heterozygosity']])});
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
  colFilters:{},      // column key -> Set of kept values, for table columns that are not facets
  openFacets:new Set(), // phenotype checkbox sections the user opened (they start collapsed)
  colKeys:null,       // the table's columns after ID, as column keys (null = the schema's default set)
  specs:new Map(),    // column key -> column description (traitColSpec), cached
  sort:null,         // {key, dir: 1 | -1}: the table's sort column (rows sort within each group)
  onlySel:false,      // "Selected only": show just the ticked lines
  mode:'page',        // 'page' (the SNPTrait tool) or 'drawer' (opened over SNPVersity, TraitDrawer)
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
    .then(j => { TRAIT.traitData[fam] = j; TRAIT.specs = new Map(); if (TRAIT.dataset === dataset) traitEnsureFacets();
      if (TRAIT.dataset === dataset && document.getElementById('traitFacets')){ traitRenderFacets(); traitRenderTable(); } })
    .catch(() => { TRAIT.traitData[fam] = null; });
}
function traitNumericTraits(){
  const tf = traitSideFile(); if (!tf || tf === 'loading') return [];
  const n = {}, dict = tf.dictionary || {};
  Object.values(tf.samples || {}).forEach(s => Object.entries(s.traits || {}).forEach(([k, v]) => {
    if (Array.isArray(v) && v.length === 5 && typeof v[1] === 'number') n[k] = (n[k] || 0) + 1; }));
  return Object.entries(n).filter(([k]) => (dict[k] || {}).scale === 'numeric')
    .map(([k, c]) => { const d = dict[k] || {};
      return {code: k, n: c, name: d.name || k, unit: d.unit || '', published: d['class'] === 'LITERATURE_BLUP', url: d.url || ''}; })
    .sort((a, b) => a.name.localeCompare(b.name));
}
function traitValue(id, code){
  const tf = traitSideFile(); if (!tf || tf === 'loading') return null;
  const s = (tf.samples || {})[id]; const v = s && s.traits && s.traits[code];
  return (Array.isArray(v) && v.length === 5 && typeof v[1] === 'number') ? v[1] : null;
}
/* code == null: read the trait, min and max from one of the sidebar's range blocks ('G' or 'P') */
function traitAddRange(code, min, max, grp){
  if (code == null){
    const g = grp || 'G', sel = document.getElementById(`trRange${g}Trait`); if (!sel) return;
    code = sel.value;
    min = document.getElementById(`trRange${g}Min`).value; max = document.getElementById(`trRange${g}Max`).value;
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
/* Two range blocks: GRIN evaluation traits (per-accession means) and published traits (BLUPs
   from a study, not GRIN), each with its own drop-down; ids trRange{G|P}Trait / Min / Max / Hint. */
const TRAIT_RANGE_GROUPS = [
  {g:'G', title:'Continuous GRIN phenotypes', published:false},
  {g:'P', title:'Continuous published phenotypes', published:true},
];
function traitRangeBlockHTML(g){
  const tf = traitSideFile();
  if (!tf) return '';
  if (tf === 'loading') return g === 'G' ? `<div class="facet-block"><h4>Continuous phenotypes</h4><em>Loading trait records…</em></div>` : '';
  const all = traitNumericTraits(); if (!all.length) return '';
  return TRAIT_RANGE_GROUPS.filter(G => !g || G.g === g).map(G => {
    const list = all.filter(t => t.published === G.published); if (!list.length) return '';
    const opts = list.map(t => `<option value="${escAttrT(t.code)}">${escT(t.name)}${t.unit ? ' (' + escT(t.unit) + ')' : ''} · ${t.n}</option>`).join('');
    const active = TRAIT.ranges.map(r => {
      const t = list.find(x => x.code === r.code); if (!t) return '';
      const lab = `${isFinite(r.min) ? traitFmt(r.min) : '−∞'} – ${isFinite(r.max) ? traitFmt(r.max) : '∞'}${t.unit ? ' ' + t.unit : ''}`;
      return `<div class="tr-chip"><span><b>${escT(t.name)}</b> ${escT(lab)}</span>
        <button class="link-btn" onclick="traitRemoveRange('${escAttrT(r.code)}')" title="Remove this range">×</button></div>`;
    }).join('');
    const urls = [...new Set(list.map(t => t.url).filter(Boolean))];
    const note = G.published
      ? `One value per line: best linear unbiased predictions (BLUPs) from a published study${urls.length ? ' (' + urls.map(u => `<a href="${escAttrT(u)}" target="_blank" rel="noopener">${escT(u.replace(/^https?:\/\/(dx\.)?doi\.org\//, 'doi:'))}</a>`).join(', ') + ')' : ''}, not GRIN evaluation records.`
      : 'Per-accession mean of GRIN evaluation records.';
    return `<div class="facet-block tr-range"><h4>${escT(G.title)}</h4>
      <select id="trRange${G.g}Trait" class="tr-sel" onchange="traitRangeHint('${G.g}')" aria-label="${escAttrT(G.title)}">${opts}</select>
      <div class="tr-mm"><input id="trRange${G.g}Min" type="number" step="any" placeholder="min" aria-label="Minimum">
        <span>–</span><input id="trRange${G.g}Max" type="number" step="any" placeholder="max" aria-label="Maximum">
        <button class="qbtn" onclick="traitAddRange(null, null, null, '${G.g}')">Add</button></div>
      <div class="tr-hint" id="trRange${G.g}Hint"></div>
      ${active}
      <div class="tr-note">${note} Lines without a value for a trait are excluded while its range is active.</div>
    </div>`;
  }).join('');
}
function traitRangeHint(grp){
  (grp ? [grp] : TRAIT_RANGE_GROUPS.map(G => G.g)).forEach(g => {
    const el = document.getElementById(`trRange${g}Hint`), sel = document.getElementById(`trRange${g}Trait`);
    if (!el || !sel) return;
    const st = traitRangeStats(sel.value);
    el.textContent = st ? `${st.n} lines with data · min ${traitFmt(st.min)} · median ${traitFmt(st.med)} · max ${traitFmt(st.max)}` : 'No lines with data';
  });
}

function traitLoad(dataset){
  TRAIT.dataset = dataset;
  const sc = TRAIT.schema = traitSchema(dataset);
  TRAIT.org = traitOrganism(dataset);
  const list = Data.accessionsFor(dataset) || [];
  // fields the current schema needs on each row
  const fields = new Set([sc.groupBy, 'strain'].concat(sc.facets.map(f=>f[0]), sc.columns.map(c=>c[0]), sc.search));
  /* every field the line carries (so any of them can be a table column); the schema's own
     fields get 'Unknown' when empty, as the facets expect */
  TRAIT.rows = list.map(a => {
    const row = Object.assign({}, a, {id:a.id, strain:a.strain || a.label || a.id});
    fields.forEach(k => { const v=a[k]; row[k] = (v!=null && v!=='') ? v : (row[k]!=null?row[k]:'Unknown'); });
    return row;
  });
  TRAIT.facets = {}; sc.facets.forEach(([k]) => TRAIT.facets[k] = new Set());
  TRAIT.q = ''; TRAIT.openGroups = new Set(); TRAIT.ranges = []; TRAIT.colFilters = {}; TRAIT.sort = null;
  TRAIT.colKeys = null; TRAIT.specs = new Map(); TRAIT._materialized = null;
  traitEnsureFacets();
  traitEnsureSideFile(dataset);
  TRAIT.hasMeta = TRAIT.rows.some(r => sc.facets.some(([k]) => r[k] && r[k]!=='Unknown' && r[k]!=='unknown'));
  /* Starts empty: it was seeded with SNPVersity's selection, so "replace" in the Send dialog
     sent those lines back unless they were unticked here. Kept between visits (traitLoad runs
     only when the dataset changes). */
  TRAIT.selected = new Set();
}

/* A row passes when, in every section with ticked values, its value is one of them: OR within a
   section, AND between sections. `skip` leaves one section out, for that section's own counts. */
function traitMatch(r){ return traitMatchExcept(r, null); }
function traitMatchExcept(r, skip){
  const sc = TRAIT.schema, f = TRAIT.facets;
  if (TRAIT.onlySel && !TRAIT.selected.has(r.id)) return false;
  for (const {key:k} of traitFacetList()){ if (k !== skip && f[k] && f[k].size && !f[k].has(r[k])) return false; }
  for (const k in TRAIT.colFilters){     // a table column's filter, for a column that is not a facet
    const cf = TRAIT.colFilters[k];
    if (k !== skip && cf && cf.size && !cf.has(String(traitCellRaw(r, traitColSpec(k))))) return false;
  }
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
  traitColMenuClose();
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
         Filter with the facets or the search box, batch-select, then send the selection to SNPVersity to build a VCF.
         Click a column heading to sort by it, or its &#9662; to sort and filter by its values, as in Excel;
         <b>Add/Remove Columns</b> changes the table (line details and GRIN traits), and <b>Download table</b> saves what you built.</p>
    </div></div>

    <div class="ds-picker" id="traitDsPicker"></div>

    <div class="runbar trait-runbar-top" id="traitRunbar"></div>

    ${TRAIT.hasMeta ? '' : `<div class="trait-note">No metadata is available for this dataset yet — showing ${escT(TRAIT.org.one)} IDs only.</div>`}

    ${traitShellHTML()}
  `;
  renderTraitDsPicker();
  traitRenderFacets();
  traitRenderTable();
  traitRenderRunbar();
}

/* Facets, search, quick selections and the grouped table: the same markup on the SNPTrait
   page and in the drawer SNPVersity opens over itself (TraitDrawer below). */
function traitShellHTML(){
  const sc = TRAIT.schema;
  return `<div class="trait-shell">
      <aside class="trait-facets" id="traitFacets"></aside>
      <div class="trait-main">
        <div class="trait-toolbar">
          <div class="trait-search">${ICONS.search}
            <input id="traitSearch" type="search" placeholder="Search ${escAttrT(sc.search.slice(1,5).join(', '))}…"
                   value="${escAttrT(TRAIT.q)}" aria-label="Search ${escAttrT(TRAIT.org.many)}"
                   oninput="TRAIT.q=this.value.trim().toLowerCase();traitRenderTable();traitRenderFacets();traitStatus()">
          </div>
          <div class="trait-qbtns">
            <button class="qbtn" onclick="traitColsPick(event)" aria-haspopup="dialog"
              data-tt="Choose which columns the table shows: line details, GRIN phenotypes and published traits. Download table can save the shown columns.">${ICONS.table} Add/Remove Columns</button>
            ${traitVisBtnHTML()}
            <button class="qbtn ${TRAIT.onlySel?'on':''}" id="traitOnlySel" aria-pressed="${TRAIT.onlySel}" onclick="traitToggleOnlySel()"
              data-tt="On/off: show only the ticked ${escAttrT(TRAIT.org.many)}, hiding the rest. Filters and search still apply on top. Nothing is ticked or unticked.">Selected only</button>
            <button class="qbtn" id="traitRandBtn" onclick="traitRandPick(event)" aria-haspopup="dialog"
              data-tt="Tick a random share (5, 10, 25 or 50%) of the ${escAttrT(TRAIT.org.many)} shown now. Replaces the ticks on the shown ${escAttrT(TRAIT.org.many)}; ticks on hidden ones stay.">Select a random % &#9662;</button>
          </div>
        </div>
        <div class="trait-grid" id="traitGrid"></div>
      </div>
    </div>`;
}
function traitToggleOnlySel(){
  TRAIT.onlySel = !TRAIT.onlySel;
  const b = document.getElementById('traitOnlySel');
  if (b){ b.classList.toggle('on', TRAIT.onlySel); b.setAttribute('aria-pressed', String(TRAIT.onlySel)); }
  traitRenderTable(); traitRenderFacets(); traitStatus();
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

/* ---- facets ----
   traitFacetList(): every checkbox section, as {key, label, section, n}:
     'line' -- the schema's facets that are not phenotypes (panel, origin, subpopulation, ...)
     'grin' -- categorical GRIN phenotypes: the schema facets that hold a decoded GRIN descriptor
               (schema.grinFields: kernel type, cob colour, ...) and every other coded or text
               GRIN descriptor in the trait side-file (key 'tr:<code>'), A-Z
     'pub'  -- categorical traits from a published study (class LITERATURE_BLUP), A-Z
   Continuous traits are not here: they are filtered by range (traitRangeBlockHTML).
   The side-file's categorical values are copied onto the rows as r['tr:<code>'] (decoded label,
   'Unknown' when the line has no record), so they filter, count, sort and show like any field. */
function traitCategoricalTraits(){
  const tf = traitSideFile(); if (!tf || tf === 'loading') return [];
  const dict = tf.dictionary || {}, same = new Set(Object.values(TRAIT.schema.grinFields || {})), n = {};
  Object.values(tf.samples || {}).forEach(s => Object.keys(s.traits || {}).forEach(k => { n[k] = (n[k] || 0) + 1; }));
  return Object.keys(dict).filter(k => n[k] && dict[k].scale && dict[k].scale !== 'numeric' && !same.has(k))
    .map(k => ({code:k, key:'tr:' + k, label:dict[k].name || k, n:n[k], published:dict[k]['class'] === 'LITERATURE_BLUP'}));
}
function traitFacetList(){
  const sc = TRAIT.schema, gf = sc.grinFields || {}, tf = traitSideFile();
  const memo = TRAIT._facetList;
  if (memo && memo.sc === sc && memo.tf === tf && memo.rows === TRAIT.rows) return memo.list;
  const nOf = k => TRAIT.rows.reduce((n, r) => n + (r[k] != null && r[k] !== '' && r[k] !== 'Unknown' ? 1 : 0), 0);
  const byLabel = (a, b) => a.label.localeCompare(b.label, undefined, {numeric:true, sensitivity:'base'});
  const cat = traitCategoricalTraits();
  const line = sc.facets.filter(([k]) => !(k in gf)).map(([k, l]) => ({key:k, label:l, section:'line', n:nOf(k)}));
  const grin = sc.facets.filter(([k]) => k in gf).map(([k, l]) => ({key:k, label:l, section:'grin', n:nOf(k)}))
    .concat(cat.filter(t => !t.published).map(t => ({key:t.key, label:t.label, section:'grin', n:t.n}))).sort(byLabel);
  const pub = cat.filter(t => t.published).map(t => ({key:t.key, label:t.label, section:'pub', n:t.n})).sort(byLabel);
  const list = line.concat(grin, pub);
  TRAIT._facetList = {sc, tf, rows:TRAIT.rows, list};
  return list;
}
/* copy the categorical traits onto the rows (once per side-file) and give every section its Set */
function traitEnsureFacets(){
  const tf = traitSideFile();
  if (tf && tf !== 'loading' && TRAIT._materialized !== tf){
    const cat = traitCategoricalTraits(), dict = tf.dictionary || {};
    TRAIT.rows.forEach(r => {
      const rec = ((tf.samples || {})[r.id] || {}).traits || {};
      cat.forEach(t => {
        const v = rec[t.code], codes = (dict[t.code] || {}).codes;
        r[t.key] = (v && v[1] != null && v[1] !== '') ? (codes && codes[v[1]] != null ? codes[v[1]] : String(v[1])) : 'Unknown';
      });
    });
    TRAIT._materialized = tf;
  }
  traitFacetList().forEach(f => { if (!(TRAIT.facets[f.key] instanceof Set)) TRAIT.facets[f.key] = new Set(); });
}
function facetCounts(rows){
  const keys = traitFacetList().map(f => f.key), z = {};
  keys.forEach(k => z[k] = {});
  rows.forEach(r => keys.forEach(k => { const v = r[k]||'Unknown'; z[k][v] = (z[k][v]||0)+1; }));
  return z;
}
/* Each section counted under every filter but its own ("disjunctive" counts). Counting a section
   under its own filter made ticking one value drop the section's other values to 0, and a 0 is
   not listed, so a section read as pick-one: to get SS or NSS you had to untick, retick and
   reselect. Now a section's other values stay, each with the lines it would add. A ticked value
   stays listed even at 0 (the other sections may exclude it), so it can always be unticked. */
function traitFacetCounts(keys){
  const z = {};
  (keys || traitFacetList().map(f => f.key)).forEach(k => {
    const c = {};
    TRAIT.rows.forEach(r => { if (traitMatchExcept(r, k)){ const v = r[k]||'Unknown'; c[v] = (c[v]||0)+1; } });
    (TRAIT.facets[k] || new Set()).forEach(v => { if (!(v in c)) c[v] = 0; });
    z[k] = c;
  });
  return z;
}
function traitClearFacet(k){
  if (TRAIT.facets[k]) TRAIT.facets[k].clear();
  traitRenderTable(); traitRenderFacets(); traitStatus();
}
/* Order of a section's values (the sidebar and the column menu): a natural order where the data
   has one, else A-Z (numbers in numeric order); 'Unknown' and blanks always last.
     - coded GRIN traits: GRIN's code order (scores 1-9, kernel type A-J, ...)
     - the grouping field: schema.groupOrder (NAM, MAP, WiDiv, Other GRIN)
     - sample heterozygosity: Inbred, Elevated, Heterozygous; yes/no: yes first
     - values that start with a number (scores "4 · Abundant spread…", "7") go by that number, so a
       score GRIN's code list leaves out (Gray Leaf Spot 3, 4, 6, 7) still sorts in place
     - values that differ only by a signed offset (DTS relative: B73-09 … B73+00 … B73+20) by it */
const TRAIT_VALUE_ORDER = {sampleQC:['Inbred', 'Elevated heterozygosity', 'Heterozygous sample'], _yesno:['yes', 'no']};
function traitValueRank(k){
  const sc = TRAIT.schema;
  if (k === sc.groupBy && (sc.groupOrder || []).length) return sc.groupOrder;
  if (TRAIT_VALUE_ORDER[k]) return TRAIT_VALUE_ORDER[k];
  const code = k.startsWith('tr:') ? k.slice(3) : (sc.grinFields || {})[k];
  const codes = code && (traitDict()[code] || {}).codes;
  if (codes) return [...new Set(Object.values(codes))];
  return null;
}
function traitValueCmp(k){
  const rank = traitValueRank(k), at = rank ? new Map(rank.map((v, i) => [v, i])) : null;
  const yn = !rank && TRAIT.rows.length && TRAIT.rows.every(r => ['yes', 'no', 'Unknown', '', undefined].includes(r[k]));
  const last = v => v === '' || v === 'Unknown' || v === 'unknown';
  const lead = v => { const m = /^(-?\d+(?:\.\d+)?)(?:\s|$)/.exec(v); return m ? +m[1] : null; };
  const off = v => { const m = /^(.*?)([+-]\d+(?:\.\d+)?)$/.exec(v); return m ? [m[1], +m[2]] : null; };
  return (a, b) => {
    if (last(a) || last(b)) return last(a) - last(b);
    if (yn) return TRAIT_VALUE_ORDER._yesno.indexOf(a) - TRAIT_VALUE_ORDER._yesno.indexOf(b);
    const na = lead(a), nb = lead(b);
    if (na != null && nb != null && na !== nb) return na - nb;
    const oa = off(a), ob = off(b);
    if (oa && ob && oa[0] === ob[0] && oa[1] !== ob[1]) return oa[1] - ob[1];
    if (at){ const ia = at.has(a) ? at.get(a) : 1e9, ib = at.has(b) ? at.get(b) : 1e9; if (ia !== ib) return ia - ib; }
    return facetLabel(k, a).localeCompare(facetLabel(k, b), undefined, {numeric:true, sensitivity:'base'});
  };
}
function facetLabel(k, v){
  const sc = TRAIT.schema;
  if (k===sc.groupBy && sc.groupNames[v]) return sc.groupNames[v];
  return v;
}
function traitRenderFacets(){
  const box=document.getElementById('traitFacets'); if(!box) return;
  traitEnsureFacets();
  const list = traitFacetList();
  /* every section starts collapsed (its name and how many lines have a value); one opens on a
     click, or by having a ticked value. Only open sections are counted. */
  const isOpen = f => TRAIT.openFacets.has(f.key) || !!(TRAIT.facets[f.key] && TRAIT.facets[f.key].size);
  const C = traitFacetCounts(list.filter(isOpen).map(f => f.key));
  const block = f => {
    const k = f.key, sel = TRAIT.facets[k], label = f.label;
    const ticked = sel.size ? `<span class="facet-any">any of ${sel.size}
        <button class="link-btn" onclick="event.stopPropagation();traitClearFacet('${escAttrT(k)}')" aria-label="Clear ${escAttrT(label)}">clear</button></span>` : '';
    if (f.section === 'line' && !f.n && !sel.size) return '';       // no line has a value
    const open = isOpen(f);
    let items = '';
    if (open){
      const cmp = traitValueCmp(k), entries = Object.entries(C[k] || {}).sort((a, b) => cmp(a[0], b[0]));
      items = entries.map(([v,n])=>`
      <label class="facet-row${n?'':' zero'}">
        <input type="checkbox" ${sel.has(v)?'checked':''}
          onchange="traitToggleFacet('${escAttrT(k)}', this.value, this.checked)" value="${escAttrT(v)}">
        <span class="fv" title="${escAttrT(facetLabel(k,v))}">${escT(facetLabel(k,v))}</span>
        <span class="fn">${sel.size && !sel.has(v) ? '+' : ''}${n}</span></label>`).join('') || '<em>No values</em>';
      /* one box ticks or unticks every value in the section; partly ticked when only some are */
      if (entries.length > 1){
        const all = entries.every(([v]) => sel.has(v)), some = !all && sel.size > 0;
        items = `<label class="facet-all"><input type="checkbox"${all?' checked':''}${some?' data-indet="1"':''}
          onchange="traitFacetAll('${escAttrT(k)}', this.checked)"><span>Select all</span></label>` + items;
      }
    }
    return `<div class="facet-block fb-col${open?' op':''}" data-facet="${escAttrT(k)}">
      <h4><button type="button" class="fb-tog" onclick="traitToggleFacetOpen('${escAttrT(k)}')" aria-expanded="${open?'true':'false'}"
        title="${f.n != null ? escAttrT(f.n + ' lines with a value') : ''}"><span class="caret">${ICONS.caret}</span><span class="fb-l">${escT(label)}</span><span class="fb-n">${f.n != null ? f.n : ''}</span></button>${ticked}</h4>${items}</div>`;
  };
  const sec = id => list.filter(f => f.section === id).map(block).filter(Boolean).join('');
  const line = sec('line'), grin = sec('grin'), pub = sec('pub');
  const rG = traitRangeBlockHTML('G'), rP = traitRangeBlockHTML('P');
  box.innerHTML = `
    <div class="facet-head">
      <span>Filters</span>
      <button class="link-btn" onclick="traitClearFacets()">Clear all</button>
    </div>
    <div class="facet-logic">Click a section to open it; the number beside it is how many lines have a value.
      Tick several values in a section to keep lines with <b>any</b> of them;
      sections combine, so a line must match <b>every</b> section you use. A <b>+n</b> count is what an
      unticked value would add.</div>
    ${traitColFilterBlockHTML()}
    ${line ? `<div class="facet-sec">Line details</div>${line}` : ''}
    ${rG || grin ? `<div class="facet-sec">GRIN phenotypes</div>${rG}${grin ? `<div class="facet-sub">Categorical GRIN phenotypes · A–Z</div>${grin}` : ''}` : ''}
    ${rP || pub ? `<div class="facet-sec">Published phenotypes <span class="facet-sec-n">not from GRIN</span></div>${rP}${pub ? `<div class="facet-sub">Categorical published phenotypes · A–Z</div>${pub}` : ''}` : ''}
    ${line || grin || pub || rG || rP ? '' : '<div class="facet-empty">No metadata to filter.</div>'}`;
  box.querySelectorAll('.facet-all input[data-indet]').forEach(i => { i.indeterminate = true; });
  traitRangeHint();
}
/* Select all: tick every value the section lists (unticked values with no line under the other
   filters too, so the box matches what is shown). Unticking clears the section. */
function traitFacetAll(k, on){
  if (!TRAIT.facets[k]) return;
  if (on) Object.keys(traitFacetCounts([k])[k] || {}).forEach(v => TRAIT.facets[k].add(v));
  else TRAIT.facets[k].clear();
  traitRenderTable(); traitRenderFacets(); traitStatus();
}
function traitToggleFacetOpen(k){
  TRAIT.openFacets.has(k) ? TRAIT.openFacets.delete(k) : TRAIT.openFacets.add(k);
  traitRenderFacets();
}
function traitToggleFacet(k,v,on){
  if(!TRAIT.facets[k]) return;
  on ? TRAIT.facets[k].add(v) : TRAIT.facets[k].delete(v);
  traitRenderTable(); traitRenderFacets(); traitStatus();
}
function traitClearFacets(){
  Object.values(TRAIT.facets).forEach(s=>s.clear()); TRAIT.ranges = []; TRAIT.colFilters = {};
  traitRenderTable(); traitRenderFacets(); traitStatus();
}

/* ---- table grouped by schema.groupBy ---- */
function traitRenderTable(){
  const grid=document.getElementById('traitGrid'); if(!grid) return;
  traitSyncVisBtn();
  const sc=TRAIT.schema, rows=traitVisible();
  const {groups, keys} = traitGroupRows(rows);
  if(!rows.length){ grid.innerHTML=`<div class="trait-empty">No ${escT(TRAIT.org.many)} match the current filters.</div>`; traitStatus(); return; }
  const cols = TRAIT.cols = traitTableCols();
  const thead = `<tr><th class="tcb stk stk0"></th>${cols.map((c,i)=>traitThHTML(c,i)).join('')}</tr>`;
  grid.innerHTML = keys.map(gv=>{
    const items=groups[gv];
    const open = TRAIT.openGroups.has(gv) || TRAIT.q || TRAIT.onlySel || TRAIT.ranges.length || Object.values(TRAIT.facets).some(s=>s.size) || traitColFilterCount();
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
    const rowsHtml = traitSortRows(items, cols).map(r=>`
      <tr class="${TRAIT.selected.has(r.id)?'on':''}" onclick="traitToggle('${escAttrT(r.id)}')">
        <td class="tcb stk stk0"><span class="cb">${ICONS.check}</span></td>
        ${cols.map((c,i)=>{ const t = traitCellText(r,c), clip = t.length > 40;     /* long text: shortened, whole on hover */
          return `<td class="${c.kind==='num'?'num':''}${i===0?' mono stk stk1':''}${clip?' clip':''}"${clip?` title="${escAttrT(t)}"`:''}>${escT(t)}</td>`; }).join('')}
      </tr>`).join('');
    /* the wrapper scrolls sideways when the columns are wider than the page; the tick box and ID stay put */
    return `<div class="tg">${head}<div class="tg-scroll"><table class="tg-table"><thead>${thead}</thead><tbody>${rowsHtml}</tbody></table></div></div>`;
  }).join('');
  traitStatus();
}

/* rows split by schema.groupBy, groups in schema.groupOrder (others after, A-Z) */
function traitGroupRows(rows){
  const sc = TRAIT.schema, gkey = sc.groupBy, groups = {};
  rows.forEach(r => { const g = r[gkey] || 'Unknown'; (groups[g] = groups[g] || []).push(r); });
  const order = sc.groupOrder || [];
  const keys = Object.keys(groups).sort((a, b) => {
    const ia = order.indexOf(a), ib = order.indexOf(b);
    return (ia < 0 ? 99 : ia) - (ib < 0 ? 99 : ib) || a.localeCompare(b);
  });
  return {groups, keys};
}
/* rows in the table's order: by group, then the table's sort within each group */
function traitOrderedRows(rows, cols){
  const {groups, keys} = traitGroupRows(rows);
  return [].concat(...keys.map(g => traitSortRows(groups[g], cols)));
}

/* ---- table columns ----
   A column is named by a key: 'id', a field of the line ('strain', 'country', 'pedigree', ...), or
   'tr:<GRIN descriptor code>' for a GRIN evaluation trait from the side-file (numeric traits show
   the per-accession mean; coded and text traits the most frequent value). The table shows ID,
   then TRAIT.colKeys (the schema's columns until the user picks others with the Columns button),
   then any trait with an active range that is not already shown.
   A heading click sorts by that column: ascending, then descending, then back to the default
   order. Rows sort within each group; blanks and "Unknown" go last either way. */
const TRAIT_FIELD_LABELS = {strain:'Name', grin:'GRIN accession', genotypeId:'Genotype ID', run:'VCF sample name',
  panel:'Primary panel', namFounder:'NAM founder', pedigree:'Pedigree', species:'Species', nTraits:'GRIN traits recorded',
  hetShare:'Heterozygous share (line QC)', sampleQC:'Sample heterozygosity', projTitle:'Project'};
const TRAIT_FIELD_SKIP = new Set(['id', 'label', 'founder', 'rep', 'reps', 'group', 'proj', 'projColor']);   // ID is always first; the rest are internal
const TRAIT_FIELD_NUM = new Set(['hetShare']);
function traitFieldLabel(k){
  const sc = TRAIT.schema;
  const col = (sc.columns || []).find(x => x[0] === k), fac = (sc.facets || []).find(x => x[0] === k);
  return col ? col[1] : fac ? fac[1] : (TRAIT_FIELD_LABELS[k] || k);
}
function traitDict(){ const tf = traitSideFile(); return (tf && tf !== 'loading' && tf.dictionary) || {}; }
function traitColSpec(key){
  if (TRAIT.specs.has(key)) return TRAIT.specs.get(key);
  let c;
  if (key === 'id') c = {key, label:'ID', kind:'text'};
  else if (key.startsWith('tr:')){
    const code = key.slice(3), d = traitDict()[code] || {};
    const num = d.scale === 'numeric', blup = d['class'] === 'LITERATURE_BLUP';   // a published BLUP is not a mean of records
    c = {key, code, trait:true, blup, scale:d.scale || '', kind:num ? 'num' : 'text', unit:d.unit || '', category:d.category || 'Other',
         name:d.name || code,
         label:(d.name || code) + (!num ? '' : blup ? (d.unit ? ` (${d.unit})` : '') : ` (mean${d.unit ? ', ' + d.unit : ''})`)};
  } else c = {key, label:traitFieldLabel(key), kind:TRAIT_FIELD_NUM.has(key) ? 'num' : 'text'};
  if (key === 'id' || !key.startsWith('tr:') || Object.keys(traitDict()).length) TRAIT.specs.set(key, c);   // trait specs once the dictionary is in
  return c;
}
function traitDefaultColKeys(){ return TRAIT.hasMeta ? TRAIT.schema.columns.map(c => c[0]) : []; }
function traitColKeys(){ return TRAIT.colKeys || traitDefaultColKeys(); }
function traitTableCols(){
  const keys = ['id'].concat(traitColKeys());
  TRAIT.ranges.forEach(rg => { if (!keys.includes('tr:' + rg.code)) keys.push('tr:' + rg.code); });
  return keys.map(traitColSpec);
}
/* a GRIN trait record: numeric [n, mean, sd, min, max]; coded [n, mode code]; text [n, mode text] */
function traitRecord(id, code){
  const tf = traitSideFile(); if (!tf || tf === 'loading') return null;
  const s = (tf.samples || {})[id]; return (s && s.traits && s.traits[code]) || null;
}
function traitCellRaw(r, c){
  if (c.trait){
    if (c.kind !== 'num' && r[c.key] !== undefined) return r[c.key];    // copied onto the row (traitEnsureFacets)
    const v = traitRecord(r.id, c.code);
    if (c.kind === 'num') return v && typeof v[1] === 'number' ? v[1] : null;
    if (!v || v[1] == null) return '';
    const codes = (traitDict()[c.code] || {}).codes;
    return codes && codes[v[1]] != null ? codes[v[1]] : String(v[1]);
  }
  const v = r[c.key];
  if (c.kind === 'num') return typeof v === 'number' ? v : (v == null || v === '' ? null : +v);
  return v == null ? '' : v;
}
function traitCellText(r, c){ return c.kind === 'num' ? traitFmt(traitCellRaw(r, c)) : String(traitCellRaw(r, c)); }
function traitSortRows(items, cols){
  const s = TRAIT.sort; if (!s) return items;
  const c = cols.find(x => x.key === s.key); if (!c) return items;   // its trait range was removed
  const blank = v => v == null || v === '' || v === 'Unknown' || v === 'unknown' || (c.kind === 'num' && !isFinite(v));
  return items.slice().sort((a, b) => {
    const va = traitCellRaw(a, c), vb = traitCellRaw(b, c), ba = blank(va), bb = blank(vb);
    if (ba || bb) return ba === bb ? String(a.id).localeCompare(String(b.id)) : (ba ? 1 : -1);
    const d = c.kind === 'num' ? va - vb : String(va).localeCompare(String(vb), undefined, {numeric:true, sensitivity:'base'});
    return d * s.dir || String(a.id).localeCompare(String(b.id));
  });
}
function traitThHTML(c, i){
  const dir = TRAIT.sort && TRAIT.sort.key === c.key ? TRAIT.sort.dir : 0;
  const ind = dir === 1 ? '▲' : dir === -1 ? '▼' : '⇅';
  const next = dir === 1 ? 'descending' : dir === -1 ? 'default order' : 'ascending';
  const filt = traitColFiltered(c);
  return `<th class="${c.kind==='num'?'num':''}${i===0?' stk stk1':''}" aria-sort="${dir===1?'ascending':dir===-1?'descending':'none'}">
    <span class="th-wrap"><button type="button" class="th-sort${dir?' on':''}" onclick="traitSortBy('${escAttrT(c.key)}')"
      title="Sort by ${escAttrT(c.label)} (${next})">${escT(c.label)}<span class="th-ind" aria-hidden="true">${ind}</span></button>
    <button type="button" class="th-filt${filt?' on':''}" onclick="traitColMenu(event, ${i})" aria-haspopup="dialog"
      aria-label="Sort and filter ${escAttrT(c.label)}${filt?' (filtered)':''}" title="Sort and filter${filt?' (filter on)':''}">${filt?TRAIT_FUNNEL:'▾'}</button></span></th>`;
}
function traitSortBy(key){
  const s = TRAIT.sort;
  TRAIT.sort = (!s || s.key !== key) ? {key, dir:1} : (s.dir === 1 ? {key, dir:-1} : null);
  traitRenderTable();
}
/* ---- column filters (Excel-style) ----
   Each heading has a ▾ button that opens a menu: sort, then the column's values with counts, a
   search box and (Select all); OK applies, Cancel / Escape / a click outside drop the edit.
   - A column that is also a facet (Country, Subpopulation, ...) filters through that facet, so
     the menu and the sidebar show and change the same filter.
   - Other columns (ID, GRIN accession, Name) keep their filter in TRAIT.colFilters.
   - A trait-mean column is a GRIN trait range: its menu edits the range's min and max.
   The values listed are those left by every OTHER filter, as in Excel; a ticked value that the
   other filters exclude stays listed (at 0) so it can be unticked. With a search typed, OK keeps
   only the ticked values that match it (Excel's "select all search results"). */
const TRAIT_FUNNEL = '<svg viewBox="0 0 24 24" fill="none" aria-hidden="true"><path d="M4 5h16l-6.2 7.4V19l-3.6-1.8v-4.8z" fill="currentColor"/></svg>';
/* the one open popover: the column menu ('traitColMenu') or the column picker ('traitColsPick') */
const TCM = {id:null, col:null, btn:null, values:[], draft:null, q:''};
function traitIsFacetKey(k){ return TRAIT.facets[k] instanceof Set; }
function traitColCurrent(k){
  const s = traitIsFacetKey(k) ? TRAIT.facets[k] : TRAIT.colFilters[k];
  return s && s.size ? s : null;
}
function traitColRangeable(c){ return !!(c.trait && c.kind === 'num'); }    // a numeric GRIN trait filters by a min-max range
function traitColFiltered(c){
  if (traitColRangeable(c)) return TRAIT.ranges.some(r => r.code === c.code);
  return c.kind === 'num' ? false : !!traitColCurrent(c.key);
}
function traitColFilterCount(){ return Object.values(TRAIT.colFilters).filter(s => s && s.size).length; }
function traitColSet(k, set){
  if (traitIsFacetKey(k)){ TRAIT.facets[k].clear(); if (set) set.forEach(v => TRAIT.facets[k].add(v)); }
  else if (set && set.size) TRAIT.colFilters[k] = set;
  else delete TRAIT.colFilters[k];
}
function traitValueLabel(k, v){ return v === '' ? '(Blanks)' : facetLabel(k, v); }
/* sidebar: the column filters that are not facets (facet ones already show in their section) */
function traitColFilterBlockHTML(){
  const on = Object.entries(TRAIT.colFilters).filter(([, s]) => s && s.size);
  if (!on.length) return '';
  /* name the shorter side: "PI 1, PI 2" when a few are kept, "all but 4: ..." when a few are dropped */
  const list = (k, vals) => vals.slice(0, 3).map(v => traitValueLabel(k, v)).join(', ') + (vals.length > 3 ? ` +${vals.length - 3} more` : '');
  return `<div class="facet-block"><h4>Column filters</h4>${on.map(([k, s]) => {
    const c = traitColSpec(k);
    const every = [...new Set(TRAIT.rows.map(r => String(traitCellRaw(r, c))))];
    const dropped = every.filter(v => !s.has(v));
    const shown = dropped.length && dropped.length < s.size ? `all but ${dropped.length}: ${list(k, dropped)}` : list(k, [...s]);
    return `<div class="tr-chip"><span><b>${escT(c.label)}</b> ${escT(shown)}</span>
      <button class="link-btn" onclick="traitColClearKey('${escAttrT(k)}')" title="Remove this filter">×</button></div>`;
  }).join('')}</div>`;
}
function traitColClearKey(k){
  traitColSet(k, null);
  traitRenderTable(); traitRenderFacets(); traitStatus();
}

/* ---- popover plumbing shared by the column menu and the column picker ---- */
function traitPopOpen(id, btn, html, label, cls){
  const m = document.createElement('div');
  m.id = id; m.className = 'tcm' + (cls ? ' ' + cls : '');
  m.setAttribute('role', 'dialog'); m.setAttribute('aria-label', label);
  m.innerHTML = html;
  document.body.appendChild(m);
  TCM.id = id; TCM.btn = btn;
  document.addEventListener('mousedown', traitColMenuOutside, true);
  document.addEventListener('touchstart', traitColMenuOutside, true);
  window.addEventListener('keydown', traitColMenuKey, true);
  window.addEventListener('resize', traitColMenuPlace);
  window.addEventListener('scroll', traitColMenuPlace, true);
  return m;
}
/* fixed position under its button, kept there while the page or the drawer scrolls; closes if
   the button is gone */
function traitColMenuPlace(){
  const m = TCM.id && document.getElementById(TCM.id); if (!m) return;
  const b = TCM.btn;
  if (!b || !document.body.contains(b)){ traitColMenuClose(); return; }
  const r = b.getBoundingClientRect(), W = window.innerWidth, H = window.innerHeight;
  const w = m.offsetWidth;
  m.style.left = Math.max(8, Math.min(r.left, W - w - 8)) + 'px';
  let top = r.bottom + 4;
  m.style.maxHeight = '';
  const h = m.offsetHeight;
  if (top + h > H - 8) top = Math.max(8, H - 8 - h);
  m.style.top = top + 'px';
  m.style.maxHeight = (H - top - 8) + 'px';
}
function traitColMenuClose(){
  const m = TCM.id && document.getElementById(TCM.id); if (m) m.remove();
  TCM.id = null; TCM.col = null; TCM.btn = null; TCM.values = []; TCM.draft = null; TCM.q = '';
  document.removeEventListener('mousedown', traitColMenuOutside, true);
  document.removeEventListener('touchstart', traitColMenuOutside, true);
  window.removeEventListener('keydown', traitColMenuKey, true);
  window.removeEventListener('resize', traitColMenuPlace);
  window.removeEventListener('scroll', traitColMenuPlace, true);
}
function traitColMenuOutside(e){
  const m = TCM.id && document.getElementById(TCM.id);
  if (!m || m.contains(e.target) || (TCM.btn && TCM.btn.contains(e.target))) return;
  traitColMenuClose();
}
/* capture phase, so Escape closes the popover and not the SNPVersity drawer behind it */
function traitColMenuKey(e){
  if (e.key === 'Escape'){ e.stopPropagation(); e.preventDefault(); traitColMenuClose(); }
  else if (e.key === 'Enter' && TCM.id === 'traitColMenu' && e.target && e.target.tagName === 'INPUT' && e.target.type !== 'checkbox'){
    e.preventDefault(); traitColApply();
  }
}

/* ---- column menu (Excel-style): sort, hide, and filter by value ----
   - A column that is also a facet (Country, Subpopulation, ...) filters through that facet, so
     the menu and the sidebar show and change the same filter.
   - Other columns (ID, GRIN accession, Name, coded GRIN traits, ...) keep their filter in
     TRAIT.colFilters.
   - A numeric GRIN trait filters by a GRIN trait range: its menu edits the range's min and max.
   The values listed are those left by every OTHER filter, as in Excel; a ticked value that the
   other filters exclude stays listed (at 0) so it can be unticked. With a search typed, OK keeps
   only the ticked values that match it (Excel's "select all search results"). */
function traitColMenu(ev, i){
  if (ev){ ev.stopPropagation(); ev.preventDefault(); }
  const c = (TRAIT.cols || [])[i]; if (!c) return;
  const btn = ev && ev.currentTarget;
  const wasThis = TCM.id === 'traitColMenu' && TCM.col && TCM.col.key === c.key;
  traitColMenuClose();
  if (wasThis) return;                                   // a second click on the same ▾ closes it
  const dir = TRAIT.sort && TRAIT.sort.key === c.key ? TRAIT.sort.dir : 0;
  const az = c.kind === 'num' ? ['Sort smallest to largest', 'Sort largest to smallest'] : ['Sort A to Z', 'Sort Z to A'];
  const ranged = traitColRangeable(c), sortOnly = c.kind === 'num' && !ranged;
  const inKeys = traitColKeys().includes(c.key);
  const rangeOn = ranged && TRAIT.ranges.some(r => r.code === c.code);
  let body = '', vals = [], draft = null;
  if (ranged){
    const rg = TRAIT.ranges.find(r => r.code === c.code) || {}, st = traitRangeStats(c.code);
    body = `<div class="tcm-lbl">Keep lines whose ${c.blup ? 'value' : 'mean'} is between</div>
      <div class="tcm-mm"><input id="tcmMin" type="number" step="any" placeholder="min" value="${isFinite(rg.min) ? rg.min : ''}">
        <span>–</span><input id="tcmMax" type="number" step="any" placeholder="max" value="${isFinite(rg.max) ? rg.max : ''}"></div>
      <div class="tcm-hint">${st ? `${st.n} lines with data · min ${traitFmt(st.min)} · median ${traitFmt(st.med)} · max ${traitFmt(st.max)}` : 'No lines with data'}</div>`;
  } else if (sortOnly){
    body = `<div class="tcm-hint">This column can be sorted but not filtered.</div>`;
  } else {
    const counts = {};
    TRAIT.rows.forEach(r => { if (traitMatchExcept(r, c.key)){ const v = String(traitCellRaw(r, c)); counts[v] = (counts[v] || 0) + 1; } });
    const cur = traitColCurrent(c.key);
    if (cur) cur.forEach(v => { if (!(v in counts)) counts[v] = 0; });
    const cmp = traitValueCmp(c.key);
    vals = Object.entries(counts).sort((a, b) => cmp(a[0], b[0]));
    draft = new Set(cur ? [...cur] : Object.keys(counts));
    body = `<div class="tcm-search">${ICONS.search}<input id="tcmQ" type="search" placeholder="Search" aria-label="Search ${escAttrT(c.label)} values"
        oninput="traitColMenuSearch(this.value)"></div>
      <label class="tcm-row tcm-all"><input type="checkbox" id="tcmAll" onchange="traitColMenuAll(this.checked)"><span id="tcmAllLbl">(Select all)</span></label>
      <div class="tcm-list" id="tcmList"></div>`;
  }
  const clearLbl = ranged ? (inKeys ? 'Clear filter' : 'Remove filter and column') : 'Clear filter';
  traitPopOpen('traitColMenu', btn, `<div class="tcm-head">${escT(c.label)}</div>
    <button type="button" class="tcm-item${dir === 1 ? ' on' : ''}" onclick="traitColMenuSort(1)">▲ ${az[0]}</button>
    <button type="button" class="tcm-item${dir === -1 ? ' on' : ''}" onclick="traitColMenuSort(-1)">▼ ${az[1]}</button>
    ${dir ? '<button type="button" class="tcm-item" onclick="traitColMenuSort(0)">⇅ Clear sort</button>' : ''}
    ${inKeys && !rangeOn ? '<button type="button" class="tcm-item" onclick="traitColHide()">✕ Hide column</button>' : ''}
    <div class="tcm-sep"></div>
    ${body}
    <div class="tcm-foot">
      ${sortOnly ? '' : `<button type="button" class="link-btn" onclick="traitColClear()"${traitColFiltered(c) ? '' : ' disabled'}>${clearLbl}</button>`}
      <span class="tcm-sp"></span>
      <button type="button" class="qbtn" onclick="traitColMenuClose()">${sortOnly ? 'Close' : 'Cancel'}</button>
      ${sortOnly ? '' : '<button type="button" class="qbtn solid" id="tcmOk" onclick="traitColApply()">OK</button>'}
    </div>`, 'Sort and filter ' + c.label);
  TCM.col = c; TCM.q = ''; TCM.values = vals; TCM.draft = draft;
  if (draft) traitColMenuList();
  traitColMenuPlace();
  const f = document.getElementById(ranged ? 'tcmMin' : 'tcmQ'); if (f) f.focus();
}
function traitColMenuVisible(){
  const q = TCM.q.toLowerCase(), k = TCM.col ? TCM.col.key : '';
  return q ? TCM.values.filter(([v]) => traitValueLabel(k, v).toLowerCase().includes(q)) : TCM.values;
}
function traitColMenuList(){
  const list = document.getElementById('tcmList'); if (!list || !TCM.col) return;
  const k = TCM.col.key, vis = traitColMenuVisible();
  list.innerHTML = vis.length ? vis.map(([v, n]) => `<label class="tcm-row${n ? '' : ' zero'}">
      <input type="checkbox" value="${escAttrT(v)}"${TCM.draft.has(v) ? ' checked' : ''} onchange="traitColMenuTick(this.value, this.checked)">
      <span class="tv" title="${escAttrT(traitValueLabel(k, v))}">${escT(traitValueLabel(k, v))}</span><span class="tn">${n}</span></label>`).join('')
    : '<div class="tcm-none">No matching values</div>';
  traitColMenuSync();
}
function traitColMenuSync(){
  const vis = traitColMenuVisible(), on = vis.filter(([v]) => TCM.draft.has(v)).length;
  const all = document.getElementById('tcmAll');
  if (all){ all.checked = vis.length > 0 && on === vis.length; all.indeterminate = on > 0 && on < vis.length; all.disabled = !vis.length; }
  const lbl = document.getElementById('tcmAllLbl'); if (lbl) lbl.textContent = TCM.q ? '(Select all search results)' : '(Select all)';
  const ok = document.getElementById('tcmOk'); if (ok) ok.disabled = !on;
}
function traitColMenuSearch(q){ TCM.q = String(q || '').trim(); traitColMenuList(); }
function traitColMenuTick(v, on){ on ? TCM.draft.add(v) : TCM.draft.delete(v); traitColMenuSync(); }
function traitColMenuAll(on){
  traitColMenuVisible().forEach(([v]) => on ? TCM.draft.add(v) : TCM.draft.delete(v));
  traitColMenuList();
}
function traitColMenuSort(dir){
  const c = TCM.col; if (!c) return;
  TRAIT.sort = dir ? {key:c.key, dir} : null;
  traitColMenuClose(); traitRenderTable();
}
function traitColApply(){
  const c = TCM.col; if (!c) return;
  if (traitColRangeable(c)){
    const mn = document.getElementById('tcmMin').value, mx = document.getElementById('tcmMax').value;
    traitColMenuClose(); traitAddRange(c.code, mn, mx); return;
  }
  if (!TCM.draft) return;
  let keep = [...TCM.draft];
  if (TCM.q){ const vis = new Set(traitColMenuVisible().map(x => x[0])); keep = keep.filter(v => vis.has(v)); }
  if (!keep.length) return;
  const kept = new Set(keep);
  const every = TCM.values.every(([v]) => kept.has(v));
  traitColSet(c.key, every ? null : kept);             // every value ticked = no filter
  traitColMenuClose();
  traitRenderTable(); traitRenderFacets(); traitStatus();
}
function traitColClear(){
  const c = TCM.col; if (!c) return;
  traitColMenuClose();
  if (traitColRangeable(c)){
    if (TRAIT.sort && TRAIT.sort.key === c.key && !traitColKeys().includes(c.key)) TRAIT.sort = null;
    traitRemoveRange(c.code); return;
  }
  traitColClearKey(c.key);
}
/* hiding a column keeps its filter (as in Excel); the sidebar still lists it */
function traitColHide(){
  const c = TCM.col; if (!c) return;
  traitColMenuClose();
  traitSetColKeys(traitColKeys().filter(k => k !== c.key));
}

/* ---- column picker (the Columns button): add or remove any column the data has ----
   Line details: every field the lines carry. GRIN evaluation traits: every descriptor in the
   trait side-file, by category. Each entry says how many lines have a value. */
function traitSetColKeys(keys){
  const def = traitDefaultColKeys();
  TRAIT.colKeys = (keys.length === def.length && keys.every((k, i) => k === def[i])) ? null : keys;
  traitRenderTable();
}
function traitColsReset(){ traitSetColKeys(traitDefaultColKeys()); traitColsPickList(); }
function traitAvailableFields(){
  const sc = TRAIT.schema, n = {};
  TRAIT.rows.forEach(r => Object.keys(r).forEach(k => {
    const v = r[k];
    if (TRAIT_FIELD_SKIP.has(k) || k.startsWith('tr:') || v == null || v === '' || v === 'Unknown' || typeof v === 'object') return;
    n[k] = (n[k] || 0) + 1;
  }));
  const order = sc.columns.map(c => c[0]).concat(sc.facets.map(f => f[0]));
  return Object.keys(n).sort((a, b) => {
    const ia = order.indexOf(a), ib = order.indexOf(b);
    return (ia < 0 ? 999 : ia) - (ib < 0 ? 999 : ib) || traitFieldLabel(a).localeCompare(traitFieldLabel(b));
  }).map(k => ({key:k, label:traitFieldLabel(k), n:n[k]}));
}
function traitAvailableTraits(){
  const tf = traitSideFile(); if (!tf || tf === 'loading') return [];
  const n = {};
  Object.values(tf.samples || {}).forEach(s => Object.keys(s.traits || {}).forEach(k => { n[k] = (n[k] || 0) + 1; }));
  return Object.keys(tf.dictionary || {}).filter(k => n[k]).map(k => {
    const c = traitColSpec('tr:' + k), d = tf.dictionary[k];
    return {key:c.key, code:k, label:c.label, category:c.category, n:n[k],
            literature:d['class'] === 'LITERATURE_BLUP', url:d.url || ''};
  }).sort(traitByLabel);
}
function traitByLabel(a, b){ return a.label.localeCompare(b.label, undefined, {numeric:true, sensitivity:'base'}); }
/* The columns on offer, in sections (the Columns picker and the "everything" download):
   - Line details: passport, panel and line-QC fields.
   - GRIN phenotypes, A-Z: every GRIN evaluation descriptor in the side-file, plus the line
     fields that hold a decoded GRIN descriptor (schema.grinFields: kernel type, cob colour, ...),
     which stand in for that descriptor (the same values, and their filter is the sidebar's).
   - Published traits: values from a study rather than GRIN (class LITERATURE_BLUP). */
function traitColumnSections(){
  const gf = TRAIT.schema.grinFields || {}, same = new Set(Object.values(gf));
  const fields = traitAvailableFields();
  const pheno = fields.filter(f => f.key in gf), lit = [];
  traitAvailableTraits().forEach(t => { if (!same.has(t.code)) (t.literature ? lit : pheno).push(t); });
  const urls = [...new Set(lit.map(t => t.url).filter(Boolean))];
  return [
    {title:'Line details', items:fields.filter(f => !(f.key in gf))},
    {title:'GRIN phenotypes', items:pheno.sort(traitByLabel)},
    {title:'Published traits (not from GRIN)', items:lit,
     note:'BLUPs from a published study' + (urls.length ? ': ' + urls.join(', ') : '') + '. Not GRIN evaluation records.'},
  ];
}
function traitColsPick(ev){
  if (ev){ ev.stopPropagation(); ev.preventDefault(); }
  const btn = ev && ev.currentTarget;
  const wasThis = TCM.id === 'traitColsPick';
  traitColMenuClose();
  if (wasThis) return;
  traitPopOpen('traitColsPick', btn, `<div class="tcm-head">Table columns</div>
    <div class="tcm-search">${ICONS.search}<input id="tcpQ" type="search" placeholder="Search columns" aria-label="Search columns"
      oninput="TCM.q=this.value.trim();traitColsPickList()"></div>
    <div class="tcm-list tcp-list" id="tcpList"></div>
    <div class="tcm-foot">
      <button type="button" class="link-btn" onclick="traitColsReset()">Reset to default</button>
      <span class="tcm-sp"></span>
      <button type="button" class="qbtn solid" onclick="traitColMenuClose()">Done</button>
    </div>`, 'Choose table columns', 'tcp');
  traitColsPickList();
  traitColMenuPlace();
  const q = document.getElementById('tcpQ'); if (q) q.focus();
}
function traitColsPickList(){
  const list = document.getElementById('tcpList'); if (!list) return;
  const q = (TCM.q || '').toLowerCase();
  const shown = new Set(traitTableCols().map(c => c.key)), keys = traitColKeys();
  const ranged = new Set(TRAIT.ranges.map(r => 'tr:' + r.code));
  const row = e => {
    const lock = ranged.has(e.key) && !keys.includes(e.key);
    return `<label class="tcm-row"${lock ? ' title="Shown while its range filter is active"' : ''}>
      <input type="checkbox" value="${escAttrT(e.key)}"${shown.has(e.key) ? ' checked' : ''}${lock ? ' disabled' : ''}
        onchange="traitColsToggle(this.value, this.checked)">
      <span class="tv" title="${escAttrT(e.label)}">${escT(e.label)}</span><span class="tn">${e.n}</span></label>`;
  };
  const match = e => !q || e.label.toLowerCase().includes(q) || (e.category || '').toLowerCase().includes(q);
  let html = '<div class="tcp-sec">Always shown: ID</div>';
  traitColumnSections().forEach(sec => {
    const items = sec.items.filter(match);
    if (!items.length) return;
    html += `<div class="tcp-sec">${escT(sec.title)}</div>${sec.note ? `<div class="tcp-note">${escT(sec.note)}</div>` : ''}${items.map(row).join('')}`;
  });
  const tf = traitSideFile();
  if (tf === 'loading') html += '<div class="tcm-none">Loading GRIN phenotypes…</div>';
  list.innerHTML = html.includes('tcm-row') || tf === 'loading' ? html : '<div class="tcm-none">No matching columns</div>';
}
function traitColsToggle(key, on){
  const keys = traitColKeys().filter(k => k !== key);
  if (on) keys.push(key);
  traitSetColKeys(keys);
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
/* One button for both: "Select visible" ticks every line shown; once they all are, it reads
   "Deselect visible" and unticks them. The label follows the table (filters, ticks), so it is
   re-drawn by traitRenderTable. */
function traitVisAllSelected(){
  const vis = traitVisible();
  return vis.length > 0 && vis.every(r => TRAIT.selected.has(r.id));
}
function traitVisBtnState(){
  const all = traitVisAllSelected(), one = TRAIT.org.one, many = TRAIT.org.many;
  return {all, label:all ? 'Deselect visible' : 'Select visible',
    tip:all ? `Untick every ${one} the filters and search leave shown (in all groups, collapsed or not). Ticks on hidden ${many} stay.`
            : `Tick every ${one} the filters and search leave shown (in all groups, collapsed or not). Hidden ${many} are not changed.`};
}
function traitVisBtnHTML(){
  const s = traitVisBtnState();
  return `<button class="qbtn${s.all ? ' on' : ''}" id="traitVisBtn" onclick="traitToggleVisible()" data-tt="${escAttrT(s.tip)}">${s.label}</button>`;
}
/* in place, so the button keeps keyboard focus and its tooltip */
function traitSyncVisBtn(){
  const b = document.getElementById('traitVisBtn'); if (!b) return;
  const s = traitVisBtnState();
  b.textContent = s.label; b.classList.toggle('on', s.all); b.dataset.tt = s.tip;
  if (typeof ttEl !== 'undefined' && ttEl === b) tt.textContent = s.tip;
}
function traitToggleVisible(){ traitSelectVisible(!traitVisAllSelected()); }
/* "Select a random %": a small menu of shares */
function traitRandPick(ev){
  if (ev){ ev.stopPropagation(); ev.preventDefault(); }
  const btn = ev && ev.currentTarget;
  const wasThis = TCM.id === 'traitRandPick';
  traitColMenuClose();
  if (wasThis) return;
  const opt = (f, l) => `<button type="button" class="tcm-item" onclick="traitColMenuClose();traitRandom(${f})">${l}</button>`;
  traitPopOpen('traitRandPick', btn, `<div class="tcm-head">Tick a random share of the shown ${escT(TRAIT.org.many)}</div>
    ${opt(.05, '5%')}${opt(.10, '10%')}${opt(.25, '25%')}${opt(.50, '50%')}`, 'Select a random percentage', 'trp');
  traitColMenuPlace();
  const first = document.querySelector('#traitRandPick .tcm-item'); if (first) first.focus();
}
function traitRandom(frac){
  const vis=traitVisible();
  if(!vis.length) return;                  // nothing matches the filters, so nothing to draw from
  const n=Math.max(1, Math.round(vis.length*frac));
  vis.forEach(r=>TRAIT.selected.delete(r.id));
  randomSample(vis, n).forEach(r=>TRAIT.selected.add(r.id));
  traitRenderTable(); traitRenderRunbar();
}

function traitStatus(){
  const el=document.getElementById('traitStatus'); if(!el) return;
  el.textContent = `Visible ${traitVisible().length} / ${TRAIT.rows.length} · Selected ${TRAIT.selected.size}`;
}

/* ---- run bar: export + hand off to SNPVersity ---- */
function traitRenderRunbar(){
  if(TRAIT.mode==='drawer'){ TraitDrawer.renderFoot(); return; }
  const rb=document.getElementById('traitRunbar'); if(!rb) return;
  const n=TRAIT.selected.size;
  rb.innerHTML=`
    <div class="summ">
      <span class="kick${n?'':' wait'}">${n?'Ready':'Select '+escT(TRAIT.org.many)+' to continue'}</span>
      <b id="traitStatus"></b>
    </div>
    <div class="spacer"></div>
    <button class="btn" onclick="traitOpenDownload()" aria-haspopup="dialog"
      title="Download the table: choose the lines, the columns and the format">${ICONS.download} Download table…</button>
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
/* ---- download ----
   traitExport(kind, opts) writes a table file.
     kind       'csv' | 'tsv' | 'json'
     opts.rows  'selected' (default) | 'visible' (the lines the filters leave) | 'all'
     opts.cols  'standard' (default): ID, group, the schema's columns and facets, line QC and the
                    active trait ranges, headed by field key (the set this export has always had;
                    not offered in the Download dialog, kept for callers such as the harness)
                'table': the columns shown in the table, in order (plus the group column, since
                    the table is grouped by it), headed by column name
                'all': every line detail and every GRIN trait, headed by column name
     opts.stats with 'table' or 'all': n, SD, min and max beside each numeric GRIN trait mean
   Rows come in the table's order (by group, then the table's sort). Values are the ones the table
   shows, with numbers at full precision; coded GRIN traits are written as their decoded label. */
function traitExportRows(which){
  const rows = which === 'all' ? TRAIT.rows : which === 'visible' ? traitVisible() : TRAIT.rows.filter(r => TRAIT.selected.has(r.id));
  return traitOrderedRows(rows, TRAIT.cols || traitTableCols());
}
function traitExportCols(which, stats){
  if (which === 'table' || which === 'all'){
    let specs;
    if (which === 'table'){
      specs = traitTableCols();
      const g = TRAIT.schema.groupBy;
      if (g && !specs.some(c => c.key === g)) specs.splice(1, 0, traitColSpec(g));
    } else {
      specs = [traitColSpec('id')].concat(...traitColumnSections().map(sec => sec.items.map(e => traitColSpec(e.key))));
    }
    const out = [];
    specs.forEach(c => {
      out.push({head:c.label, get:r => traitCellRaw(r, c)});
      if (stats && traitColRangeable(c)){
        const at = (i, lbl) => out.push({head:`${c.name} ${lbl}`, get:r => { const v = traitRecord(r.id, c.code); return v && typeof v[i] === 'number' ? v[i] : (i === 0 && v ? v[0] : null); }});
        at(0, 'n'); at(2, 'SD'); at(3, 'min'); at(4, 'max');
      }
    });
    const seen = {};                                      // unique headings
    out.forEach(o => { const h = o.head; seen[h] = (seen[h] || 0) + 1; if (seen[h] > 1) o.head = `${h} (${seen[h]})`; });
    return out;
  }
  const sc = TRAIT.schema;
  const keys = [...new Set(['id', sc.groupBy].concat(sc.columns.map(c => c[0]), sc.facets.map(f => f[0])))];
  if (TRAIT.rows.some(r => r.hetShare != null)) keys.push('hetShare');   // line QC: share of clean sites heterozygous
  return keys.map(k => ({head:k, get:r => r[k]}))
    .concat(TRAIT.ranges.map(rg => ({head:rg.code, get:r => traitValue(r.id, rg.code)})));
}
function traitExport(kind, opts){
  opts = opts || {};
  const which = opts.rows || 'selected';
  const rows = traitExportRows(which);
  if (!rows.length) return;
  const cols = traitExportCols(opts.cols || 'standard', !!opts.stats);
  /* a missing value is written as an empty cell (null in JSON). 'Unknown' is the placeholder the
     table and filters give a line with no value (traitLoad, traitEnsureFacets); it is not data. */
  const missing = v => v == null || v === '' || v === 'Unknown' || v === 'unknown' || (typeof v === 'number' && !isFinite(v));
  const cell = v => missing(v) ? '' : String(v);
  const many = TRAIT.org.many;
  const base = which === 'all' ? 'all_' + many : which === 'visible' ? 'shown_' + many : TRAIT.org.exportBase;
  let blob, name;
  if (kind === 'json'){
    const data = rows.map(r => { const o = {}; cols.forEach(c => { const v = c.get(r); o[c.head] = missing(v) ? null : v; }); return o; });
    blob = new Blob([JSON.stringify(data, null, 2)], {type:'application/json'}); name = base + '.json';
  } else if (kind === 'tsv'){
    const esc = v => cell(v).replace(/[\t\r\n]+/g, ' ');
    const txt = [cols.map(c => esc(c.head)).join('\t')].concat(rows.map(r => cols.map(c => esc(c.get(r))).join('\t'))).join('\n');
    blob = new Blob(['﻿' + txt], {type:'text/tab-separated-values'}); name = base + '.tsv';
  } else {
    const esc = v => { v = cell(v); return /[",\r\n]/.test(v) ? `"${v.replace(/"/g, '""')}"` : v; };
    const txt = [cols.map(c => esc(c.head)).join(',')].concat(rows.map(r => cols.map(c => esc(c.get(r))).join(','))).join('\n');
    /* a byte-order mark, so Excel reads the file as UTF-8 (accented place names) */
    blob = new Blob([(opts.cols && opts.cols !== 'standard' ? '﻿' : '') + txt], {type:'text/csv'}); name = base + '.csv';
  }
  const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = name;
  document.body.appendChild(a); a.click(); a.remove(); setTimeout(() => URL.revokeObjectURL(a.href), 1500);
}
/* the Download dialog: which lines, which columns, which format */
function traitOpenDownload(){
  traitCloseDownload(); traitColMenuClose();
  const nSel = TRAIT.selected.size, nVis = traitVisible().length, nAll = TRAIT.rows.length;
  const unit = n => escT(traitUnit(n));
  const nTable = traitExportCols('table').length, nEvery = traitExportCols('all').length;
  const radio = (name, val, on, dis, html) => `<label class="trait-send-check${dis ? ' dis' : ''}"><input type="radio" name="${name}" value="${val}"${on ? ' checked' : ''}${dis ? ' disabled' : ''}><span>${html}</span></label>`;
  const host = document.createElement('div'); host.id = 'traitDlDlg';
  host.innerHTML = `<div class="trait-send-backdrop" onclick="traitCloseDownload()"></div>
    <div class="trait-send" role="dialog" aria-modal="true" aria-labelledby="traitDlTitle">
      <div class="trait-send-head"><h3 id="traitDlTitle">Download table</h3>
        <button type="button" class="trait-send-x" onclick="traitCloseDownload()" aria-label="Close">&times;</button></div>
      <div class="trait-send-body">
        <fieldset class="tdl-set"><legend>Rows</legend>
          ${radio('tdlRows', 'selected', nSel > 0, !nSel, `Selected ${unit(nSel)} <span class="trait-send-n">(${nSel})</span>`)}
          ${radio('tdlRows', 'visible', !nSel, !nVis, `Current filter results <span class="trait-send-n">(${nVis})</span>`)}
          ${radio('tdlRows', 'all', false, false, `All ${unit(nAll)} in the dataset <span class="trait-send-n">(${nAll})</span>`)}
        </fieldset>
        <fieldset class="tdl-set"><legend>Columns</legend>
          ${radio('tdlCols', 'table', true, false, `Current table view <span class="trait-send-n">(${nTable} columns)</span>`)}
          ${radio('tdlCols', 'all', false, false, `Everything available <span class="trait-send-n">(${nEvery} columns)</span>`)}
          <label class="trait-send-check tdl-sub"><input type="checkbox" id="tdlStats"><span>Add n, SD, min and max beside each numeric GRIN trait mean</span></label>
        </fieldset>
        <fieldset class="tdl-set"><legend>Format</legend>
          <div class="tdl-fmt">${radio('tdlFmt', 'csv', true, false, 'CSV')}${radio('tdlFmt', 'tsv', false, false, 'TSV')}${radio('tdlFmt', 'json', false, false, 'JSON')}</div>
        </fieldset>
      </div>
      <div class="trait-send-foot">
        <button type="button" class="btn" onclick="traitCloseDownload()">Cancel</button>
        <button type="button" class="btn primary" id="tdlGo" onclick="traitDownloadGo()">${ICONS.download} Download</button>
      </div>
    </div>`;
  document.body.appendChild(host);
  document.addEventListener('keydown', traitDlKey);
  const go = document.getElementById('tdlGo'); if (go) go.focus();
}
function traitDlKey(e){ if (e.key === 'Escape') traitCloseDownload(); }
function traitCloseDownload(){
  const el = document.getElementById('traitDlDlg'); if (el) el.remove();
  document.removeEventListener('keydown', traitDlKey);
}
function traitDownloadGo(){
  const pick = n => { const i = document.querySelector(`input[name="${n}"]:checked`); return i ? i.value : null; };
  const opts = {rows:pick('tdlRows') || 'visible', cols:pick('tdlCols') || 'table', stats:!!(document.getElementById('tdlStats') || {}).checked};
  const fmt = pick('tdlFmt') || 'csv';
  traitCloseDownload();
  traitExport(fmt, opts);
}

/* =====================================================================
 *  TraitDrawer -- SNPTrait opened over SNPVersity to edit its selection in
 *  place: the same facets, trait ranges, search and grouped table, on a DRAFT
 *  of the selection it is given. Apply hands the draft to opts.onApply; Cancel,
 *  Escape and the backdrop drop it. While the drawer is open the SNPTrait
 *  page's own state (its selection and filters) is set aside and put back on
 *  close, and the drawer keeps its own filters between openings (per dataset).
 *
 *    TraitDrawer.open({dataset, selected:[ids], onApply(ids), title})
 * ===================================================================== */
const TraitDrawer = (function(){
  const KEYS = ['selected','q','facets','ranges','openGroups','onlySel','sort','colFilters'];
  let pageState = null, opts = null, base = new Set(), returnFocus = null, pageNodes = null;
  const drawerState = {};                              // dataset -> its filters, kept between openings
  const take = () => { const o = {}; KEYS.forEach(k => { o[k] = TRAIT[k]; }); return o; };
  const put  = o => KEYS.forEach(k => { TRAIT[k] = o[k]; });
  function freshFilters(){
    const f = {}; traitFacetList().forEach(({key}) => { f[key] = new Set(); });
    return {q:'', facets:f, ranges:[], openGroups:new Set(), onlySel:false, sort:null, colFilters:{}};
  }
  function isOpen(){ return !!document.getElementById('traitDrawer'); }

  function open(o){
    if (isOpen()) return;
    traitColMenuClose();
    opts = o || {};
    injectTraitCSS();
    const ds = opts.dataset || S.dataset;
    if (TRAIT.dataset !== ds) traitLoad(ds);
    pageState = take();
    /* Opened from the SNPTrait page itself (the top bar's selection chip): the page and the
       drawer use the same element ids, so the page's nodes wait outside the document. */
    const page = document.getElementById('page');
    if (page && page.querySelector('#traitGrid')){
      pageNodes = document.createDocumentFragment();
      while (page.firstChild) pageNodes.appendChild(page.firstChild);
    }
    const known = new Set(TRAIT.rows.map(r => r.id));
    base = new Set((opts.selected || []).filter(id => known.has(id)));
    put(Object.assign({}, drawerState[ds] || freshFilters(), {selected: new Set(base)}));
    TRAIT.mode = 'drawer';
    returnFocus = document.activeElement;
    const dsName = (Data.datasets().find(d => d.id === ds) || {}).name || ds;
    const host = document.createElement('div'); host.id = 'traitDrawer';
    host.innerHTML = `<div class="tdr-backdrop" onclick="TraitDrawer.close(false)"></div>
      <section class="tdr" role="dialog" aria-modal="true" aria-labelledby="tdrTitle">
        <header class="tdr-head">
          <div><h2 id="tdrTitle">${escT(opts.title || 'Choose ' + TRAIT.org.many)}</h2>
            <div class="tdr-sub">${escT(dsName)} · filter by panel, origin, subpopulation and GRIN traits, then tick
              ${escT(TRAIT.org.many)}. The selection in SNPVersity changes only when you apply.</div></div>
          <button type="button" class="tdr-x" onclick="TraitDrawer.close(false)" aria-label="Close without applying">&times;</button>
        </header>
        <div class="tdr-body">
          ${TRAIT.hasMeta ? '' : `<div class="trait-note">No metadata is available for this dataset yet — showing ${escT(TRAIT.org.one)} IDs only.</div>`}
          ${traitShellHTML()}
        </div>
        <footer class="tdr-foot" id="tdrFoot"></footer>
      </section>`;
    document.body.appendChild(host);
    document.body.classList.add('tdr-lock');
    document.addEventListener('keydown', onKey);
    traitRenderFacets(); traitRenderTable(); renderFoot();
    const q = document.getElementById('traitSearch'); if (q) q.focus();
  }
  function onKey(e){ if (e.key === 'Escape'){ e.preventDefault(); close(false); } }
  function renderFoot(){
    const f = document.getElementById('tdrFoot'); if (!f) return;
    const n = TRAIT.selected.size;
    let add = 0, rem = 0;
    TRAIT.selected.forEach(id => { if (!base.has(id)) add++; });
    base.forEach(id => { if (!TRAIT.selected.has(id)) rem++; });
    const change = (add || rem)
      ? `${add ? `<b class="tdr-add">+${add}</b>` : ''}${add && rem ? ' · ' : ''}${rem ? `<b class="tdr-rem">−${rem}</b>` : ''} vs the ${base.size} selected in SNPVersity`
      : `same as the ${base.size} selected in SNPVersity`;
    f.innerHTML = `<div class="tdr-summ"><span class="tdr-n">${n}</span> ${escT(traitUnit(n))} selected
        <span class="tdr-delta">${change}</span><div class="tdr-status" id="traitStatus"></div></div>
      <div class="tdr-acts">
        <button type="button" class="btn" onclick="TraitDrawer.close(false)">Cancel</button>
        <button type="button" class="btn primary" id="tdrApply" onclick="TraitDrawer.close(true)">Apply to SNPVersity</button>
      </div>`;
    traitStatus();
  }
  function close(apply){
    if (!isOpen()) return;
    const chosen = [...TRAIT.selected];
    const st = take(); delete st.selected; drawerState[TRAIT.dataset] = st;
    put(pageState); pageState = null; TRAIT.mode = 'page';
    traitColMenuClose();
    document.getElementById('traitDrawer').remove();
    if (pageNodes){ const page = document.getElementById('page'); if (page) page.appendChild(pageNodes); pageNodes = null; }
    document.body.classList.remove('tdr-lock');
    document.removeEventListener('keydown', onKey);
    const cb = opts && opts.onApply; opts = null;
    if (apply && typeof cb === 'function') cb(chosen);
    if (returnFocus && document.body.contains(returnFocus) && returnFocus.focus) returnFocus.focus();
    returnFocus = null;
  }
  return {open, close, renderFoot, isOpen};
})();

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
  .trait-shell{display:grid;grid-template-columns:230px minmax(0,1fr);gap:16px;margin-top:14px;align-items:start}
  /* minmax(0,1fr), not 1fr: a 1fr track is at least as wide as its content, and a filtered table
     made the whole page 730-1,016 px wide on a 375 px phone */
  @media(max-width:820px){.trait-shell{grid-template-columns:minmax(0,1fr)}}
  .trait-shell>*{min-width:0}
  .trait-facets{border:1px solid var(--line);border-radius:12px;background:#fff;padding:10px 12px;position:sticky;top:10px;max-height:80vh;overflow:auto}
  .facet-head{display:flex;align-items:center;justify-content:space-between;font-weight:700;font-size:12.5px;margin-bottom:6px}
  .facet-block{margin:10px 0}
  .facet-block h4{margin:0 0 5px;font-size:11px;letter-spacing:.04em;text-transform:uppercase;color:var(--muted)}
  .facet-row{display:flex;align-items:center;gap:7px;font-size:12px;padding:2px 0;cursor:pointer}
  .facet-row .fv{flex:1;min-width:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
  .facet-row .fn{color:var(--faint);font-family:var(--mono);font-size:11px}
  .facet-all{display:flex;align-items:center;gap:7px;font-size:12px;font-weight:600;padding:2px 0 4px;margin-bottom:2px;border-bottom:1px solid var(--line);cursor:pointer}
  .facet-empty,.facet-block em{color:var(--faint);font-size:12px}
  .facet-logic{font-size:11.5px;line-height:1.45;color:var(--muted);background:#f6f8fc;border:1px solid var(--line);
    border-radius:8px;padding:7px 9px;margin:2px 0 6px}
  .facet-logic b{color:var(--ink);font-weight:600}
  .facet-any{float:right;text-transform:none;letter-spacing:0;font-weight:600;color:var(--blue-600,#2563eb)}
  .facet-any .link-btn{font-size:11px;margin-left:4px}
  .facet-row.zero .fv,.facet-row.zero .fn{color:var(--faint)}
  /* sidebar sections: Line details / GRIN phenotypes / Published phenotypes */
  .facet-sec{margin:16px 0 4px;padding:6px 0 4px;border-top:2px solid var(--line);font-weight:700;font-size:12.5px;color:var(--ink)}
  .facet-sec-n{font-weight:500;font-size:11px;color:var(--muted);margin-left:4px}
  /* same look as the .facet-block h4 headings (e.g. GRIN CONTINUOUS PHENOTYPES) */
  .facet-sub{margin:12px 0 2px;font-size:11px;font-weight:700;letter-spacing:.04em;text-transform:uppercase;color:var(--muted)}
  .facet-block.fb-col{margin:2px 0}
  .fb-col h4{display:flex;align-items:center;gap:4px;margin:0}
  .fb-tog{all:unset;cursor:pointer;display:flex;align-items:center;gap:4px;flex:1;min-width:0;padding:3px 0;
    font-size:11px;letter-spacing:.04em;text-transform:uppercase;color:var(--muted)}
  .fb-tog:hover{color:var(--ink)}
  .fb-tog:focus-visible{outline:2px solid var(--blue-600,#2563eb);outline-offset:1px}
  .fb-tog .caret{display:inline-flex;transition:transform .15s}
  .fb-tog .caret svg{width:11px;height:11px}
  .fb-col.op .fb-tog .caret{transform:rotate(90deg)}
  .fb-tog .fb-l{flex:1;min-width:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
  .fb-tog .fb-n{font-family:var(--mono);font-size:10.5px;text-transform:none;letter-spacing:0;color:var(--faint)}
  .fb-col.op{margin-bottom:8px}
  .trait-toolbar{display:flex;gap:10px;flex-wrap:wrap;align-items:center;margin-bottom:10px}
  .trait-search{flex:1;min-width:200px;display:flex;align-items:center;gap:7px;border:1px solid var(--line);border-radius:9px;padding:7px 10px;background:#fff}
  .trait-search svg{width:16px;height:16px;color:var(--muted)}
  .trait-search input{border:0;outline:0;flex:1;font-size:13px;background:transparent;color:var(--ink)}
  .trait-qbtns{display:flex;gap:6px;flex-wrap:wrap}
  .qbtn{border:1px solid var(--line);background:#fff;border-radius:7px;padding:6px 10px;font:600 12px/1 var(--body);color:var(--ink);cursor:pointer}
  .qbtn.solid{background:var(--blue-600,#2563eb);border-color:var(--blue-600,#2563eb);color:#fff}
  .trait-grid{display:flex;flex-direction:column;gap:8px}
  .ds-pill.nometa{opacity:.55}
  .trait-runbar-top{margin-top:12px}            /* under the dataset picker, above the filters and table */
  .tr-range .tr-sel{width:100%;font-size:12px;padding:5px;border:1px solid var(--line);border-radius:7px;background:#fff}
  .tr-range .tr-mm{display:flex;gap:5px;align-items:center;margin-top:6px}
  .tr-range .tr-mm input{width:70px;font-size:12px;padding:4px 6px;border:1px solid var(--line);border-radius:6px}
  .tr-range .tr-hint,.tr-range .tr-note{font-size:11px;color:var(--muted);margin-top:5px;line-height:1.35}
  .tr-chip{display:flex;justify-content:space-between;align-items:center;gap:6px;margin-top:6px;padding:4px 8px;border:1px solid #cfe3d6;background:#f3faf5;border-radius:8px;font-size:12px}
  .tg-table td.num{text-align:right;font-family:var(--mono)}
  .tg-table td.clip{max-width:22em;overflow:hidden;text-overflow:ellipsis}
  .tg{border:1px solid var(--line);border-radius:10px;overflow:hidden;background:#fff}
  .tg-head{display:flex;flex-wrap:wrap;align-items:center;gap:9px;padding:9px 11px;cursor:pointer;background:#fbfcfe}
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
  /* sideways scrolling: cells do not wrap, so a wide table scrolls inside its group instead of
     running off the page; the tick box and ID columns stay pinned at the left while it scrolls */
  .tg-scroll{overflow-x:auto;-webkit-overflow-scrolling:touch}
  .tg-table th,.tg-table td{white-space:nowrap}
  .tg-table .tcb{width:16px;min-width:16px;padding-left:10px;padding-right:6px}
  .tg-table .stk{position:sticky;z-index:1;background:#fff}
  .tg-table thead .stk{background:#fafbfd}
  .tg-table tbody tr:hover .stk{background:#f5f8ff}
  .tg-table tbody tr.on .stk{background:#eaf2ff}
  .tg-table .stk0{left:0}
  .tg-table .stk1{left:32px;box-shadow:inset -1px 0 0 #e3e8f0}
  .tg-table th.num{text-align:right}
  .th-sort{all:unset;cursor:pointer;display:inline-flex;align-items:center;gap:5px;font:inherit;color:inherit;
    text-transform:inherit;letter-spacing:inherit;border-radius:4px}
  .th-sort:hover{color:var(--ink)}
  .th-sort:focus-visible{outline:2px solid var(--blue-600,#2563eb);outline-offset:2px}
  .th-sort .th-ind{font-size:9px;opacity:.45}
  .th-sort.on{color:var(--blue-600,#2563eb)}
  .th-sort.on .th-ind{opacity:1}
  .th-wrap{display:inline-flex;align-items:center;gap:4px}
  .th-filt{border:1px solid transparent;background:none;border-radius:4px;width:18px;height:18px;padding:0;
    display:inline-flex;align-items:center;justify-content:center;font-size:10px;line-height:1;color:var(--muted);cursor:pointer}
  .th-filt:hover{border-color:var(--line);background:#fff;color:var(--ink)}
  .th-filt:focus-visible{outline:2px solid var(--blue-600,#2563eb);outline-offset:1px}
  .th-filt.on{color:#fff;background:var(--blue-600,#2563eb);border-color:var(--blue-600,#2563eb)}
  .th-filt svg{width:11px;height:11px}
  /* column menu (Excel-style sort and filter), above the drawer (z 61) */
  .tcm{position:fixed;z-index:70;width:260px;max-width:calc(100vw - 16px);display:flex;flex-direction:column;
    background:#fff;border:1px solid var(--line);border-radius:10px;box-shadow:var(--shadow-lg,0 10px 30px rgba(10,15,28,.18));
    padding:6px 0;font-size:12.5px;color:var(--ink);text-transform:none;letter-spacing:0}
  .tcm-head{padding:4px 12px 6px;font-weight:700;font-size:12px;color:var(--muted);white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
  .tcm-item{all:unset;display:block;padding:6px 12px;cursor:pointer;font-size:12.5px}
  .tcm-item:hover,.tcm-item:focus-visible{background:#f2f6fd}
  .tcm-item.on{color:var(--blue-600,#2563eb);font-weight:600}
  .tcm.trp{width:220px}
  .tcm-sep{height:1px;background:var(--line);margin:5px 0}
  .tcm-search{display:flex;align-items:center;gap:6px;margin:4px 12px 6px;border:1px solid var(--line);border-radius:7px;padding:5px 8px}
  .tcm-search svg{width:14px;height:14px;color:var(--muted);flex:0 0 auto}
  .tcm-search input{border:0;outline:0;flex:1;min-width:0;font-size:12.5px;background:transparent;color:var(--ink)}
  .tcm-row{display:flex;align-items:center;gap:7px;padding:3px 12px;cursor:pointer}
  .tcm-row:hover{background:#f7f9fd}
  .tcm-row .tv{flex:1;min-width:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
  .tcm-row .tn{color:var(--faint);font-family:var(--mono);font-size:11px}
  .tcm-row.zero .tv,.tcm-row.zero .tn{color:var(--faint)}
  .tcm-all{font-weight:600}
  .tcm-list{flex:1 1 auto;min-height:60px;max-height:260px;overflow:auto;border-top:1px solid #eef1f6;border-bottom:1px solid #eef1f6;margin:2px 0}
  .tcm-none{padding:10px 12px;color:var(--faint)}
  .tcm-lbl{padding:4px 12px;color:var(--muted)}
  .tcm-mm{display:flex;gap:6px;align-items:center;padding:2px 12px}
  .tcm-mm input{width:90px;font-size:12px;padding:4px 6px;border:1px solid var(--line);border-radius:6px}
  .tcm-hint{padding:5px 12px;font-size:11px;color:var(--muted);line-height:1.35}
  .tcm-foot{display:flex;align-items:center;gap:6px;padding:8px 12px 4px}
  .tcm-foot .tcm-sp{flex:1}
  .tcm-foot .link-btn[disabled]{opacity:.4;cursor:default}
  .tcm-foot .qbtn[disabled]{opacity:.5;cursor:default}
  /* column picker */
  .tcm.tcp{width:330px}
  .tcp-list{max-height:380px}
  .tcp-sec{padding:8px 12px 3px;font-size:10.5px;font-weight:700;letter-spacing:.04em;text-transform:uppercase;color:var(--muted);
    background:#fbfcfe;position:sticky;top:0;z-index:1}
  .tcp-list .tcm-row input[disabled]+.tv{color:var(--muted)}
  .tcp-note{padding:2px 12px 5px;font-size:11px;line-height:1.35;color:var(--muted);word-break:break-word}
  .trait-qbtns .qbtn svg{width:13px;height:13px;vertical-align:-2px}
  /* download dialog */
  .tdl-set{border:0;margin:0;padding:0;display:flex;flex-direction:column;gap:7px}
  .tdl-set legend{font-size:11px;font-weight:700;letter-spacing:.04em;text-transform:uppercase;color:var(--muted);padding:0;margin-bottom:6px}
  .tdl-fmt{display:flex;gap:18px}
  .tdl-sub{margin-top:6px;padding-top:9px;border-top:1px solid var(--line)}   /* applies to any column choice, so not indented under one */
  .trait-send-check.dis{opacity:.5;cursor:default}
  .tg-table .cb{display:inline-flex;width:15px;height:15px;border:1.5px solid #c3ccda;border-radius:4px;color:#fff}
  .tg-table .cb svg{width:12px;height:12px;opacity:0}
  .tg-table tr.on .cb{background:var(--blue-600);border-color:var(--blue-600)}
  .tg-table tr.on .cb svg{opacity:1}
  .trait-empty{padding:26px;text-align:center;color:var(--faint);border:1px dashed var(--line);border-radius:10px}
  .qbtn.on{background:#eef4ff;border-color:var(--blue-600,#2563eb);color:var(--blue-600,#2563eb)}
  /* drawer: SNPTrait over SNPVersity (TraitDrawer) */
  body.tdr-lock{overflow:hidden}
  .tdr-backdrop{position:fixed;inset:0;background:rgba(10,15,28,.38);z-index:60}
  .tdr{position:fixed;top:0;right:0;bottom:0;width:min(1180px,calc(100vw - 56px));background:var(--paper);z-index:61;
    display:flex;flex-direction:column;box-shadow:var(--shadow-lg);animation:tdrIn .18s ease-out}
  @keyframes tdrIn{from{transform:translateX(28px);opacity:.5}to{transform:none;opacity:1}}
  @media (prefers-reduced-motion:reduce){.tdr{animation:none}}
  .tdr-head{display:flex;align-items:flex-start;justify-content:space-between;gap:14px;padding:16px 22px 13px;background:#fff;border-bottom:1px solid var(--line)}
  .tdr-head h2{font-family:var(--disp);font-size:18px;font-weight:600;margin:0;letter-spacing:-.2px}
  .tdr-sub{font-size:12.5px;color:var(--muted);margin-top:3px;max-width:75ch}
  .tdr-x{border:0;background:none;font-size:26px;line-height:1;color:var(--muted);cursor:pointer;padding:0 4px}
  .tdr-x:hover{color:var(--ink)}
  .tdr-body{flex:1;overflow:auto;padding:14px 22px 18px}
  .tdr-body .trait-shell{margin-top:4px}
  .tdr-body .trait-facets{top:0}
  .tdr-foot{display:flex;align-items:center;gap:14px;flex-wrap:wrap;padding:12px 22px;background:#fff;border-top:1px solid var(--line)}
  .tdr-summ{flex:1;min-width:220px;font-size:13px;color:var(--ink)}
  .tdr-n{font-family:var(--mono);font-size:18px;font-weight:700}
  .tdr-delta{color:var(--muted);margin-left:8px}
  .tdr-add{color:#176c3a} .tdr-rem{color:#b42318}
  .tdr-status{font-size:11.5px;color:var(--faint);margin-top:2px}
  .tdr-acts{display:flex;gap:8px}
  @media (max-width:700px){ .tdr{left:0;width:auto} .tdr-head,.tdr-body,.tdr-foot{padding-left:14px;padding-right:14px} }
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
