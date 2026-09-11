/* =====================================================================
 *  snpgeo.js — SNPGeo: geographic distribution of SNPs.
 *
 *  Maps the geographic distribution of SNPs across countries and sub-regions.
 *  Shows the percentage of isolates carrying each variant in every country,
 *  with interactive hover stats and regional breakdowns.
 *
 *  Receives variant data from SNPVersity via S.geoInput:
 *    rows: array of variant objects (variant, consequence, residue, scores)
 *    accs: array of accession objects (id, label, gt columns)
 *    chr, start, end, dataset, datasetName
 *
 *  Geographic metadata is pulled from:
 *    window.SNPGEO_REGIONS[isolateId] = {country, iso3, state, county, ...}
 *    window.SNPGEO_COUNTRY_ISO = {countryName -> iso3}
 *    data/geo/countries.geo.json = GeoJSON of world countries
 *
 *  Depends on d3.v5 (loaded for PanEffect), core.js (S, ICONS), snphandoff.js
 * ===================================================================== */

const GEO = {
  dataset: null,
  input: null,           // S.geoInput snapshot
  rows: [],              // variant rows
  accs: [],              // accession array
  snpIndex: 0,           // currently active SNP row index
  selectedCountry: null, // user-selected country (iso3)
  snpStats: {},          // snpIndex -> { countryIso3: {count, total, pct, called, alt, ref, altFrac, isolates, states} }
  geoJson: null,         // loaded countries.geo.json
  mapSvg: null,          // d3 selection of map container
  colorMode: 'freq',     // 'freq' | 'refalt' — see GEO_MODES
};

/* ================= HELPERS ================= */
function escGeo(s){
  return String(s==null?'':s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
}
function escGeoAttr(s){
  return escGeo(s).replace(/"/g,'&quot;').replace(/'/g,'&#39;');
}

/* ================= PAGE RENDER ================= */
/* ---- dataset chooser: mirrors SNPFold's datasetChooser()/datasetCardsHTML() exactly ---- */
function geoDatasetCardsHTML(){
  const list = (typeof Data !== 'undefined' && typeof Data.datasets === 'function') ? Data.datasets() : [];
  return list.map(d => {
    const id = d.id;
    const sel = S.dataset === id;
    const name = d.name || id;
    const sub = d.sub || '';
    return `<button type="button" class="geo-ds ${sel?'sel':''}" onclick="S.dataset='${escGeoAttr(id)}';geoInitDataset()">
      <span class="geo-ds-dot"></span>
      <span class="geo-ds-txt"><span class="geo-ds-name">${escGeo(name)}</span>${sub?`<span class="geo-ds-sub">${escGeo(sub)}</span>`:''}</span>
    </button>`;
  }).join('');
}
function geoDatasetChooser(){
  const cards = geoDatasetCardsHTML();
  if (!cards) return '';
  return `<div class="card pad geo-ds-card" style="margin-bottom:16px">
    <div class="geo-ds-head">Dataset</div>
    <div class="geo-ds-grid" id="geoDsGrid">${cards}</div>
  </div>`;
}

/* Pre-filled into the gene box so the page is runnable with one click on
   Search. It is only a default: any gene the user has already looked up wins,
   and the empty-state focus selects the text so typing replaces it outright. */
const GEO_DEFAULT_GENE = 'FGSG_00025';

/* ---- search bar: mirrors SNPFold's searchBar() exactly, minus the model-source picker ---- */
function geoSearchBar(){
  const gene = (S.geoInput && S.geoInput.gene) || GEO_DEFAULT_GENE;
  return `<div class="card pad" style="margin-bottom:16px;display:flex;gap:10px;align-items:center;flex-wrap:wrap">
    <span style="font-size:10.5px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.4px">Gene model</span>
    <input id="geoGeneInput" value="${escGeoAttr(gene)}" spellcheck="false"
      style="flex:1;min-width:280px;border:1px solid var(--line);border-radius:9px;padding:9px 11px;font-family:var(--mono);font-size:13px"
      onkeydown="if(event.key==='Enter')geoLookupGene()"
      placeholder="e.g. FGSG_00777">
    <button class="btn" onclick="geoLookupGene()">Search</button>
  </div>`;
}

function renderGeo(page){
  const hasResults = S.geoInput && S.geoInput.rows && S.geoInput.rows.length;
  page.className = 'page fade';
  injectGeoCSS();

  const header = geoDatasetChooser() + geoSearchBar()
    + `<div id="geoLookupStatus" style="margin:-8px 0 16px;font-size:12px;color:var(--muted)"></div>`;

  const results = hasResults ? `
    <div class="card pad" style="margin-bottom:16px">
      <div class="n" style="margin:0">${S.geoInput.gene ? 'Gene: ' + escGeo(S.geoInput.gene) : 'Region: chr' + S.geoInput.chr + ':' + S.geoInput.start.toLocaleString() + '-' + S.geoInput.end.toLocaleString()}</div>
      <div style="color:var(--muted);font-size:13px;margin-top:4px">${S.geoInput.rows.length} variants</div>
    </div>
    <div class="geo-shell">
      <div class="geo-map-wrap">
        <div class="geo-map-tools">
          <select id="geoColorMode" class="geo-mode-sel" onchange="geoSetColorMode(this.value)"
                  data-tt="Choose what the country colours encode.">${geoModeOptionsHTML()}</select>
          <button class="geo-tool-btn" onclick="geoExportPNG()" data-tt="Download the map as a PNG image (2x resolution).">PNG</button>
          <button class="geo-tool-btn" onclick="geoExportSVG()" data-tt="Download the map as a vector SVG.">SVG</button>
        </div>
        <div id="geoMap" class="geo-map"></div>
        <div id="geoLegend" class="geo-legend"></div>
      </div>
      <div id="geoDetail" class="geo-detail"></div>
    </div>
    <div class="geo-table" id="geoTableContainer">
      <div class="geo-hint" id="geoKeyHint">
        Click a row, then use <kbd>&uarr;</kbd><kbd>&darr;</kbd> to step through the variants
      </div>
      <div id="geoTable"></div>
    </div>
  ` : `
    <div class="card pad" style="text-align:center;color:var(--muted);padding:48px">
      Enter a gene name above, or load data from SNPVersity
    </div>
  `;

  page.innerHTML = header + results;

  if (!hasResults) {
    // Focus gene input when no results yet
    setTimeout(() => {
      const inp = document.getElementById('geoGeneInput');
      // select() so the pre-filled default is replaced by whatever they type.
      if (inp){ inp.focus(); inp.select(); }
    }, 100);
    return;
  }

  // Results loaded: render map + table in main area
  geoLoad();

  // Load GeoJSON and render
  loadCountriesGeoJSON(() => {
    geoRenderTable();
    geoRenderMap();
    geoRenderOverview();   // right pane opens on the country overview
  });
}

// Register the tool immediately after renderGeo is defined
if (typeof SNPTools !== 'undefined' && SNPTools.register) {
  SNPTools.register('snpgeo', { render: renderGeo });
}

/* ================= GENE LOOKUP ================= */
async function geoLookupGene(){
  const inp = document.getElementById('geoGeneInput');
  const gene = inp ? inp.value.trim() : '';
  const statusEl = document.getElementById('geoLookupStatus');
  
  if (!gene){
    if (statusEl) statusEl.innerHTML = '<span style="color:#f87171">Please enter a gene name</span>';
    return;
  }
  
  if (statusEl) statusEl.innerHTML = '<span style="color:var(--muted)">Searching...</span>';
  
  try {
    // Query the backend for gene coordinates and variants.
    // IMPORTANT: unlike SNPFold (which only needs a representative sample to map
    // onto one structure), SNPGeo needs EVERY isolate in the dataset queried —
    // the map's per-country percentages are meaningless over a partial sample.
    // Data.queryVariantsByGene() defaults to a ~12-accession sample when `ids` is
    // omitted, so pass the full accession list explicitly (same pattern SNPFunction
    // uses for whole-dataset queries).
    const allIds = (Data.accessionsFor(S.dataset) || []).map(a => a.id);
    const result = await Data.queryVariantsByGene(S.dataset, gene, allIds);
    
    if (!result || !result.rows || result.rows.length === 0){
      if (statusEl) statusEl.innerHTML = '<span style="color:#f87171">No variants found for this gene</span>';
      return;
    }
    
    // Populate S.geoInput with the results
    S.geoInput = {
      rows: result.rows,
      accs: result.accs || [],
      chr: result.chr || 'unknown',
      start: result.start || 0,
      end: result.end || 0,
      dataset: S.dataset,
      datasetName: (Data.datasets().find(d=>d.id===S.dataset)||{}).name || S.dataset,
      gene: gene,  // track the gene name
      vcfUrl: result.vcfUrl
    };
    
    if (statusEl) statusEl.innerHTML = `<span style="color:#10b981">${result.rows.length} variants found</span>`;
    
    // Re-render SNPGeo with results
    setTimeout(() => {
      const page = document.querySelector('.page');
      if (page) renderGeo(page);
    }, 200);
    
  } catch (err){
    console.error('[SNPGeo] Gene lookup error:', err);
    if (statusEl) statusEl.innerHTML = '<span style="color:#f87171">Error: ' + escGeo(err.message) + '</span>';
  }
}

/* ================= DATASET SWITCHING ================= */
function geoInitDataset(){
  // Called when user changes dataset in sidebar
  // Reset previous results and re-render
  S.geoInput = null;
  const page = document.querySelector('.page');
  if (page) {
    renderGeo(page);
  }
}

/* ================= DATA LOAD & PREP ================= */
function geoLoad(){
  GEO.input = S.geoInput;
  GEO.rows = S.geoInput.rows || [];
  GEO.accs = S.geoInput.accs || [];
  GEO.dataset = S.geoInput.dataset;
  GEO.snpIndex = 0;
  GEO.selectedCountry = null;
  GEO.snpStats = {};
  
  // Pre-compute stats for all SNPs (lazily computed on first access)
}

/* Classify one genotype call. Fusarium is HAPLOID — a sample column is a
   single allele index ("0", "1", ".") — but the diploid spellings are
   accepted too so this stays correct if the store ever emits them.

   An uncalled genotype is MISSING, not an alt call. The previous test
   (`gt && gt!=='0' && gt!=='0|0' && gt!=='0/0'`) scored "." as carrying the
   variant, so every isolate with no call at a site inflated that country's
   carrier count, its map colour and its region table. */
function geoGtClass(gt){
  if (gt == null) return 'missing';
  const s = String(gt).trim();
  if (s === '' || s === '.' || s === './.' || s === '.|.') return 'missing';
  if (s === '0' || s === '0/0' || s === '0|0') return 'ref';
  return 'alt';
}

/* Per-country statistics for one variant.

   Two denominators, deliberately kept apart:
     total  — every isolate the geographic metadata places in the country,
              queried or not. Drives the "x/y" and % in the overview table
              and the frequency map (a country's variant frequency is over
              its whole isolate collection).
     called — queried isolates in the country with a non-missing genotype at
              this site. Drives the reference/alternative view, where the
              question is what the *observed* alleles look like; a country
              with called === 0 has no answer and is drawn as "no calls"
              rather than being painted at the 100%-reference end. */
function aggregateGeoData(snpIndex){
  // If cached, return it
  if (GEO.snpStats[snpIndex]) return GEO.snpStats[snpIndex];
  
  const row = GEO.rows[snpIndex];
  if (!row) return {};
  
  const result = {};
  const REGIONS = window.SNPGEO_REGIONS || {};
  
  // STEP 1: Count ALL isolates per country (from geographic metadata)
  Object.entries(REGIONS).forEach(([accId, geo]) => {
    const iso3 = geo.iso3 || '???';
    const state = geo.state || geo.county || 'Unknown';
    
    if (!result[iso3]){
      result[iso3] = {
        country: geo.country || 'Unknown',
        iso3: iso3,
        count: 0,       // isolates carrying the alt allele
        total: 0,       // all isolates known from this country
        pct: 0,         // count / total
        called: 0,      // queried isolates with a genotype call here
        alt: 0,         // = count, named for the ref/alt view
        ref: 0,         // called isolates carrying the reference allele
        missing: 0,     // queried isolates with no call here
        altFrac: null,  // alt / called, or null when called === 0
        isolates: [],
        states: {}
      };
    }
    
    // Count total isolates in this country (regardless of whether they were queried)
    result[iso3].total += 1;
    
    if (!result[iso3].states[state]){
      result[iso3].states[state] = {count: 0, total: 0, called: 0, alt: 0, ref: 0, isolates: []};
    }
    result[iso3].states[state].total += 1;
  });
  
  // STEP 2: Classify every queried isolate's call at this site
  GEO.accs.forEach((acc, accIdx) => {
    const accId = String(acc.id);
    const geo = REGIONS[accId];
    if (!geo) return; // unmapped isolate
    
    const iso3 = geo.iso3 || '???';
    const state = geo.state || geo.county || 'Unknown';
    const c = result[iso3];
    if (!c) return;
    const st = c.states[state];
    
    // Genotypes may be in row.gts (array, one entry per accession) or in
    // flat row.gt_N columns.
    const gt = (row.gts && Array.isArray(row.gts)) ? row.gts[accIdx] : row[`gt_${accIdx}`];
    const cls = geoGtClass(gt);
    
    if (cls === 'missing'){ c.missing += 1; return; }
    
    c.called += 1;
    if (st) st.called += 1;
    
    if (cls === 'alt'){
      c.count += 1;
      c.alt += 1;
      c.isolates.push(accId);
      if (st){ st.count += 1; st.alt += 1; st.isolates.push(accId); }
    } else {
      c.ref += 1;
      if (st) st.ref += 1;
    }
  });
  
  // Derived fractions
  Object.values(result).forEach(r => {
    r.pct = r.total > 0 ? (r.count / r.total) : 0;
    r.altFrac = r.called > 0 ? (r.alt / r.called) : null;
    Object.values(r.states).forEach(s => {
      s.altFrac = s.called > 0 ? (s.alt / s.called) : null;
    });
  });
  
  GEO.snpStats[snpIndex] = result;
  return result;
}

/* ================= TABLE RENDERING ================= */
/* Residue position + "R234L"-style residue-substitution label, both derived
   from row.sub (the raw INFO SUB token, e.g. "T209A"/"K90R"/".") via
   Data.parseSub. row.resi/row.variant are NOT fields queryVariants ever
   produces — only rankImpact()'s output (a different code path, used by
   SNPImpact) has them, so reading row.resi directly always fell back to
   "—" here regardless of consequence. Mirrors Data's own (unexported)
   hgvsProtein() special-casing for frameshift/stop/start-loss so the label
   matches what SNPFold would show for the same variant. */
function geoResiduePair(row){
  const p = Data.parseSub ? Data.parseSub(row.sub) : null;
  if (!p || p.resi == null) return {resi: null, label: null};
  const cls = Data.classifyConsequence ? Data.classifyConsequence(row.effect) : null;
  const ref = p.ref || '';
  let label;
  if (cls){
    switch (cls.label){
      case 'Frameshift':  label = ref + p.resi + 'fs'; break;
      case 'Stop gained': label = ref + p.resi + '*'; break;
      case 'Start lost':  label = (ref || 'M') + p.resi + '?'; break;
      case 'Stop lost':   label = '*' + p.resi + (p.alt && p.alt !== '*' ? p.alt : 'ext'); break;
      default:
        label = cls.klass === 'indel'
          ? ref + p.resi + (cls.label.indexOf('insertion') >= 0 ? 'ins' : 'del')
          : ref + p.resi + (p.alt || '');
    }
  } else {
    // Non-structural (synonymous, intron, UTR, …) or unrecognized effect —
    // still show a plain ref+resi+alt label when parseSub gave us one.
    label = (ref || p.alt) ? ref + p.resi + (p.alt || '') : null;
  }
  return {resi: p.resi, label};
}

function geoRenderTable(){
  const tbl = document.getElementById('geoTable');
  if (!tbl) return;
  
  const rows = GEO.rows;
  if (!rows.length){
    tbl.innerHTML = '<div class="geo-empty">No variants to display.</div>';
    return;
  }
  
  // Score columns share SNPFold's colour ramp and pill form via window.ScorePill
  // (defined in snpfold.js), so the SNPFold and SNPGeo tables render scores
  // identically. ScorePill.cell() returns the inner <span>; wrap it in the
  // numeric <td> this table expects. Fallback keeps the column readable if
  // snpfold.js somehow hasn't loaded.
  const models = Data.scoreModels(GEO.dataset);
  const pill = (v) => (window.ScorePill ? ScorePill.cell(v) : (v == null ? '—' : v));
  const sc = (v) => `<td class="num">${pill(v)}</td>`;

  const thead = `
    <tr>
      <th class="geo-cb"></th>
      <th data-tt="Variant identifier — protein substitution when coding, else position and alleles.">Variant</th>
      <th data-tt="Predicted consequence — e.g. missense, synonymous, intron, frameshift.">Consequence</th>
      <th class="num" data-tt="Affected residue position in the protein, when applicable.">Residue</th>
      <th class="num" data-tt="Residue-level substitution — wild-type residue, position, mutant residue (e.g. R234L).">Residue Variant</th>
      ${models.map(m=>`<th class="num ${m.kind==='protein'?'lm-prot':'lm-dna'}" data-tt="${escGeoAttr(m.tip)}">${m.label}</th>`).join('')}
    </tr>`;
  
  const tbody = rows.map((row, idx) => {
    const active = idx === GEO.snpIndex ? 'active' : '';
    // Construct variant string from pos/ref/alt if not already present
    const variant = row.variant || (row.pos && row.ref && row.alt ? `${row.pos} ${row.ref}>${row.alt}` : '—');
    const impact = row.impact ? String(row.impact).toLowerCase() : '';
    const consequence = row.consequence || row.consClass || row.effect || '—';
    const {resi, label: resiVariant} = geoResiduePair(row);
    return `
      <tr class="${active}" onclick="geoSelectSnp(${idx})">
        <td class="geo-cb"><input type="radio" name="geoVariantSel" class="geo-radio"
          tabindex="-1" ${active ? 'checked' : ''}
          aria-label="Show ${escGeoAttr(variant)} on the map"></td>
        <td class="c-mono">${escGeo(variant)}</td>
        <td class="effect-cell">${escGeo(consequence)}${impact?` <span class="pill ${impact}">${escGeo(row.impact)}</span>`:''}</td>
        <td class="num">${resi != null ? resi : '—'}</td>
        <td class="c-mono">${resiVariant ? escGeo(resiVariant) : '—'}</td>
        ${models.map(m=>sc(row[m.key])).join('')}
      </tr>`;
  }).join('');
  
  tbl.innerHTML = `<div class="tbl-wrap"><table class="vcf geo-vcf">
    <thead>${thead}</thead>
    <tbody>${tbody}</tbody>
  </table></div>`;
  geoFitTableHeight();
  attachTT();
}

/* Size the variant list's scrollport to whatever viewport is left underneath
   it, so the list scrolls inside its own box and the map above never moves.
   The shared .tbl-wrap default is 72vh, which is taller than the space left
   once the dataset chooser, search bar and 400px map are stacked above it —
   so the page itself scrolled and took the map off-screen with it.
   Measured rather than hardcoded: the chooser and search bar re-wrap at
   different window widths, so the offset above the table is not a constant. */
const GEO_VISIBLE_ROWS = 10;   // variants shown at once; the rest are scrolled to
const GEO_ROW_H_FALLBACK = 29; // used only where layout is unavailable (headless)
const GEO_HEAD_H_FALLBACK = 31;

function geoFitTableHeight(){
  const wrap = document.querySelector('#geoTable .tbl-wrap');
  if (!wrap) return;
  const head = wrap.querySelector('thead');
  const row = wrap.querySelector('tbody tr');
  if (!row) return;
  // Measured from a real row rather than assumed, so the frozen height stays
  // exactly GEO_VISIBLE_ROWS deep at any browser zoom or font size.
  const rowH = (row.getBoundingClientRect && row.getBoundingClientRect().height) || GEO_ROW_H_FALLBACK;
  const headH = (head && head.getBoundingClientRect && head.getBoundingClientRect().height) || GEO_HEAD_H_FALLBACK;
  wrap.style.maxHeight = Math.round(headH + rowH * GEO_VISIBLE_ROWS) + 'px';
}

function geoSelectSnp(idx){
  GEO.snpIndex = idx;
  GEO.selectedCountry = null;
  geoRenderTable();
  geoRenderMap();
  geoRenderOverview();   // a different SNP means different numbers — back to the overview
}

/* ================= KEYBOARD NAVIGATION =================
   Up/Down step through the variant list, so a whole gene can be scanned
   without going back to the mouse for every row. */

/* Bring the active row into view by moving ONLY the table's own scrollport.
   scrollIntoView() was wrong here: it walks every scrollable ancestor up to
   the document, so stepping past the last visible row scrolled the whole page
   and pushed the map off the top. Writing wrap.scrollTop directly cannot move
   anything outside the table, which is the guarantee we actually want. */
function geoScrollActiveIntoView(){
  const wrap = document.querySelector('#geoTable .tbl-wrap');
  const tr = wrap && wrap.querySelector('tr.active');
  if (!wrap || !tr || !wrap.getBoundingClientRect) return;
  const wr = wrap.getBoundingClientRect();
  const rr = tr.getBoundingClientRect();
  // The thead is sticky, so it overlays the top of the scrollport: a row
  // scrolled flush to wr.top would sit underneath it.
  const head = wrap.querySelector('thead');
  const headH = (head && head.getBoundingClientRect().height) || 0;
  if (rr.top < wr.top + headH)      wrap.scrollTop -= (wr.top + headH - rr.top);
  else if (rr.bottom > wr.bottom)   wrap.scrollTop += (rr.bottom - wr.bottom);
}

function geoKeyNav(ev){
  if (ev.key !== 'ArrowDown' && ev.key !== 'ArrowUp') return;
  if (ev.metaKey || ev.ctrlKey || ev.altKey || ev.shiftKey) return;
  // Leave the arrows alone while the user is in the gene box, a dropdown, or
  // any other text entry — there they mean caret/option movement.
  const t = ev.target;
  if (t && (/^(INPUT|TEXTAREA|SELECT)$/.test(t.tagName) || t.isContentEditable)) return;
  // #geoTable only exists while the SNPGeo page is the one rendered, so its
  // presence doubles as the "are we on this page" test.
  const tbl = document.getElementById('geoTable');
  if (!tbl) return;
  const n = (GEO.rows || []).length;
  if (!n) return;
  const next = ev.key === 'ArrowDown'
    ? Math.min(GEO.snpIndex + 1, n - 1)
    : Math.max(GEO.snpIndex - 1, 0);
  // Swallow the key even at the ends, otherwise the page scrolls instead and
  // the list appears to jump away from the user.
  ev.preventDefault();
  if (next === GEO.snpIndex) return;
  geoSelectSnp(next);
  geoScrollActiveIntoView();
}

/* Bound once on the document: geoRenderTable() rewrites its own innerHTML on
   every selection, so a listener attached to the table would not survive the
   first keypress. */
function geoBindKeyNav(){
  if (typeof document === 'undefined' || geoBindKeyNav.done) return;
  document.addEventListener('keydown', geoKeyNav);
  // The fitted table height is viewport-relative, so it has to be recomputed
  // when the window changes. Debounced — resize fires continuously on drag.
  if (typeof window !== 'undefined' && window.addEventListener){
    let t = null;
    window.addEventListener('resize', () => {
      clearTimeout(t);
      t = setTimeout(geoFitTableHeight, 120);
    });
  }
  geoBindKeyNav.done = true;
}
geoBindKeyNav();

/* ================= COLOUR MODES =================
   Two ways to paint the choropleth, chosen from the dropdown over the map.

   freq   — sequential, "what fraction of this country's isolates carry the
            variant". Denominator is every isolate known from the country.
   refalt — diverging, "which allele does this country look like". A football
            field: 100% reference at one end, 50/50 in the middle, 100%
            alternative at the other, over the isolates actually called here.
            The palette is Okabe-Ito blue -> neutral -> orange/vermillion,
            which stays separable under deuteranopia, protanopia and
            tritanopia (red/green diverging schemes do not).

   Countries with no call at this site are painted GEO_NODATA in both modes
   rather than at a scale endpoint — "nothing observed" is not "0% carriers"
   and certainly not "100% reference". */
const GEO_NODATA = '#e4e8ee';

const GEO_MODES = {
  freq: {
    label: 'Variant frequency',
    legendTitle: '% carrying variant',
    stops: [[0, '#f0fdf4'], [1, '#059669']],
    ends: ['0%', '50%', '100%'],
    value: st => (st && st.called > 0) ? st.pct : null,
  },
  refalt: {
    label: 'Reference ↔ Alternative allele',
    legendTitle: 'Allele composition (called isolates)',
    stops: [
      [0.00, '#0072B2'],   // Okabe-Ito blue        — all reference
      [0.25, '#56B4E9'],   // Okabe-Ito sky blue
      [0.50, '#F2F2F2'],   // neutral               — 50 / 50
      [0.75, '#E69F00'],   // Okabe-Ito orange
      [1.00, '#D55E00'],   // Okabe-Ito vermillion  — all alternative
    ],
    ends: ['100% ref', '50 / 50', '100% alt'],
    value: st => (st && st.altFrac != null) ? st.altFrac : null,
  },
};
function geoMode(){ return GEO_MODES[GEO.colorMode] || GEO_MODES.freq; }

function geoHexRGB(h){
  const s = h.replace('#','');
  const n = s.length === 3 ? s.split('').map(c => c + c).join('') : s;
  return [parseInt(n.slice(0,2),16), parseInt(n.slice(2,4),16), parseInt(n.slice(4,6),16)];
}
/* Piecewise-linear interpolation across a mode's stops, t in [0,1]. */
function geoRamp(stops, t){
  const x = Math.max(0, Math.min(1, t));
  for (let i = 1; i < stops.length; i++){
    const a = stops[i-1], b = stops[i];
    if (x <= b[0]){
      const u = (b[0] === a[0]) ? 0 : (x - a[0]) / (b[0] - a[0]);
      const ca = geoHexRGB(a[1]), cb = geoHexRGB(b[1]);
      return `rgb(${[0,1,2].map(k => Math.round(ca[k] + u * (cb[k] - ca[k]))).join(',')})`;
    }
  }
  return stops[stops.length-1][1];
}
function geoFill(stat){
  const m = geoMode();
  const v = m.value(stat);
  return (v == null) ? GEO_NODATA : geoRamp(m.stops, v);
}
function geoSetColorMode(mode){
  if (!GEO_MODES[mode]) return;
  GEO.colorMode = mode;
  geoRenderMap();
}
function geoModeOptionsHTML(){
  return Object.entries(GEO_MODES).map(([k, m]) =>
    `<option value="${k}"${GEO.colorMode === k ? ' selected' : ''}>${escGeo(m.label)}</option>`).join('');
}

/* Label for the variant currently driving the map. */
function geoVariantLabel(row){
  if (!row) return '';
  if (row.variant) return String(row.variant);
  const {label} = geoResiduePair(row);
  if (label) return label;
  return (row.pos && row.ref && row.alt) ? `${row.pos} ${row.ref}>${row.alt}` : 'variant ' + (GEO.snpIndex + 1);
}
function geoScopeLabel(){
  const i = GEO.input || {};
  if (i.gene) return String(i.gene);
  if (i.chr) return `chr${i.chr}:${i.start}-${i.end}`;
  return 'region';
}

/* Tooltip body — mode-aware, so the number the colour encodes is the number
   the tooltip leads with. */
function geoTipHTML(stat){
  const pct = (stat.pct * 100).toFixed(1);
  let body = `${stat.count}/${stat.total} isolates carry the variant (${pct}%)`;
  if (stat.called > 0){
    const af = (stat.altFrac * 100).toFixed(1);
    body += `<br/>Called: ${stat.called} · ref ${stat.ref} / alt ${stat.alt} (${af}% alt)`;
  } else {
    body += `<br/><span style="opacity:.7">No genotype calls at this site</span>`;
  }
  return `<b>${escGeo(stat.country)}</b><br/>${body}`;
}

/* ================= MAP RENDERING ================= */
function geoRenderMap(){
  const mapDiv = document.getElementById('geoMap');
  if (!mapDiv) return;
  
  const stats = aggregateGeoData(GEO.snpIndex);
  
  if (!GEO.geoJson){
    mapDiv.innerHTML = '<div class="geo-loading">Loading map...</div>';
    return;
  }
  
  mapDiv.innerHTML = '';
  const width = mapDiv.clientWidth || 600;
  const height = mapDiv.clientHeight || 400;
  
  // Create SVG
  const svg = d3.select(mapDiv).append('svg')
    .attr('width', width)
    .attr('height', height)
    .attr('viewBox', `0 0 ${width} ${height}`);
  GEO.mapSvg = svg;
  
  // Simple Mercator projection
  const projection = d3.geoMercator()
    .fitSize([width, height], GEO.geoJson);
  
  const pathGenerator = d3.geoPath().projection(projection);

  // Draw countries
  svg.selectAll('path.country')
    .data(GEO.geoJson.features)
    .enter()
    .append('path')
    .attr('class', 'country')
    .attr('d', pathGenerator)
    .style('fill', d => geoFill(stats[d.properties.iso3]))
    .style('stroke', '#ddd')
    .style('stroke-width', '0.5px')
    .style('cursor', 'pointer')
    .on('mouseenter', function(d){
      const stat = stats[d.properties.iso3];
      if (!stat) return;
      d3.select(this).style('stroke', '#000').style('stroke-width', '1.5px');
      const tt = document.getElementById('tt');
      if (tt){
        tt.innerHTML = geoTipHTML(stat);
        tt.classList.add('show');
      }
    })
    .on('mousemove', function(){
      /* d3 v5 passes (datum, index, nodes) to listeners — the event is on
         d3.event. The old `function(e)` read clientX off the datum, so the
         tooltip never tracked the pointer. */
      const ev = d3.event;
      const tt = document.getElementById('tt');
      if (tt && ev){
        tt.style.left = (ev.clientX + 12) + 'px';
        tt.style.top = (ev.clientY + 14) + 'px';
      }
    })
    .on('mouseleave', function(){
      d3.select(this).style('stroke', '#ddd').style('stroke-width', '0.5px');
      const tt = document.getElementById('tt');
      if (tt) tt.classList.remove('show');
    })
    .on('click', function(d){
      const iso3 = d.properties.iso3;
      if (!stats[iso3]) return;
      GEO.selectedCountry = iso3;
      geoShowDetail(iso3);
    });
  
  // Legend
  geoRenderLegend();
}

function geoRenderLegend(){
  const leg = document.getElementById('geoLegend');
  if (!leg) return;
  const m = geoMode();
  // Continuous bar rather than discrete swatches: the scale is continuous,
  // and for refalt the football-field reading depends on seeing the middle.
  const grad = m.stops.map(s => `${s[1]} ${(s[0]*100).toFixed(0)}%`).join(',');
  leg.innerHTML = `
    <div class="geo-leg-title">${escGeo(m.legendTitle)}</div>
    <div class="geo-leg-bar" style="background:linear-gradient(to right,${grad})"></div>
    <div class="geo-leg-labels">
      <span>${escGeo(m.ends[0])}</span>
      <span>${escGeo(m.ends[1])}</span>
      <span>${escGeo(m.ends[2])}</span>
    </div>
    <div class="geo-leg-nd"><span class="geo-leg-nd-sw"></span>No calls</div>`;
}

/* ================= MAP EXPORT (PNG / SVG) =================
   Builds a standalone figure from the live map: title block, the country
   paths (their fills are already inline style attributes, so nothing depends
   on the page stylesheet), and a legend drawn in SVG rather than the HTML
   overlay the page uses. PNG is that same SVG rasterised through a canvas at
   2x. No external references, so the canvas is never tainted. */
const GEO_SVG_NS = 'http://www.w3.org/2000/svg';
const GEO_FONT = 'system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif';

function geoSvgEl(tag, attrs, text){
  const el = document.createElementNS(GEO_SVG_NS, tag);
  Object.entries(attrs || {}).forEach(([k, v]) => el.setAttribute(k, v));
  if (text != null) el.textContent = text;
  return el;
}

function geoLegendSVG(x, y, w){
  const m = geoMode();
  const g = geoSvgEl('g', {transform: `translate(${x},${y})`});
  const gid = 'geoLegGrad' + Date.now();
  const defs = geoSvgEl('defs');
  const lg = geoSvgEl('linearGradient', {id: gid, x1: '0%', y1: '0%', x2: '100%', y2: '0%'});
  m.stops.forEach(s => lg.appendChild(geoSvgEl('stop', {offset: (s[0]*100) + '%', 'stop-color': s[1]})));
  defs.appendChild(lg);
  g.appendChild(defs);

  g.appendChild(geoSvgEl('text', {x: 0, y: 0, 'font-family': GEO_FONT, 'font-size': 11,
    'font-weight': 600, fill: '#4b5563'}, m.legendTitle));
  g.appendChild(geoSvgEl('rect', {x: 0, y: 8, width: w, height: 12, rx: 2,
    fill: `url(#${gid})`, stroke: '#d1d5db', 'stroke-width': .8}));
  const labs = m.ends;
  const anchors = ['start', 'middle', 'end'];
  [0, w/2, w].forEach((lx, i) => {
    g.appendChild(geoSvgEl('text', {x: lx, y: 33, 'font-family': GEO_FONT, 'font-size': 10.5,
      'text-anchor': anchors[i], fill: '#6b7280'}, labs[i]));
  });
  // no-data key
  g.appendChild(geoSvgEl('rect', {x: 0, y: 42, width: 11, height: 11, rx: 2,
    fill: GEO_NODATA, stroke: '#d1d5db', 'stroke-width': .8}));
  g.appendChild(geoSvgEl('text', {x: 17, y: 51, 'font-family': GEO_FONT, 'font-size': 10.5,
    fill: '#6b7280'}, 'No genotype calls'));
  return g;
}

function geoBuildExportSVG(){
  const live = document.querySelector('#geoMap svg');
  if (!live) return null;
  const mapW = parseFloat(live.getAttribute('width')) || live.clientWidth || 600;
  const mapH = parseFloat(live.getAttribute('height')) || live.clientHeight || 400;
  const padT = 52, padB = 78, padX = 16;
  const W = Math.round(mapW), H = Math.round(mapH + padT + padB);

  const svg = geoSvgEl('svg', {
    xmlns: GEO_SVG_NS, 'xmlns:xlink': 'http://www.w3.org/1999/xlink',
    width: W, height: H, viewBox: `0 0 ${W} ${H}`
  });
  svg.appendChild(geoSvgEl('rect', {x: 0, y: 0, width: W, height: H, fill: '#ffffff'}));

  const row = GEO.rows[GEO.snpIndex];
  svg.appendChild(geoSvgEl('text', {x: padX, y: 24, 'font-family': GEO_FONT, 'font-size': 14,
    'font-weight': 700, fill: '#111827'}, `SNPGeo — ${geoScopeLabel()} · ${geoVariantLabel(row)}`));
  const sub = `${geoMode().label} · ${(GEO.input && GEO.input.datasetName) || GEO.dataset || ''}`;
  svg.appendChild(geoSvgEl('text', {x: padX, y: 41, 'font-family': GEO_FONT, 'font-size': 11,
    fill: '#6b7280'}, sub));

  const g = geoSvgEl('g', {transform: `translate(0,${padT})`});
  Array.from(live.childNodes).forEach(n => g.appendChild(n.cloneNode(true)));
  svg.appendChild(g);

  svg.appendChild(geoLegendSVG(padX, padT + mapH + 24, Math.min(300, W - 2*padX)));
  return {svg, W, H};
}

function geoExportName(ext){
  const parts = ['SNPGeo', geoScopeLabel(), geoVariantLabel(GEO.rows[GEO.snpIndex]), GEO.colorMode];
  return parts.join('_').replace(/[^A-Za-z0-9._-]+/g, '_').replace(/_+/g, '_') + '.' + ext;
}
function geoDownloadBlob(blob, name){
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = name;
  document.body.appendChild(a);
  a.click();
  setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 0);
}

function geoExportSVG(){
  const built = geoBuildExportSVG();
  if (!built){ console.warn('[SNPGeo] No map to export'); return; }
  const str = '<?xml version="1.0" encoding="UTF-8"?>\n'
            + new XMLSerializer().serializeToString(built.svg);
  geoDownloadBlob(new Blob([str], {type: 'image/svg+xml;charset=utf-8'}), geoExportName('svg'));
}

function geoExportPNG(scale){
  const built = geoBuildExportSVG();
  if (!built){ console.warn('[SNPGeo] No map to export'); return; }
  const k = scale || 2;
  const str = new XMLSerializer().serializeToString(built.svg);
  const url = URL.createObjectURL(new Blob([str], {type: 'image/svg+xml;charset=utf-8'}));
  const img = new Image();
  img.onload = () => {
    const cv = document.createElement('canvas');
    cv.width  = Math.round(built.W * k);
    cv.height = Math.round(built.H * k);
    const ctx = cv.getContext('2d');
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, cv.width, cv.height);
    ctx.drawImage(img, 0, 0, cv.width, cv.height);
    URL.revokeObjectURL(url);
    cv.toBlob(b => { if (b) geoDownloadBlob(b, geoExportName('png')); }, 'image/png');
  };
  img.onerror = () => {
    URL.revokeObjectURL(url);
    console.error('[SNPGeo] PNG export failed — falling back to SVG');
    geoExportSVG();
  };
  img.src = url;
}

/* ================= RIGHT PANE: OVERVIEW =================
   Default state of the right-hand pane: one row per country with the number
   of isolates carrying the SNP over the country's total isolate count, and
   that as a percentage. The country name is a button into the per-country
   view (geoShowDetail), which is also what a click on the map opens. */
function geoRenderOverview(){
  const det = document.getElementById('geoDetail');
  if (!det) return;
  GEO.selectedCountry = null;

  const stats = aggregateGeoData(GEO.snpIndex);
  const list = Object.values(stats);
  if (!list.length){
    det.innerHTML = `<div class="geo-pane"><div class="geo-pane-body">
      <div class="geo-empty">No geographic metadata for this dataset.</div></div></div>`;
    return;
  }

  // Carriers first (highest % at the top), then by collection size, then name.
  list.sort((a, b) => (b.pct - a.pct) || (b.total - a.total) || a.country.localeCompare(b.country));

  const totIso = list.reduce((s, r) => s + r.total, 0);
  const totCar = list.reduce((s, r) => s + r.count, 0);
  const nWith  = list.filter(r => r.count > 0).length;

  const rows = list.map(r => `
    <tr class="${r.count ? '' : 'geo-ov-zero'}">
      <td><button type="button" class="geo-ov-link" onclick="geoShowDetail('${escGeoAttr(r.iso3)}')"
            data-tt="Open the ${escGeoAttr(r.country)} breakdown">${escGeo(r.country)}</button></td>
      <td class="geo-cell-num">${r.count}/${r.total}</td>
      <td class="geo-cell-num">${(r.pct * 100).toFixed(1)}%</td>
    </tr>`).join('');

  det.innerHTML = `
    <div class="geo-pane">
      <div class="geo-pane-head">
        <h4 style="margin:0">Overview</h4>
        <span class="geo-pane-sub">${escGeo(geoVariantLabel(GEO.rows[GEO.snpIndex]))}</span>
      </div>
      <div class="geo-pane-body">
        <div class="geo-stat-block">
          <div class="geo-stat-row">
            <span class="geo-stat-label">Isolates with variant:</span>
            <span class="geo-stat-value">${totCar}/${totIso}</span>
          </div>
          <div class="geo-stat-row">
            <span class="geo-stat-label">Countries represented:</span>
            <span class="geo-stat-value">${nWith}/${list.length}</span>
          </div>
        </div>
        <table class="geo-region-table geo-ov-table">
          <thead>
            <tr>
              <th data-tt="Click a country to open its regional breakdown.">Country</th>
              <th class="geo-cell-num" data-tt="Isolates carrying the variant over all isolates known from that country.">Isolates</th>
              <th class="geo-cell-num" data-tt="Carriers as a percentage of the country's isolates.">%</th>
            </tr>
          </thead>
          <tbody>${rows}</tbody>
        </table>
      </div>
    </div>`;
  attachTT();
}

function geoShowDetail(iso3){
  const det = document.getElementById('geoDetail');
  if (!det) return;
  
  const stats = aggregateGeoData(GEO.snpIndex);
  const stat = stats[iso3];
  if (!stat) return;
  GEO.selectedCountry = iso3;
  
  const stateRows = Object.entries(stat.states)
    .sort((a, b) => b[1].count - a[1].count)
    .map(([state, s]) => `
      <tr>
        <td>${escGeo(state)}</td>
        <td class="geo-cell-num">${s.count}</td>
        <td class="geo-cell-num">${s.total}</td>
        <td class="geo-cell-num">${(s.total > 0 ? (s.count/s.total*100).toFixed(1) : 0)}%</td>
      </tr>`).join('');
  
  const isoPills = stat.isolates.length
    ? stat.isolates.map(id => `<span class="geo-iso-pill">${escGeo(id)}</span>`).join('')
    : '<span style="color:var(--faint);font-size:12px">None</span>';

  det.innerHTML = `
    <div class="geo-pane">
      <div class="geo-pane-head">
        <button class="geo-back-btn" onclick="geoCloseDetail()" data-tt="Back to the country overview">←</button>
        <h4 style="margin:0;flex:1">${escGeo(stat.country)}</h4>
      </div>
      <div class="geo-pane-body">
        <div class="geo-stat-block">
          <div class="geo-stat-row">
            <span class="geo-stat-label">Isolates with variant:</span>
            <span class="geo-stat-value">${stat.count}/${stat.total}</span>
          </div>
          <div class="geo-stat-row">
            <span class="geo-stat-label">Percentage:</span>
            <span class="geo-stat-value">${(stat.pct * 100).toFixed(1)}%</span>
          </div>
          <div class="geo-stat-row">
            <span class="geo-stat-label">Genotyped here:</span>
            <span class="geo-stat-value">${stat.called}${stat.missing ? ` (+${stat.missing} no call)` : ''}</span>
          </div>
          <div class="geo-stat-row">
            <span class="geo-stat-label">Ref / alt (called):</span>
            <span class="geo-stat-value">${stat.called ? `${stat.ref} / ${stat.alt} · ${(stat.altFrac*100).toFixed(1)}% alt` : '—'}</span>
          </div>
        </div>

        <h5 class="geo-subhead">Isolates carrying this variant</h5>
        <div class="geo-iso-pills">${isoPills}</div>

        <h5 class="geo-subhead">By region</h5>
        <table class="geo-region-table">
          <thead>
            <tr><th>Region</th><th>Count</th><th>Total</th><th>%</th></tr>
          </thead>
          <tbody>${stateRows}</tbody>
        </table>
      </div>
    </div>`;
  attachTT();
}

/* The pane is never empty now — leaving a country returns to the overview,
   which is the pane's resting state. */
function geoCloseDetail(){
  GEO.selectedCountry = null;
  geoRenderOverview();
}

/* ================= GEOJSON LOADING ================= */

/* Some countries in countries.geo.json carry non-contiguous overseas parts:
   France's MultiPolygon includes French Guiana, so shading FRA painted a patch
   of South America. Isolate counts are attributed per country, not per
   department, so the overseas rings carry no separate value — they just put the
   country's colour somewhere the reader doesn't expect. Restrict those
   countries to the polygons falling inside their mainland window.

   Keyed by iso3 -> [minLon, maxLon, minLat, maxLat]. Only add an entry where
   the overseas part is genuinely misleading; CAN/USA/RUS are deliberately
   absent because their outlying polygons (Alaska, Aleutians, Arctic islands)
   are part of how those countries are normally drawn. */
const GEO_MAINLAND = {
  FRA: [-10, 12, 40, 52],   // metropolitan France + Corsica; drops French Guiana
};

function geoTrimOverseas(geojson){
  (geojson.features || []).forEach(ft => {
    const box = GEO_MAINLAND[ft.properties && ft.properties.iso3];
    const gm = ft.geometry;
    if (!box || !gm || gm.type !== 'MultiPolygon') return;
    const [wLon, eLon, sLat, nLat] = box;
    const kept = gm.coordinates.filter(poly => {
      const ring = poly[0] || [];
      if (!ring.length) return false;
      // Keep the polygon if its centre sits inside the mainland window.
      let lon = 0, lat = 0;
      ring.forEach(pt => { lon += pt[0]; lat += pt[1]; });
      lon /= ring.length; lat /= ring.length;
      return lon >= wLon && lon <= eLon && lat >= sLat && lat <= nLat;
    });
    // Never blank a country out: if the window matched nothing, leave it alone.
    if (kept.length && kept.length < gm.coordinates.length) gm.coordinates = kept;
  });
  return geojson;
}

function loadCountriesGeoJSON(callback){
  if (GEO.geoJson){
    callback();
    return;
  }
  
  fetch('./data/geo/countries.geo.json')
    .then(r => r.json())
    .then(geojson => {
      // Trim before anything reads GEO.geoJson — the Mercator .fitSize() in
      // geoRenderMap() derives world bounds from it.
      GEO.geoJson = geoTrimOverseas(geojson);
      callback();
    })
    .catch(err => {
      console.error('[SNPGeo] Failed to load geojson:', err);
      const mapDiv = document.getElementById('geoMap');
      if (mapDiv) mapDiv.innerHTML = '<div class="geo-error">Failed to load geographic data</div>';
    });
}

/* ================= STYLES ================= */
function injectGeoCSS(){
  if (document.getElementById('snpgeo-css')) return;
  const s = document.createElement('style');
  s.id = 'snpgeo-css';
  s.textContent = `
  /* dataset chooser + search bar — mirrors SNPFold's .fold-ds-* / search bar exactly */
  .geo-ds-card{padding:12px 14px}
  .geo-ds-head{font-size:10.5px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.4px;margin-bottom:9px}
  .geo-ds-grid{display:flex;flex-wrap:wrap;gap:10px}
  .geo-ds{display:flex;align-items:center;gap:9px;text-align:left;cursor:pointer;background:#fff;
    border:1px solid var(--line);border-radius:10px;padding:9px 13px;min-width:210px;transition:border-color .12s,box-shadow .12s}
  .geo-ds:hover{border-color:#c3cee0}
  .geo-ds.sel{border-color:#9db4dd;box-shadow:0 0 0 2px rgba(47,106,208,.12);background:#f7faff}
  .geo-ds-dot{flex:0 0 auto;width:10px;height:10px;border-radius:50%;border:2px solid #c3cee0;background:#fff}
  .geo-ds.sel .geo-ds-dot{border-color:#2f6ad0;background:#2f6ad0;box-shadow:inset 0 0 0 2px #fff}
  .geo-ds-txt{display:flex;flex-direction:column;line-height:1.25;min-width:0}
  .geo-ds-name{font-weight:600;font-size:13px;color:var(--ink)}
  .geo-ds-sub{font-size:11.5px;color:var(--muted)}

  /* Map + detail card sit side by side at a fixed height (--geo-h): neither
     grows with the other's content. Grid's default align-items:stretch would
     otherwise size both to match whichever is tallest — e.g. a country with
     many regions/isolates would stretch the map's box (and leave dead white
     space below the actual leaflet canvas, itself hard-pinned at --geo-h).
     align-items:start plus an explicit height on both items pins the map's
     window size and pushes overflow in the detail card to its own internal
     scroll (.geo-pane-body already has overflow:auto; it just never had a
     bounded parent height to scroll within before). */
  .geo-shell { display: grid; grid-template-columns: 1fr 280px; gap: 14px; margin: 14px 0; align-items: start; --geo-h: 400px; }
  @media(max-width: 1000px) { .geo-shell { grid-template-columns: 1fr; } }
  
  .geo-map-wrap { position: relative; border: 1px solid var(--line); border-radius: 10px; background: #fff; overflow: hidden; height: var(--geo-h); }
  .geo-map { width: 100%; height: 100%; }
  .geo-loading, .geo-error { padding: 40px 20px; text-align: center; color: var(--muted); font-size: 14px; }
  
  /* Map toolbar: colour-mode picker + export buttons, over the map's top-right. */
  .geo-map-tools { position: absolute; top: 10px; right: 10px; z-index: 3; display: flex; gap: 6px; align-items: center; }
  .geo-mode-sel { font: inherit; font-size: 12px; padding: 5px 8px; border: 1px solid var(--line); border-radius: 7px;
    background: rgba(255,255,255,.96); color: var(--ink); cursor: pointer; max-width: 230px; }
  .geo-mode-sel:hover { border-color: #c3cee0; }
  .geo-tool-btn { font: inherit; font-size: 11px; font-weight: 600; letter-spacing: .02em; padding: 5px 9px;
    border: 1px solid var(--line); border-radius: 7px; background: rgba(255,255,255,.96); color: var(--muted); cursor: pointer; }
  .geo-tool-btn:hover { border-color: #c3cee0; color: var(--ink); }

  .geo-legend { position: absolute; bottom: 12px; left: 12px; background: rgba(255,255,255,0.95); border: 1px solid var(--line); border-radius: 8px; padding: 8px 12px; font-size: 12px; }
  .geo-leg-title { font-weight: 600; margin-bottom: 6px; color: var(--ink); }
  /* Continuous bar: the scale is continuous, and the reference/alternative
     view is only readable as a football field if the middle is visible. */
  .geo-leg-bar { width: 190px; height: 14px; border: 1px solid #d1d5db; border-radius: 3px; margin-bottom: 4px; }
  .geo-leg-labels { display: flex; justify-content: space-between; font-size: 11px; color: var(--muted); }
  .geo-leg-nd { display: flex; align-items: center; gap: 5px; margin-top: 5px; font-size: 11px; color: var(--muted); }
  .geo-leg-nd-sw { width: 11px; height: 11px; border-radius: 3px; background: #e4e8ee; border: 1px solid #d1d5db; }
  
  .geo-detail { border: 1px solid var(--line); border-radius: 10px; background: #fff; position: relative; height: var(--geo-h); overflow: hidden; }
  .geo-pane { height: 100%; display: flex; flex-direction: column; overflow: hidden; }
  .geo-pane-head { padding: 10px 12px; border-bottom: 1px solid var(--line); display: flex; align-items: center; justify-content: space-between; gap: 8px; }
  .geo-pane-sub { font-family: var(--mono); font-size: 11px; color: var(--muted); white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
  .geo-close-btn { background: none; border: none; font-size: 20px; cursor: pointer; color: var(--muted); }
  .geo-back-btn { background: none; border: 1px solid var(--line); border-radius: 6px; width: 24px; height: 24px; line-height: 1;
    font-size: 14px; cursor: pointer; color: var(--muted); flex: 0 0 auto; }
  .geo-back-btn:hover { border-color: #c3cee0; color: var(--ink); }

  /* Overview table — country name is the affordance into the country view. */
  .geo-ov-table td { vertical-align: middle; }
  .geo-ov-link { background: none; border: none; padding: 0; font: inherit; font-size: 12px; color: var(--blue-600);
    cursor: pointer; text-align: left; }
  .geo-ov-link:hover { text-decoration: underline; }
  .geo-ov-zero td { color: var(--faint); }
  .geo-ov-zero .geo-ov-link { color: var(--muted); }
  .geo-pane-body { flex: 1; overflow: auto; padding: 12px; }
  
  .geo-stat-block { background: #f8fafc; border-radius: 8px; padding: 10px; margin-bottom: 12px; }
  .geo-stat-row { display: flex; justify-content: space-between; padding: 4px 0; font-size: 13px; }
  .geo-stat-label { color: var(--muted); }
  .geo-stat-value { font-weight: 600; color: var(--blue-600); font-family: var(--mono); }

  .geo-subhead { margin: 12px 0 6px; font-size: 12px; color: var(--muted); text-transform: uppercase; letter-spacing: .03em; }
  .geo-iso-pills { display: flex; flex-wrap: wrap; gap: 6px; margin-bottom: 4px; }
  .geo-iso-pill { display: inline-block; padding: 3px 10px; border-radius: 20px; font-size: 11px; font-weight: 600; font-family: var(--mono); background: #e6eeff; color: #3358bf; white-space: nowrap; }

  .geo-region-table { width: 100%; border-collapse: collapse; font-size: 12px; }
  .geo-region-table th { text-align: left; padding: 6px 4px; border-bottom: 1px solid var(--line); color: var(--muted); font-size: 11px; text-transform: uppercase; letter-spacing: .03em; }
  .geo-region-table td { padding: 5px 4px; border-bottom: 1px solid var(--line-2); }
  .geo-cell-num { text-align: right; font-family: var(--mono); font-size: 11px; }
  
  /* SNP table: reuses SNPVersity's global .tbl-wrap / table.vcf theme (navy
     sticky header, .c-mono, .pill, .score gradient cells) so it looks and
     feels identical. These few rules only add what table.vcf doesn't already
     define: the row-select indicator column and a persistent "active" tint
     for whichever SNP's map is currently shown (table.vcf's own :hover is
     transient and isn't enough on its own to mark the selected row). */
  .geo-table { margin-top: 14px; }
  .geo-hint { display: flex; align-items: center; gap: 4px; margin: 0 2px 7px;
    font-size: 11.5px; color: var(--muted); }
  .geo-hint kbd { display: inline-flex; align-items: center; justify-content: center;
    min-width: 17px; height: 17px; padding: 0 3px; border: 1px solid var(--line);
    border-bottom-width: 2px; border-radius: 4px; background: #fff;
    font-family: var(--mono); font-size: 11px; line-height: 1; color: var(--navy-600); }
  /* Stop the page (and so the map) from scrolling on once the variant list
     hits its own top or bottom. Height itself is set by geoFitTableHeight(). */
  .geo-table .tbl-wrap { overscroll-behavior: contain; }
  table.vcf.geo-vcf tbody tr { cursor: pointer; }
  table.vcf.geo-vcf tbody tr.active td { background: #eaf2ff; }
  table.vcf.geo-vcf .geo-cb { width: 28px; text-align: center; }
  /* A radio, not a checkbox: the map shows exactly one variant at a time, and
     a square box wrongly implied several could be ticked at once. tabindex=-1
     keeps Tab from walking through every row — up/down arrows do that. */
  table.vcf.geo-vcf .geo-radio { width: 14px; height: 14px; margin: 0; vertical-align: middle;
    accent-color: var(--blue-600); cursor: pointer; }

  .geo-empty { padding: 26px; text-align: center; color: var(--faint); }
  `;
  document.head.appendChild(s);
}
