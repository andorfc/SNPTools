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
  snpStats: {},          // snpIndex -> { countryIso3: {count, total, pct, isolates, states} }
  gene: null,            // gene currently in the search box (example until chosen)
  geneUser: false,       // true once the user types/loads a specific gene
  geoJson: null,         // loaded countries.geo.json
  mapSvg: null,          // d3 selection of map container
};

/* ================= HELPERS ================= */
function escGeo(s){
  return String(s==null?'':s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
}
function escGeoAttr(s){
  return escGeo(s).replace(/"/g,'&quot;').replace(/'/g,'&#39;');
}

/* Default example gene for the CURRENT reference (graminearum -> FGSG_,
   7600 -> FVEG_, MRC826 -> FVERT4_), via Data.exampleGenes() — the same source
   SNPFold uses. Mirrors SNPFold's behaviour: the search box is PREFILLED with a
   real, editable value instead of a greyed-out placeholder, so the page is one
   click away from a working example. */
function geoDefaultExampleGene(){
  const ex = (typeof Data !== 'undefined' && Data.exampleGenes) ? Data.exampleGenes(S.dataset) : null;
  return (ex && ex[0]) || 'FGSG_00025';
}
/* Keeps the prefill following the active reference. A gene the user actually
   chose (GEO.geneUser) is kept — EXCEPT when it belongs to a different reference
   than the one now selected: an FGSG_ id left in the box after switching to a
   verticillioides reference can only ever return "not found", so the prefill
   falls back to that reference's own example. Prefix (FGSG_ / FVEG_ / FVERT4_)
   is what identifies the reference here, same convention the ID schemes use. */
function geoIdPrefix(id){
  const m = /^([A-Za-z0-9]+_)/.exec(String(id || ''));
  return m ? m[1].toUpperCase() : '';
}
function geoEnsureDefaultGene(){
  const ex = geoDefaultExampleGene();
  if (GEO.geneUser && geoIdPrefix(GEO.gene) !== geoIdPrefix(ex)) GEO.geneUser = false;
  if (!GEO.geneUser) GEO.gene = ex;
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

/* ---- search bar: mirrors SNPFold's searchBar() exactly, minus the model-source picker ---- */
function geoSearchBar(){
  const gene = (S.geoInput && S.geoInput.gene) || GEO.gene || geoDefaultExampleGene();
  return `<div class="card pad" style="margin-bottom:16px;display:flex;gap:10px;align-items:center;flex-wrap:wrap">
    <span style="font-size:10.5px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.4px">Gene model</span>
    <input id="geoGeneInput" value="${escGeoAttr(gene)}" spellcheck="false"
      style="flex:1;min-width:280px;border:1px solid var(--line);border-radius:9px;padding:9px 11px;font-family:var(--mono);font-size:13px"
      onkeydown="if(event.key==='Enter')geoLookupGene()"
      placeholder="e.g. ${escGeoAttr(geoDefaultExampleGene())}">
    <button class="btn" onclick="geoLookupGene()">Search</button>
  </div>`;
}

function renderGeo(page){
  const hasResults = S.geoInput && S.geoInput.rows && S.geoInput.rows.length;
  page.className = 'page fade';
  injectGeoCSS();
  geoEnsureDefaultGene();   // prefill the search box with this reference's example gene

  const header = geoDatasetChooser() + geoSearchBar()
    + `<div id="geoLookupStatus" style="margin:-8px 0 16px;font-size:12px;color:var(--muted)"></div>`;

  const results = hasResults ? `
    <div class="card pad" style="margin-bottom:16px">
      <div class="n" style="margin:0">${S.geoInput.gene ? 'Gene: ' + escGeo(S.geoInput.gene) : 'Region: chr' + S.geoInput.chr + ':' + S.geoInput.start.toLocaleString() + '-' + S.geoInput.end.toLocaleString()}</div>
      <div style="color:var(--muted);font-size:13px;margin-top:4px">${S.geoInput.rows.length} variants · click a row, then use <kbd class="geo-kbd">&uarr;</kbd><kbd class="geo-kbd">&darr;</kbd> to step through them</div>
    </div>
    <div class="geo-shell">
      <div class="geo-map-wrap">
        <div id="geoMap" class="geo-map"></div>
        <div id="geoLegend" class="geo-legend"></div>
      </div>
      <div id="geoDetail" class="geo-detail"></div>
    </div>
    <div class="geo-table" id="geoTableContainer">
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
      /* select(), not just focus(): the box is prefilled, so selecting the text
         means typing a different gene replaces it instead of appending to it. */
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
  
  GEO.gene = gene; GEO.geneUser = true;   // an explicit choice — stop tracking the example
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
  geoEnsureDefaultGene();   // follow the new reference's example gene (FGSG_/FVEG_/FVERT4_)
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
        count: 0,
        total: 0,
        pct: 0,
        isolates: [],
        states: {}
      };
    }
    
    // Count total isolates in this country (regardless of whether they were queried)
    result[iso3].total += 1;
    
    if (!result[iso3].states[state]){
      result[iso3].states[state] = {count: 0, total: 0, isolates: []};
    }
    result[iso3].states[state].total += 1;
  });
  
  // STEP 2: For queried isolates that carry the variant, mark them
  GEO.accs.forEach((acc, accIdx) => {
    const accId = String(acc.id);
    const geo = REGIONS[accId];
    if (!geo) return; // unmapped isolate
    
    const iso3 = geo.iso3 || '???';
    const state = geo.state || geo.county || 'Unknown';
    
    // Check if this accession carries the variant
    // Variants may be in row.gts (array) or row.gt_N columns
    let hasVar = false;
    if (row.gts && Array.isArray(row.gts)) {
      const gt = row.gts[accIdx];
      hasVar = gt && gt !== '0' && gt !== '0|0' && gt !== '0/0';
    } else {
      const gtCol = `gt_${accIdx}`;
      const gt = row[gtCol];
      hasVar = gt && gt !== '0' && gt !== '0|0' && gt !== '0/0';
    }
    
    if (hasVar){
      result[iso3].count += 1;
      result[iso3].isolates.push(accId);
      result[iso3].states[state].count += 1;
      result[iso3].states[state].isolates.push(accId);
    }
  });
  
  // Calculate percentages
  Object.values(result).forEach(r => {
    r.pct = r.total > 0 ? (r.count / r.total) : 0;
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
  
  // Score columns mirror SNPVersity exactly: same model list/labels/tooltips
  // (Data.scoreModels), same DNA-vs-protein tinting, same gColor() gradient.
  const models = Data.scoreModels(GEO.dataset);
  const sc = (v) => `<td class="score ${v==null?'na':''}" style="${v==null?'':'background:'+gColor(v)}">${v==null?'N/A':v}</td>`;

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
        <td class="geo-cb"><span class="rb" role="radio" tabindex="${active?'0':'-1'}" aria-checked="${active?'true':'false'}" aria-label="Map ${escGeoAttr(variant)}"></span></td>
        <td class="c-mono">${escGeo(variant)}</td>
        <td class="effect-cell">${escGeo(consequence)}${impact?` <span class="pill ${impact}">${escGeo(row.impact)}</span>`:''}</td>
        <td class="num">${resi != null ? resi : '—'}</td>
        <td class="c-mono">${resiVariant ? escGeo(resiVariant) : '—'}</td>
        ${models.map(m=>sc(row[m.key])).join('')}
      </tr>`;
  }).join('');
  
  tbl.innerHTML = `<div class="tbl-wrap"><table class="vcf geo-vcf">
    <thead>${thead}</thead>
    <tbody role="radiogroup" aria-label="Select a variant to map" onkeydown="geoTableKey(event)">${tbody}</tbody>
  </table></div>`;
  attachTT();
}

/* The map SVG is sized from its container at render time, so a window resize
   (or a --geo-h breakpoint crossing) needs a re-render to re-fit the
   projection. Debounced, and a no-op unless the map is actually on screen. */
let _geoResizeT = null;
window.addEventListener('resize', () => {
  if (!document.getElementById('geoMap')) return;
  clearTimeout(_geoResizeT);
  _geoResizeT = setTimeout(() => { if (document.getElementById('geoMap')) geoRenderMap(); }, 180);
});

/* Move the "active" marker between two rows WITHOUT rebuilding the tbody.
   geoRenderTable() re-writes innerHTML, which would destroy focus on every
   keystroke (and re-render hundreds of rows per arrow press). Swapping the
   class/aria/tabindex on just the outgoing and incoming rows keeps focus alive
   and makes held-down arrow keys smooth. Returns the new row, or null if the
   table isn't on the page. */
function geoSetActiveRow(idx){
  const tb = document.querySelector('table.geo-vcf tbody');
  if (!tb || !tb.rows.length) return null;
  const mark = (tr, on) => {
    if (!tr) return;
    tr.classList.toggle('active', on);
    const rb = tr.querySelector('.rb');
    if (rb){ rb.setAttribute('aria-checked', on ? 'true' : 'false'); rb.tabIndex = on ? 0 : -1; }
  };
  mark(tb.rows[GEO.snpIndex], false);
  const next = tb.rows[idx];
  mark(next, true);
  return next;
}

/* opts.viaKeyboard — scroll the row into view inside .tbl-wrap and keep focus on
   its radio so the next arrow press still lands here. */
function geoSelectSnp(idx, opts){
  opts = opts || {};
  if (!GEO.rows.length) return;
  idx = Math.max(0, Math.min(GEO.rows.length - 1, idx));
  const row = geoSetActiveRow(idx);
  GEO.snpIndex = idx;
  GEO.selectedCountry = null;
  if (!row) geoRenderTable();          // table not built yet — fall back to a full render
  else {
    /* focus the radio so arrow keys work straight after a mouse click, and keep
       it inside the scroller as the selection walks off-screen. preventScroll
       stops the browser also yanking the page. */
    const rb = row.querySelector('.rb');
    if (rb) rb.focus({preventScroll: true});
    if (opts.viaKeyboard) row.scrollIntoView({block: 'nearest'});
  }
  geoRenderMap();
  geoCloseDetail();
}

/* Standard ARIA radiogroup keyboard model: arrows step the selection, Home/End
   jump to the ends, PageUp/PageDown move by ten. Every handled key calls
   preventDefault so the page (and the .tbl-wrap scroller) doesn't scroll at the
   same time the selection moves. */
function geoTableKey(e){
  const n = GEO.rows.length;
  if (!n) return;
  const i = GEO.snpIndex;
  let next;
  switch (e.key){
    case 'ArrowDown': case 'ArrowRight': next = i + 1; break;
    case 'ArrowUp':   case 'ArrowLeft':  next = i - 1; break;
    case 'PageDown':  next = i + 10; break;
    case 'PageUp':    next = i - 10; break;
    case 'Home':      next = 0; break;
    case 'End':       next = n - 1; break;
    default: return;
  }
  e.preventDefault();
  next = Math.max(0, Math.min(n - 1, next));
  if (next !== i) geoSelectSnp(next, {viaKeyboard: true});
}

/* ================= MAP RENDERING (PLACEHOLDER) ================= */
function geoRenderMap(){
  const mapDiv = document.getElementById('geoMap');
  if (!mapDiv) return;
  
  const stats = aggregateGeoData(GEO.snpIndex);
  const row = GEO.rows[GEO.snpIndex];
  const variant = row.variant || 'variant ' + GEO.snpIndex;
  
  if (!GEO.geoJson){
    mapDiv.innerHTML = '<div class="geo-loading">Loading map...</div>';
    return;
  }
  
  mapDiv.innerHTML = '';
  const width = mapDiv.clientWidth || 600;
  /* Height comes from the container (.geo-map is 100% of .geo-map-wrap, which
     is height:var(--geo-h)) so the responsive --geo-h steps below actually
     resize the projection instead of cropping a fixed 400px SVG. */
  const height = mapDiv.clientHeight || 400;
  
  // Create SVG
  const svg = d3.select(mapDiv).append('svg')
    .attr('width', width)
    .attr('height', height)
    .attr('viewBox', `0 0 ${width} ${height}`);
  
  // Simple Mercator projection
  const projection = d3.geoMercator()
    .fitSize([width, height], GEO.geoJson);
  
  const pathGenerator = d3.geoPath().projection(projection);
  
  // Color scale: white (0%) to dark green (100%)
  const colorScale = d3.scaleLinear()
    .domain([0, 1])
    .range(['#f0fdf4', '#059669']);
  
  // Draw countries
  svg.selectAll('path.country')
    .data(GEO.geoJson.features)
    .enter()
    .append('path')
    .attr('class', 'country')
    .attr('d', pathGenerator)
    .style('fill', d => {
      const iso3 = d.properties.iso3;
      const stat = stats[iso3];
      const pct = stat ? stat.pct : 0;
      return colorScale(pct);
    })
    .style('stroke', '#ddd')
    .style('stroke-width', '0.5px')
    .style('cursor', 'pointer')
    .on('mouseenter', function(d){
      const iso3 = d.properties.iso3;
      const stat = stats[iso3];
      if (!stat) return;
      d3.select(this).style('stroke', '#000').style('stroke-width', '1.5px');
      const pct = (stat.pct * 100).toFixed(1);
      const tt = document.getElementById('tt');
      if (tt){
        tt.innerHTML = `<b>${escGeo(stat.country)}</b><br/>${stat.count}/${stat.total} (${pct}%)`;
        tt.classList.add('show');
      }
    })
    .on('mousemove', function(e){
      const tt = document.getElementById('tt');
      if (tt){
        tt.style.left = (e.clientX + 12) + 'px';
        tt.style.top = (e.clientY + 14) + 'px';
      }
    })
    .on('mouseleave', function(){
      d3.select(this).style('stroke', '#ddd').style('stroke-width', '0.5px');
      const tt = document.getElementById('tt');
      if (tt) tt.classList.remove('show');
    })
    .on('click', function(d){
      const iso3 = d.properties.iso3;
      GEO.selectedCountry = iso3;
      geoShowDetail(iso3);
    });
  
  // Legend
  geoRenderLegend();
}

function geoRenderLegend(){
  const leg = document.getElementById('geoLegend');
  if (!leg) return;
  
  const stats = aggregateGeoData(GEO.snpIndex);
  const minPct = stats && Object.values(stats).length > 0 
    ? Math.min(...Object.values(stats).map(s => s.pct * 100))
    : 0;
  const maxPct = stats && Object.values(stats).length > 0
    ? Math.max(...Object.values(stats).map(s => s.pct * 100))
    : 100;
  
  const colorScale = d3.scaleLinear()
    .domain([0, 1])
    .range(['#f0fdf4', '#059669']);
  
  leg.innerHTML = `
    <div class="geo-leg-title">% carrying variant</div>
    <div class="geo-leg-bar">
      ${[0, 0.25, 0.5, 0.75, 1].map(v => `
        <div class="geo-leg-swatch" style="background:${colorScale(v)}"></div>
      `).join('')}
    </div>
    <div class="geo-leg-labels">
      <span>0%</span>
      <span>50%</span>
      <span>100%</span>
    </div>`;
}

function geoShowDetail(iso3){
  const det = document.getElementById('geoDetail');
  if (!det) return;
  
  const stats = aggregateGeoData(GEO.snpIndex);
  const stat = stats[iso3];
  if (!stat) return;
  
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
        <h4 style="margin:0">${escGeo(stat.country)}</h4>
        <button class="geo-close-btn" onclick="geoCloseDetail()">×</button>
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
}

function geoCloseDetail(){
  const det = document.getElementById('geoDetail');
  if (det) det.innerHTML = '';
  GEO.selectedCountry = null;
}

/* ================= GEOJSON LOADING ================= */
function loadCountriesGeoJSON(callback){
  if (GEO.geoJson){
    callback();
    return;
  }
  
  fetch('./data/geo/countries.geo.json')
    .then(r => r.json())
    .then(geojson => {
      GEO.geoJson = geojson;
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
  /* The map is what this page is for, so it is pinned directly beneath the
     sticky .topbar (60px) and stays fully visible while the variant table is
     scrolled. The opaque --paper background is what stops table rows showing
     through the pinned band. --geo-h / --geo-top live on :root so the table
     rule below can size itself against them. */
  /* --geo-gap = everything stacked below the pinned map that still has to fit
     on screen: .geo-shell's 14px top+bottom padding, .geo-table's 14px margin,
     .page's 60px bottom padding. Without it the page out-scrolls the pinned
     block and the table's own header slides up underneath the map. */
  :root { --geo-h: 400px; --geo-top: 60px; --geo-gap: 116px; }
  /* On shorter windows a 400px map would leave only a few table rows under it,
     so step the map down rather than squeezing the table to nothing. */
  @media(max-height: 900px) { :root { --geo-h: 340px; } }
  @media(max-height: 820px) { :root { --geo-h: 300px; } }
  .geo-shell { display: grid; grid-template-columns: 1fr 280px; gap: 14px; align-items: start;
    margin: 0 0 14px; padding: 14px 0;
    position: sticky; top: var(--geo-top); z-index: 20; background: var(--paper); }
  /* Bound the table to whatever viewport is left under the pinned map, so
     table.vcf's navy sticky thead parks right below the map rather than
     partway down the page. */
  .geo-table .tbl-wrap { max-height: calc(100vh - var(--geo-top) - var(--geo-h) - var(--geo-gap)); min-height: 180px; }
  /* Stacked (map above detail pane) or a short window: the pinned block would
     eat the whole viewport, so fall back to normal flow. */
  @media(max-width: 1000px) { .geo-shell { grid-template-columns: 1fr; position: static; }
    .geo-table .tbl-wrap { max-height: 72vh; } }
  @media(max-height: 700px) { .geo-shell { position: static; }
    .geo-table .tbl-wrap { max-height: 72vh; } }
  
  .geo-map-wrap { position: relative; border: 1px solid var(--line); border-radius: 10px; background: #fff; overflow: hidden; height: var(--geo-h); }
  .geo-map { width: 100%; height: 100%; }
  .geo-loading, .geo-error { padding: 40px 20px; text-align: center; color: var(--muted); font-size: 14px; }
  
  .geo-legend { position: absolute; bottom: 12px; left: 12px; background: rgba(255,255,255,0.95); border: 1px solid var(--line); border-radius: 8px; padding: 8px 12px; font-size: 12px; }
  .geo-leg-title { font-weight: 600; margin-bottom: 6px; color: var(--ink); }
  .geo-leg-bar { display: flex; gap: 2px; height: 20px; margin-bottom: 4px; }
  .geo-leg-swatch { flex: 1; border: 1px solid #ddd; }
  .geo-leg-labels { display: flex; justify-content: space-between; font-size: 11px; color: var(--muted); }
  
  .geo-detail { border: 1px solid var(--line); border-radius: 10px; background: #fff; position: relative; height: var(--geo-h); overflow: hidden; }
  .geo-pane { height: 100%; display: flex; flex-direction: column; overflow: hidden; }
  .geo-pane-head { padding: 10px 12px; border-bottom: 1px solid var(--line); display: flex; align-items: center; justify-content: space-between; }
  .geo-close-btn { background: none; border: none; font-size: 20px; cursor: pointer; color: var(--muted); }
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
  table.vcf.geo-vcf tbody tr { cursor: pointer; }
  table.vcf.geo-vcf tbody tr.active td { background: #eaf2ff; }
  table.vcf.geo-vcf .geo-cb { width: 28px; text-align: center; }
  /* Exactly one variant is mapped at a time, so the row selector is a radio,
     not a checkbox: circular, with the classic filled-dot-in-a-ring checked
     state. The <tbody> carries role="radiogroup" and each cell role="radio". */
  table.vcf.geo-vcf .rb { display: inline-block; width: 14px; height: 14px; box-sizing: border-box;
    border: 1.5px solid #c3ccda; border-radius: 50%; background: #fff; vertical-align: middle; transition: border-color .12s, background .12s; }
  table.vcf.geo-vcf tbody tr:hover .rb { border-color: #9db4dd; }
  /* The selected radio holds focus so arrow keys keep working; show that. */
  table.vcf.geo-vcf .rb:focus-visible { outline: 2px solid var(--blue-600); outline-offset: 2px; }
  table.vcf.geo-vcf tr.active .rb { background: var(--blue-600); border-color: var(--blue-600); box-shadow: inset 0 0 0 2.5px #fff; }

  .geo-kbd { display: inline-block; min-width: 15px; text-align: center; padding: 0 4px; margin: 0 1px;
    border: 1px solid var(--line); border-bottom-width: 2px; border-radius: 4px; background: #fff;
    font-family: var(--mono); font-size: 11px; line-height: 16px; color: var(--muted); }

  .geo-empty { padding: 26px; text-align: center; color: var(--faint); }
  `;
  document.head.appendChild(s);
}
