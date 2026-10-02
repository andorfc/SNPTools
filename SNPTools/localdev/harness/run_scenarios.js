/* run_scenarios.js <siteRoot> <outDir> -- end-to-end checks of SNPTrait and SNPGeo
 * on the zmgrin2026_imp dataset; writes results.json, map PNG/SVG renders
 * (the app's own geoBuildExportSVG() output, rasterised with resvg) and static
 * HTML snapshots of the rendered pages. */
const fs = require('fs'), path = require('path');
const {Resvg} = require('@resvg/resvg-js');
const {openSite} = require('./ui_harness');
const [,, ROOT, OUT] = process.argv;
fs.mkdirSync(OUT, {recursive: true});
const R = {};

function snapshot(site, name, title){
  const w = site.window, d = w.document;
  const css = ['css/main.css', 'css/paneffect-native.css'].map(f => fs.readFileSync(path.join(ROOT, f), 'utf8')).join('\n');
  const inj = [...d.querySelectorAll('style')].map(s => s.textContent).join('\n');
  const body = d.body.cloneNode(true);
  body.querySelectorAll('script').forEach(s => s.remove());
  const html = `<!DOCTYPE html><html><head><meta charset="utf-8"><title>${title}</title>
<!-- Static DOM snapshot produced by harness/run_scenarios.js (jsdom); scripts removed. -->
<style>${css}\n${inj}</style></head>${body.outerHTML}</html>`;
  fs.writeFileSync(path.join(OUT, name + '.html'), html);
}
function mapPNG(site, name){
  const built = site.$eval('geoBuildExportSVG()');
  if (!built) return null;
  const str = new site.window.XMLSerializer().serializeToString(built.svg);
  fs.writeFileSync(path.join(OUT, name + '.svg'), str);
  const png = new Resvg(str, {fitTo: {mode: 'zoom', value: 2}, font: {loadSystemFonts: true}, background: '#ffffff'}).render().asPng();
  fs.writeFileSync(path.join(OUT, name + '.png'), png);
  return {w: built.W, h: built.H, bytes: png.length};
}
/* SNPVersity annotation columns as rendered: header labels + status, and per column the
   number of body cells with a value (every page of the current filter, via rowHTML). */
function annotStats(site){
  return site.$eval(`(()=>{
    const th=[...document.querySelectorAll('#rtBody table.vcf thead th')].map(t=>t.textContent.trim());
    const cols=[...document.querySelectorAll('#rtBody table.vcf thead th[data-col]')].map(t=>({key:t.dataset.col, label:t.textContent.trim(), status:t.dataset.status}));
    const first=th.indexOf('Gene model');
    const host=document.createElement('tbody');
    host.innerHTML=tableView(S.results.rows).fr.map(r=>rowHTML(r)).join('');
    const filled={};
    cols.forEach((c,j)=>{ let n=0; host.querySelectorAll('tr').forEach(tr=>{ const td=tr.children[first+j];
      const v=td?td.textContent.trim():''; if(v && !/^(—|N\\/A|NA|\\.|intergenic)$/.test(v)) n++; }); filled[c.key]=n; });
    return {dataset:S.dataset, headers:th.slice(0, first+cols.length), cols, rows:S.results.rows.length, filled,
            note:(document.getElementById('annotNote')||{}).textContent||''};
  })()`);
}
async function until(site, expr, ms = 8000){
  const t0 = Date.now();
  while (Date.now() - t0 < ms){ if (site.$eval(expr)) return true; await site.wait(50); }
  return false;
}

(async () => {
  const site = await openSite(ROOT);
  const $ = site.$eval;
  /* which zmgrin2026 stores are installed: the 3,495-site demo store or a full chromosome */
  try {
    R.store = JSON.parse(require('child_process').execFileSync(process.env.PYTHON_PATH || 'python3', ['-c',
      'import h5py,json,glob,os,sys\nr={}\nfor f in sorted(glob.glob(sys.argv[1]+"/hdf5/version3/zmgrin2026_chr*_impute.h5")):\n  h=h5py.File(f,"r"); r[os.path.basename(f).split("_")[1]]=int(h["POS"].shape[0]); h.close()\nprint(json.dumps(r))', ROOT]).toString());
  } catch (e) { R.store = {error: String(e.message || e).slice(0, 200)}; }
  R.load = {registry: $('Object.keys(SNPTools.registry)'), nav: $('[...document.querySelectorAll("#nav .navitem")].map(b=>b.textContent.trim().replace(/\\s+/g," "))'),
            datasets: $('Data.datasets().map(d=>d.id)'), zmgrinN: $('Data.accessionsFor("zmgrin2026_imp").length'),
            otherFamilies: $('Object.keys(window.SNP_CATALOG.families)'),
            mgdb2026N: $('Data.accessionsFor("mgdb2026_hq").length')};

  /* ---------- SNPTrait ---------- */
  $('S.dataset="zmgrin2026_imp"; go("snptrait")');
  const T = R.snptrait = {};
  T.header = $('document.querySelector("#page h2").textContent');
  T.schemaGroupBy = $('TRAIT.schema.groupBy'); T.hasMeta = $('TRAIT.hasMeta'); T.rows = $('TRAIT.rows.length');
  T.facetKeys = $('TRAIT.schema.facets.map(f=>f[0])');
  T.facetCounts = $('JSON.parse(JSON.stringify(facetCounts(traitVisible())))');
  T.domGroups = $('[...document.querySelectorAll("#traitGrid .tg-head b")].map(b=>b.textContent)');
  T.domFacetBlocks = $('document.querySelectorAll("#traitFacets .facet-block").length');
  T.statusInitial = $('document.getElementById("traitStatus").textContent');
  snapshot(site, 'snptrait_initial', 'SNPTrait - zmgrin2026_imp (initial)');
  $('traitToggleFacet("subpop","SS",true); traitToggleFacet("inAmes282","yes",true); traitToggleFacet("kernelType","Dent",true)');
  T.filterSS_Ames_Dent = $('traitVisible().map(r=>r.id).sort()');
  T.domRowsFiltered = $('document.querySelectorAll("#traitGrid tbody tr").length');
  T.statusFiltered = $('document.getElementById("traitStatus").textContent');
  $('TRAIT.q="iowa"; traitRenderTable(); traitRenderFacets(); traitStatus()');
  T.filterPlusIowa = $('traitVisible().map(r=>r.id).sort()');
  T.statusIowa = $('document.getElementById("traitStatus").textContent');
  // numeric range facet from data/traits/zmgrin2026.traits.json (lazy-loaded)
  $('TRAIT.q=""; traitClearFacets()');
  await until(site, 'traitSideFile() && traitSideFile()!=="loading"');
  T.rangeTraits = $('traitNumericTraits().length');
  $('traitAddRange("KERNEL-WEIGHT-1000", 250, 300)');
  T.rangeKW = $('traitVisible().map(r=>r.id).sort()');
  T.rangeKWHeader = $('[...document.querySelectorAll("#traitGrid thead th")].map(t=>t.textContent).filter(t=>/mean/.test(t))[0]||null');
  T.rangeKWStatus = $('document.getElementById("traitStatus").textContent');
  $('traitToggleFacet("panel","Ames282",true)');
  T.rangeKW_Ames = $('traitVisible().length');
  snapshot(site, 'snptrait_range', 'SNPTrait - 1000-kernel weight 250-300 g x Ames282');
  $('traitClearFacets(); traitToggleFacet("subpop","SS",true); traitToggleFacet("inAmes282","yes",true); traitToggleFacet("kernelType","Dent",true)');
  $('TRAIT.q=""; traitSelectVisible(true)');
  T.selectedAfterSelectVisible = $('TRAIT.selected.size');
  T.runbar = $('document.getElementById("traitRunbar").textContent.replace(/\\s+/g," ").trim()');
  snapshot(site, 'snptrait_filtered', 'SNPTrait - SS x Ames282 x Dent');
  $('traitSendToVersity()'); await site.wait(300);
  T.handoff = {tool: $('S.tool'), selected: $('S.selected.size'), dataset: $('S.dataset')};
  // a family without a generated schema falls back to the neutral schema
  $('S.dataset="mgdb2026_hq"; go("snptrait")');
  T.mgdb2026 = {groupBy: $('TRAIT.schema.groupBy'), facets: $('TRAIT.schema.facets.map(f=>f[0])'), rows: $('TRAIT.rows.length'), header: $('document.querySelector("#page h2").textContent')};

  /* ---------- Help page ---------- */
  $('openHelp()'); await site.wait(150);
  R.help = {hasSNPGeo: $('/SNPGeo/.test(document.getElementById("page").textContent)'),
            traitLive: $('(()=>{const t=document.getElementById("page").textContent; return /Select lines by passport and trait metadata/.test(t);})()')};

  /* ---------- SNPGeo: gene search ---------- */
  $('S.dataset="zmgrin2026_imp"; S.geoInput=null; go("snpgeo")');
  const G = R.snpgeo_gene = {};
  G.emptyState = $('document.getElementById("geoGeneInput").value');
  $('document.getElementById("geoGeneInput").value="Zm00001eb374090"');
  await $('geoLookupGene()');
  G.status = $('(document.getElementById("geoLookupStatus")||{}).textContent');
  await until(site, 'GEO.rows && GEO.rows.length && document.querySelector("#geoMap svg")');
  await site.wait(300);
  G.header = $('document.querySelector(".page .card .n").textContent');
  G.coverage = $('[...document.querySelectorAll(".page .card")].map(c=>c.textContent).find(t=>/genotyped/.test(t)).replace(/\\s+/g," ").trim()');
  G.nRows = $('GEO.rows.length'); G.nAccs = $('GEO.accs.length'); G.gtsType = $('Object.prototype.toString.call(GEO.rows[0].gts)');
  G.view = $('geoView()'); G.admin1Paths = $('document.querySelectorAll("#geoMap path.admin1").length');
  G.tableRows = $('document.querySelectorAll("#geoTable tbody tr").length');
  G.tableHead = $('[...document.querySelectorAll("#geoTable thead th")].map(t=>t.textContent.trim())');
  // pick the site with the most carriers below half of the called samples for display
  G.pick = $(`(()=>{let best=-1,bi=0; GEO.rows.forEach((r,i)=>{let c=0,n=0; for(const g of r.gts){ if(g!==3){n++; if(g>0)c++;} } const f=c/n; if(f<0.5 && c>best){best=c;bi=i;} }); return {idx:bi, pos:GEO.rows[bi].pos, carriers:best};})()`);
  $(`geoSelectSnp(${G.pick.idx})`); await site.wait(100);
  G.overview = $('[...document.querySelectorAll("#geoDetail tbody tr")].slice(0,8).map(tr=>[...tr.children].map(td=>td.textContent.trim()))');
  G.pngNA_freq = mapPNG(site, 'snpgeo_gene_Zm00001eb374090_NA_freq');
  $('geoSetColorMode("af")'); G.pngNA_af = mapPNG(site, 'snpgeo_gene_Zm00001eb374090_NA_af');
  G.fillsNA_af = $('[...document.querySelectorAll("#geoMap path.admin1")].filter(p=>p.style.fill && p.style.fill!=="rgb(228, 232, 238)" && p.style.fill!=="#e4e8ee").length');
  $('geoSetMapView("world"); geoSetColorMode("freq")'); await site.wait(50);
  G.worldPaths = $('document.querySelectorAll("#geoMap path.country").length');
  G.pngWorld_freq = mapPNG(site, 'snpgeo_gene_Zm00001eb374090_world_freq');
  $('geoShowDetail("USA")');
  G.usaDetail = $('document.getElementById("geoDetail").textContent.replace(/\\s+/g," ").trim().slice(0,600)');
  G.usaStates = $('[...document.querySelectorAll("#geoDetail .geo-region-table tbody tr")].slice(0,6).map(tr=>[...tr.children].map(td=>td.textContent.trim()))');
  snapshot(site, 'snpgeo_gene_Zm00001eb374090', 'SNPGeo - Zm00001eb374090');
  // full per-site statistics for the independent recomputation
  G.allStats = $(`JSON.parse(JSON.stringify(GEO.rows.map((r,i)=>{GEO.snpStats={}; const st=aggregateGeoData(i); const o={pos:r.pos, c:{}}; for(const [k,v] of Object.entries(st)) o.c[k]={total:v.total,called:v.called,count:v.count,ref:v.ref,missing:v.missing,het:v.het,altAlleles:v.altAlleles,calledAlleles:v.calledAlleles, states:Object.fromEntries(Object.entries(v.states).map(([s,x])=>[s,{total:x.total,called:x.called,count:x.count,het:x.het,altAlleles:x.altAlleles,calledAlleles:x.calledAlleles,code:x.code}]))}; return o;})))`);
  G.accOrder = $('GEO.accs.map(a=>a.id)');

  /* ---------- SNPGeo: gene with no store (chr3) -> graceful error ---------- */
  $('S.geoInput=null; go("snpgeo"); document.getElementById("geoGeneInput").value="Zm00001eb118470"');
  await $('geoLookupGene()');
  R.snpgeo_nostore = $('document.getElementById("geoLookupStatus").textContent');

  /* ---------- SNPVersity region (26 NAM lines, B73 + 25 founders) -> Send to SNPGeo ---------- */
  const V = R.versity_to_geo = {};
  $('S.dataset="zmgrin2026_imp"; S.chr="chr10"; S.start=9788000; S.end=9826500; S.selected=new Set(Data.defaultSelectionFor("zmgrin2026_imp")); go("snpversity")');
  await site.wait(100);
  await $('runQuery()');
  V.versityRows = $('S.results && S.results.rows.length'); V.versityAccs = $('S.results && S.results.accs.length');
  V.versityTableRows = $('document.querySelectorAll("#rtBody tbody tr").length');
  V.sendButton = $('[...document.querySelectorAll("button")].some(b=>/Send to SNPGeo/.test(b.textContent))');
  snapshot(site, 'snpversity_chr10_region', 'SNPVersity - chr10:9,788,000-9,826,500, 26 NAM lines');
  await $('sendToGeo()');
  await until(site, 'S.tool==="snpgeo" && document.querySelector("#geoMap svg")'); await site.wait(300);
  V.geo = {requeried: $('S.geoInput.requeried'), accs: $('S.geoInput.accs.length'), rows: $('GEO.rows.length'),
           header: $('document.querySelector(".page .card .n").textContent'), partialBanner: $('!!document.querySelector(".geo-partial")')};
  V.png = mapPNG(site, 'snpgeo_from_snpversity_chr10_NA_freq');
  snapshot(site, 'snpgeo_from_snpversity_chr10', 'SNPGeo - hand-off from SNPVersity');
  // partial hand-off (no re-query) must show the warning banner
  $('S.geoInput = Object.assign({}, S.geoInput, {accs: S.results.accs, rows: S.results.rows, requeried:false}); renderGeo(document.getElementById("page"))');
  await site.wait(300);
  V.partialBannerWhenPartial = $('!!document.querySelector(".geo-partial")');
  V.partialBannerText = $('(document.querySelector(".geo-partial")||{}).textContent||""').replace(/\s+/g, ' ').trim();

  /* ---------- GWAS Explorer -> SNPVersity (VCF-set dropdown) ---------- */
  const W = R.gwas = {};
  $('go("snpgwas")');
  await until(site, 'document.querySelector(".gwx-trait-dd-opt")', 15000);
  $('[...document.querySelectorAll(".gwx-trait-dd-opt")].find(b=>b.getAttribute("data-value")==="Tassel Branch Number").click()');
  await site.wait(100);
  await until(site, 'document.querySelector(\'button.gwx-chip[data-entry="tibbscortes2024_nam_tpbn_intcp"]\')');
  $('document.querySelector(\'button.gwx-chip[data-entry="tibbscortes2024_nam_tpbn_intcp"]\').click()');
  await until(site, 'SNPTools.registry.snpgwas.activeEntryId()==="tibbscortes2024_nam_tpbn_intcp" && document.getElementById("gwxOpenSendBtn")', 20000);
  W.entry = $('SNPTools.registry.snpgwas.activeEntryId()');
  const openPopup = async () => {
    $('SNPTools.registry.snpgwas.selectRegion(2, 4491424, 4499434)'); await site.wait(50);
    $('document.getElementById("gwxOpenSendBtn").click()'); await site.wait(50);
  };
  await openPopup();
  W.region = $('document.getElementById("gwxRegionCoord").textContent');
  W.regionStats = {total: $('document.getElementById("gwxStatTotal").textContent'), sig: $('document.getElementById("gwxStatSig").textContent')};
  W.options = $('[...document.querySelectorAll("#gwxSendDataset option")].map(o=>[o.value,o.textContent,o.selected])');
  W.defaultMapText = $('document.getElementById("gwxSendMap").textContent');
  W.replacePreticked = $('document.getElementById("gwxSendAccReplaceChk").checked');
  W.perSet = {};
  for (const ds of ['mgdb2026_hq', 'mgdb2026_hc', 'zmgrin2026_imp']){
    $(`(()=>{const s=document.getElementById("gwxSendDataset"); s.value="${ds}"; s.dispatchEvent(new Event("change"));})()`);
    W.perSet[ds] = {map: $('document.getElementById("gwxSendMap").textContent'),
                    missing: $('(document.querySelector("#gwxSendMap .gwx-ho-miss")||{}).textContent||""'),
                    replaceChecked: $('document.getElementById("gwxSendAccReplaceChk").checked')};
  }
  $('(()=>{const s=document.getElementById("gwxSendDataset"); s.value="mgdb2026_hq"; s.dispatchEvent(new Event("change"));})()');
  // (a) default VCF set: behaviour unchanged (all NAM runs of the MaizeGDB 2026 catalogue)
  $('(()=>{const c=document.getElementById("gwxSendAccReplaceChk"); c.checked=true; c.dispatchEvent(new Event("change")); document.getElementById("gwxSendConfirmBtn").click();})()');
  await site.wait(300);
  W.defaultSend = {tool: $('S.tool'), dataset: $('S.dataset'), selected: $('S.selected.size'), ids: $('[...S.selected].sort()'), chr: $('S.chr'), start: $('S.start'), end: $('S.end'),
                   banner: $('(document.getElementById("inboundBanner")||{}).textContent||""').replace(/\s+/g, ' ').trim().slice(0, 400),
                   founders: $('[...new Set([...S.selected].map(id=>(Data.accessionById(id)||{}).founder))].sort()')};
  // (b) GRIN-linked release: NAM names translated to ZmG_* sample ids
  $('go("snpgwas")'); await until(site, 'document.getElementById("gwxOpenSendBtn")', 10000);
  await openPopup();
  $('(()=>{const s=document.getElementById("gwxSendDataset"); s.value="zmgrin2026_imp"; s.dispatchEvent(new Event("change"));})()');
  W.grinMapText = $('document.getElementById("gwxSendMap").textContent');
  W.grinMissing = $('(document.querySelector("#gwxSendMap .gwx-ho-miss")||{}).textContent||""');
  $('(()=>{const c=document.getElementById("gwxSendAccReplaceChk"); c.checked=true; c.dispatchEvent(new Event("change"));})()');
  snapshot(site, 'gwas_send_popup_grin', 'GWAS Explorer - Send to SNPVersity (GRIN-linked VCF set)');
  $('document.getElementById("gwxSendConfirmBtn").click()'); await site.wait(300);
  W.grinSend = {tool: $('S.tool'), dataset: $('S.dataset'), selected: $('[...S.selected].sort()'), chr: $('S.chr'), start: $('S.start'), end: $('S.end'),
                banner: $('(document.getElementById("inboundBanner")||{}).textContent||""').replace(/\s+/g, ' ').trim().slice(0, 400)};
  W.pageHasNotAvailable = $('/not available in this set/i.test(document.getElementById("page").textContent)');
  await $('runQuery()');
  W.versity = {rows: $('S.results && S.results.rows.length'), accs: $('S.results && S.results.accs.length')};
  R.annot = {zmgrin2026_imp: annotStats(site)};
  snapshot(site, 'snpversity_from_gwas_grin', 'SNPVersity - GWAS region chr2:4,491,424-4,499,434, 26 NAM lines (GRIN-linked, release v1.4)');
  await $('sendToGeo()');
  await until(site, 'S.tool==="snpgeo" && document.querySelector("#geoMap svg")'); await site.wait(300);
  const idx = $('GEO.rows.findIndex(r=>r.pos===4494625)');
  W.geo = {requeried: $('S.geoInput.requeried'), accs: $('S.geoInput.accs.length'), rows: $('GEO.rows.length'), idx494625: idx};
  if (idx >= 0){ $(`geoSelectSnp(${idx})`); await site.wait(50);
    W.geo.usa = $('(()=>{const s=aggregateGeoData(GEO.snpIndex).USA; return s && {count:s.count,total:s.total,het:s.het,af:s.alleleFreq};})()');
    W.geo.nam = $('(()=>{const ids=new Set(Data.accessionsFor("zmgrin2026_imp").filter(a=>a.namFounder).map(a=>a.id)); const r=GEO.rows[GEO.snpIndex]; let c=0,n=0; GEO.accs.forEach((a,i)=>{ if(ids.has(a.id)&&r.gts[i]!==3){n++; if(r.gts[i]>0)c++;} }); return {carriers:c, called:n};})()');
    W.geo.png = mapPNG(site, 'snpgeo_from_gwas_tpbn_chr2_4494625_NA_freq'); }

  /* ---------- fresh page: GWAS -> SNPVersity on the default set (the path a user takes first) ---------- */
  {
    const f = await openSite(ROOT); const $f = f.$eval;
    const F = R.gwas_fresh = {selectedBefore: $f('S.selected.size'), datasetBefore: $f('S.dataset'),
      b73InDefault: $f('[...S.selected].some(id=>/^B73_/.test(id))')};
    $f('go("snpgwas")');
    const u = async (e, ms = 20000) => { const t0 = Date.now(); while (Date.now() - t0 < ms){ if ($f(e)) return; await f.wait(50); } };
    await u('document.querySelector(".gwx-trait-dd-opt")');
    $f('[...document.querySelectorAll(".gwx-trait-dd-opt")].find(b=>b.getAttribute("data-value")==="Tassel Branch Number").click()');
    await u('document.querySelector(\'button.gwx-chip[data-entry="tibbscortes2024_nam_tpbn_intcp"]\')');
    $f('document.querySelector(\'button.gwx-chip[data-entry="tibbscortes2024_nam_tpbn_intcp"]\').click()');
    await u('SNPTools.registry.snpgwas.activeEntryId()==="tibbscortes2024_nam_tpbn_intcp" && document.getElementById("gwxOpenSendBtn")');
    $f('SNPTools.registry.snpgwas.selectRegion(2, 4491424, 4499434)'); await f.wait(50);
    $f('document.getElementById("gwxOpenSendBtn").click()'); await f.wait(50);
    F.popupDefault = {value: $f('document.getElementById("gwxSendDataset").value'), replace: $f('document.getElementById("gwxSendAccReplaceChk").checked'),
                      map: $f('document.getElementById("gwxSendMap").textContent'),
                      options: $f('[...document.querySelectorAll("#gwxSendDataset option")].map(o=>o.textContent)')};
    $f('document.getElementById("gwxSendConfirmBtn").click()'); await f.wait(300);
    F.after = {tool: $f('S.tool'), dataset: $f('S.dataset'), selected: $f('S.selected.size'),
               b73: $f('[...S.selected].filter(id=>/^B73_/.test(id)).length'),
               lines: $f('new Set([...S.selected].map(id=>(Data.accessionById(id)||{}).founder)).size')};
    f.window.close();
  }

  /* ---------- SNPVersity annotation columns for MaizeGDB 2026 HQ / HC ----------
     No MaizeGDB 2026 HDF5 store exists locally, so the table is rendered from the real INFO
     of the MaizeGDB 2026 VCFs at the same chr2 window (sites only; fixtures/mgdb2026_*). The
     rows go through the same Data.parseVcf() + renderResults() path as a store query. */
  for (const [ds, fx] of [['mgdb2026_hq', 'mgdb2026_hq_chr2_4491424_4499434.sites.vcf.gz'],
                          ['mgdb2026_hc', 'mgdb2026_hc_chr2_4491424_4499434.sites.vcf.gz']]){
    const text = require('zlib').gunzipSync(fs.readFileSync(path.join(ROOT, 'localdev/fixtures', fx))).toString('utf8');
    site.window.__fixtureVcf = text;
    $(`go("snpversity"); S.dataset=${JSON.stringify(ds)}; S.chr="chr2"; S.start=4491424; S.end=4499434;
       (()=>{ const ids=Data.defaultSelectionFor(S.dataset).slice(0,5);
              const acc=new Map(Data.accessionsFor(S.dataset).map(a=>[a.id,a]));
              S.results={rows:Data.parseVcf(window.__fixtureVcf, ids).rows, accs:ids.map(id=>acc.get(id)||{id}), vcfUrl:'fixture:${ds}'};
              S.page=1; renderResults(); })()`);
    R.annot[ds] = annotStats(site);
    snapshot(site, 'snpversity_annotation_' + ds, 'SNPVersity - annotation columns, ' + ds + ' (MaizeGDB 2026 INFO, chr2 window)');
  }

  /* ---------- SNPCompare genome-wide scope + SNPTree genome-wide Newick ----------
     Synthetic distance files (make_synthetic_distance.py) in a temp dir, served by the real
     ibsCompare.php through SNPTOOLS_DISTANCE_DIR: zmgrin2026 (933 ids, all + SNP-only
     matrices, 2 trees) and the unchanged MaizeGDB 2026 layout (60 ids). Then the same page
     with an empty distance dir (genome-wide scope must be disabled for zmgrin2026). */
  {
    const os = require('os'), {execFileSync} = require('child_process');
    const dd = fs.mkdtempSync(path.join(os.tmpdir(), 'snpt-dist-')), empty = fs.mkdtempSync(path.join(os.tmpdir(), 'snpt-nodist-'));
    const gen = JSON.parse(execFileSync(process.env.PYTHON_PATH || 'python3', [path.join(__dirname, 'make_synthetic_distance.py'), dd, ROOT]).toString());
    const prevDD = process.env.SNPTOOLS_DISTANCE_DIR;
    const C = R.compare = {gen, distanceDir: dd};
    const cmpRun = async (s, ds, focal, sites) => {
      const $s = s.$eval;
      $s(`S.dataset=${JSON.stringify(ds)}; S.compareInput=null; go("snpcompare")`);
      await until(s, 'document.getElementById("cmpFocal")', 20000);
      const out = {gAvail: $s('SNPCompare.globalAvailable()'),
        scopeDisabled: $s('[...document.querySelectorAll("button")].filter(b=>b.textContent.trim()==="Genome-wide").map(b=>b.disabled)[0]'),
        sitesSelect: $s('!!document.getElementById("cmpGSites")'),
        note: $s('[...document.querySelectorAll(".mtx-note")].map(n=>n.textContent.replace(/\\s+/g," ").trim()).join(" | ")')};
      if (out.gAvail && focal){
        if (sites) $s(`SNPCompare.setGSites(${JSON.stringify(sites)})`);
        $s(`SNPCompare.pickFocal(${JSON.stringify(focal)}); SNPCompare.setMode("global"); SNPCompare.setView("table")`);
        await until(s, 'SNPCompare._ST.ran', 20000);
        out.focal = focal;
        out.rows = $s('SNPCompare._ST.allRows.map(r=>[r.id, r.gsim, r.gmiss])');
        out.tableRows = $s('document.querySelectorAll("#cmpViewWrap tbody tr").length');
        out.count = $s('(document.getElementById("cmpCount")||{}).textContent||""');
      }
      return out;
    };
    process.env.SNPTOOLS_DISTANCE_DIR = dd;
    const s2 = await openSite(ROOT);
    C.zmgrin_all = await cmpRun(s2, 'zmgrin2026_imp', 'ZmG_B73');
    C.zmgrin_snp = await cmpRun(s2, 'zmgrin2026_imp', 'ZmG_CML103', 'snp');
    snapshot(s2, 'snpcompare_zmgrin2026_global', 'SNPCompare - zmgrin2026_imp genome-wide (synthetic matrices)');
    const mFocal = fs.readFileSync(path.join(dd, 'ids.txt'), 'utf8').split('\n')[3];
    C.mgdb = await cmpRun(s2, 'mgdb2026_hq', mFocal);
    s2.$eval('S.dataset="zmgrin2026_imp"; S.treeInput=null; go("snptree")');
    await until(s2, 'document.querySelector("#treeGW a")', 20000);
    C.tree = {links: s2.$eval('[...document.querySelectorAll("#treeGW a")].map(a=>[a.textContent.trim(), a.getAttribute("href"), a.getAttribute("download")])'),
              nj: await s2.$eval('fetch(document.getElementById("treeGW_nj").getAttribute("href")).then(r=>r.text())'),
              upgma: await s2.$eval('fetch(document.getElementById("treeGW_upgma").getAttribute("href")).then(r=>r.text())')};
    C.requests = s2.log.filter(x => x.url === 'ibsCompare.php').map(x => x.query);
    s2.window.close();
    process.env.SNPTOOLS_DISTANCE_DIR = empty;
    const s3 = await openSite(ROOT);
    C.zmgrin_nofiles = await cmpRun(s3, 'zmgrin2026_imp', null);
    s3.$eval('S.dataset="zmgrin2026_imp"; S.treeInput=null; go("snptree")'); await s3.wait(300);
    C.tree_nofiles = s3.$eval('(document.getElementById("treeGW")||{}).innerHTML||""');
    s3.window.close();
    if (prevDD === undefined) delete process.env.SNPTOOLS_DISTANCE_DIR; else process.env.SNPTOOLS_DISTANCE_DIR = prevDD;
    // the installed files, when present (./distance/zmgrin2026/, not part of the tree)
    if (fs.existsSync(path.join(ROOT, 'distance/zmgrin2026/similarity.csv'))){
      const s4 = await openSite(ROOT);
      C.real = await cmpRun(s4, 'zmgrin2026_imp', 'ZmG_B73');
      C.real_mo17 = await cmpRun(s4, 'zmgrin2026_imp', 'ZmG_MO17');
      C.real_snp = await cmpRun(s4, 'zmgrin2026_imp', 'ZmG_B73', 'snp');
      snapshot(s4, 'snpcompare_zmgrin2026_global_real', 'SNPCompare - zmgrin2026_imp genome-wide, ZmG_B73 (installed matrices, SNPs only)');
      s4.$eval('S.dataset="zmgrin2026_imp"; S.treeInput=null; go("snptree")');
      await until(s4, 'document.querySelector("#treeGW a")', 20000);
      C.real_tree = {links: s4.$eval('[...document.querySelectorAll("#treeGW a")].map(a=>a.getAttribute("href"))'),
                     njBytes: (await s4.$eval('fetch("ibsCompare.php?tree=nj&dataset=zmgrin2026").then(r=>r.text())')).length};
      s4.window.close();
    }
  }

  /* ---------- query timing on a full-chromosome store (chr2), when installed ---------- */
  if ((R.store || {}).chr2 > 100000){
    const s5 = await openSite(ROOT); const $t = s5.$eval;
    const Tm = R.timing = {};
    const nam = $t('Data.accessionsFor("zmgrin2026_imp").filter(a=>a.namFounder).map(a=>a.id)');
    for (const [name, lo, hi] of [['gene_Zm00001eb067740', 4493424, 4497434], ['region_1Mb', 4000000, 5000000]]){
      for (const [sel, ids] of [['26_NAM', nam], ['all_933', null]]){
        const t0 = Date.now();
        const r = await $t(`Data.queryVariants("zmgrin2026_imp","chr2",${lo},${hi},${ids ? JSON.stringify(ids) : 'Data.accessionsFor("zmgrin2026_imp").map(a=>a.id)'})
          .then(r=>({rows:r.rows?r.rows.length:0, accs:r.accs?r.accs.length:0, wide:!!r.wide, variants:r.variants||null,
                     pc1:r.rows?r.rows.filter(x=>x.pc1!=null).length:0}))`);
        const t1 = Date.now();
        const php = s5.log.filter(x => x.url === 'processForm.php').slice(-1)[0] || {};
        Tm[name + '_' + sel] = Object.assign(r, {total_ms: t1 - t0, php_ms: php.ms, vcf_bytes: php.bytes});
      }
    }
    // table render for the 1-Mb x 26 NAM result (what SNPVersity draws)
    $t(`S.dataset="zmgrin2026_imp"; go("snpversity")`);
    const t2 = Date.now();
    await $t(`Data.queryVariants("zmgrin2026_imp","chr2",4000000,5000000,${JSON.stringify(nam)}).then(r=>{S.chr="chr2";S.start=4000000;S.end=5000000;S.results=r;S.page=1;renderResults();})`);
    Tm.region_1Mb_26_NAM_with_table_ms = Date.now() - t2;
    s5.window.close();
  }

  R.requests = site.log.filter(x => x.url && x.url.endsWith('.php')).map(x => ({url: x.url, ms: x.ms, n_genotypes: x.n_genotypes, reply: x.reply && x.reply.slice(0, 160)}));
  R.missingStatic = site.log.filter(x => x.status === 404).map(x => x.url);
  R.consoleErrors = site.errors;
  fs.writeFileSync(path.join(OUT, 'results.json'), JSON.stringify(R));
  const brief = JSON.parse(JSON.stringify(R)); delete brief.snpgeo_gene.allStats; delete brief.snpgeo_gene.accOrder;
  brief.gwas.grinSend.selected = brief.gwas.grinSend.selected.length; delete brief.gwas.defaultSend.ids;
  brief.snptrait.filterSS_Ames_Dent = brief.snptrait.filterSS_Ames_Dent.length; brief.snptrait.rangeKW = brief.snptrait.rangeKW.length; delete brief.snptrait.facetCounts; brief.snptrait.filterPlusIowa = brief.snptrait.filterPlusIowa.length;
  if (brief.compare) Object.values(brief.compare).forEach(v => { if (v && Array.isArray(v.rows)) { v.top3 = v.rows.slice().sort((x, y) => y[1] - x[1]).slice(0, 3); v.rows = v.rows.length; } });
  fs.writeFileSync(path.join(OUT, 'results_brief.json'), JSON.stringify(brief, null, 1));
  site.window.close();
  console.log('done');
})().catch(e => { console.error('HARNESS FAILURE', e); process.exit(1); });
