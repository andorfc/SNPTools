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
async function until(site, expr, ms = 8000){
  const t0 = Date.now();
  while (Date.now() - t0 < ms){ if (site.$eval(expr)) return true; await site.wait(50); }
  return false;
}

(async () => {
  const site = await openSite(ROOT);
  const $ = site.$eval;
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

  /* ---------- SNPVersity region (25 NAM founders) -> Send to SNPGeo ---------- */
  const V = R.versity_to_geo = {};
  $('S.dataset="zmgrin2026_imp"; S.chr="chr10"; S.start=9788000; S.end=9826500; S.selected=new Set(Data.defaultSelectionFor("zmgrin2026_imp")); go("snpversity")');
  await site.wait(100);
  await $('runQuery()');
  V.versityRows = $('S.results && S.results.rows.length'); V.versityAccs = $('S.results && S.results.accs.length');
  V.versityTableRows = $('document.querySelectorAll("#rtBody tbody tr").length');
  V.sendButton = $('[...document.querySelectorAll("button")].some(b=>/Send to SNPGeo/.test(b.textContent))');
  snapshot(site, 'snpversity_chr10_region', 'SNPVersity - chr10:9,788,000-9,826,500, 25 NAM founders');
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

  R.requests = site.log.filter(x => x.url && x.url.endsWith('.php')).map(x => ({url: x.url, ms: x.ms, n_genotypes: x.n_genotypes, reply: x.reply && x.reply.slice(0, 160)}));
  R.missingStatic = site.log.filter(x => x.status === 404).map(x => x.url);
  R.consoleErrors = site.errors;
  fs.writeFileSync(path.join(OUT, 'results.json'), JSON.stringify(R));
  const brief = JSON.parse(JSON.stringify(R)); delete brief.snpgeo_gene.allStats; delete brief.snpgeo_gene.accOrder;
  brief.snptrait.filterSS_Ames_Dent = brief.snptrait.filterSS_Ames_Dent.length; brief.snptrait.rangeKW = brief.snptrait.rangeKW.length; delete brief.snptrait.facetCounts; brief.snptrait.filterPlusIowa = brief.snptrait.filterPlusIowa.length;
  fs.writeFileSync(path.join(OUT, 'results_brief.json'), JSON.stringify(brief, null, 1));
  site.window.close();
  console.log('done');
})().catch(e => { console.error('HARNESS FAILURE', e); process.exit(1); });
