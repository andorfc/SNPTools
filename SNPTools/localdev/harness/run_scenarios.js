/* run_scenarios.js <siteRoot> <outDir> -- end-to-end checks of SNPTrait and SNPGeo
 * on the zmgrin2026_imp dataset (the only set offered in the initial release); writes results.json, map PNG/SVG renders
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
    const filled={}, decimals={}, sample={};
    cols.forEach((c,j)=>{ let n=0, d=0; const ex=[]; host.querySelectorAll('tr').forEach(tr=>{ const td=tr.children[first+j];
      const v=td?td.textContent.trim():''; if(v && !/^(—|N\\/A|NA|\\.|intergenic)$/.test(v)){ n++;
        if(/^-?[0-9.]+$/.test(v)){ d=Math.max(d,(v.split('.')[1]||'').length); if(ex.length<5) ex.push(v); } } });
      filled[c.key]=n; decimals[c.key]=d; sample[c.key]=ex; });
    return {dataset:S.dataset, headers:th.slice(0, first+cols.length), cols, rows:S.results.rows.length, filled, decimals, sample,
            fields:Data.annotationFields(S.dataset).map(f=>({key:f.key, label:f.label, status:f.status})),
            cellsPerRow:host.querySelector('tr') ? host.querySelector('tr').children.length : 0, headerCells:document.querySelectorAll('#rtBody table.vcf thead th').length,
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
            defaultDataset: $('S.dataset')};   // initial release: zmgrin2026_imp only (MaizeGDB 2026 commented out)

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
  /* the Send dialog, as in GWAS Explorer: replace pre-ticked; "add" keeps SNPVersity's
     selection (here 5 lines SNPTrait has not selected, so add and replace differ) and adds
     the SNPTrait lines to it */
  $('go("snptrait"); S.selected=new Set(Data.accessionsFor("zmgrin2026_imp").map(a=>a.id).filter(id=>!TRAIT.selected.has(id)).slice(0,5)); traitRenderRunbar()');
  T.dialogBefore = {versity: $('[...S.selected].sort()'), trait: $('[...TRAIT.selected].sort()'),
                    button: $('[...document.querySelectorAll("#traitRunbar button")].map(b=>b.textContent.trim()).find(t=>/SNPVersity/.test(t))')};
  $('[...document.querySelectorAll("#traitRunbar button")].find(b=>/SNPVersity/.test(b.textContent)).click()');
  T.dialog = {open: $('!!document.querySelector("#traitSendDlg .trait-send")'), text: $('document.querySelector(".trait-send").textContent.replace(/\\s+/g," ").trim()'),
              replace: $('document.getElementById("traitSendReplace").checked'), add: $('document.getElementById("traitSendAdd").checked')};
  $('document.getElementById("traitSendAdd").click()');
  T.dialog.afterAdd = {replace: $('document.getElementById("traitSendReplace").checked'), add: $('document.getElementById("traitSendAdd").checked'),
                       sendDisabled: $('document.getElementById("traitSendGo").disabled')};
  snapshot(site, 'snptrait_send_dialog', 'SNPTrait - Send to SNPVersity (add to the current selection)');
  $('document.getElementById("traitSendGo").click()'); await site.wait(300);
  T.handoffAdd = {tool: $('S.tool'), dialogClosed: $('!document.getElementById("traitSendDlg")'), selected: $('[...S.selected].sort()'),
                  banner: $('(document.getElementById("inboundBanner")||{}).textContent||""').replace(/\s+/g, ' ').trim().slice(0, 300)};
  /* A family without a generated schema falls back to the neutral schema: no listed dataset
     lacks one while only the GRIN-linked set is offered, so the MaizeGDB 2026 case is off.
  $('S.dataset="mgdb2026_hq"; go("snptrait")');
  T.mgdb2026 = {groupBy: $('TRAIT.schema.groupBy'), facets: $('TRAIT.schema.facets.map(f=>f[0])'), rows: $('TRAIT.rows.length'), header: $('document.querySelector("#page h2").textContent')};
  */

  /* ---------- SNPVersity step 3: quick picks, the SNPTrait drawer, a pasted list ----------
     Step 3 shows the selection (count, panel mix, chips) with one-click picks and a pasted or
     uploaded list; "Browse & filter lines…" opens SNPTrait over the page on a draft of the
     selection (TraitDrawer). check_snptrait.py recomputes every expected set from the catalogue. */
  {
    const Q = R.step3 = {};
    // the SNPTrait page's own selection, which the drawer must leave alone
    $('go("snptrait"); traitClearFacets(); TRAIT.q=""; TRAIT.onlySel=false; TRAIT.selected=new Set(["ZmG_B73"]); traitRenderTable()');
    Q.pageSel = $('[...TRAIT.selected]');
    $('S.selected=new Set(Data.defaultSelectionFor("zmgrin2026_imp")); go("snpversity")'); await site.wait(100);
    Q.before = $('[...S.selected].sort()');
    Q.initial = {count: $('document.getElementById("selCount").textContent'), oldPicker: $('!!document.getElementById("accList")'),
                 picks: $('[...document.querySelectorAll("#accPicks .qbtn")].map(b=>[b.textContent.replace(/\\s+/g," ").trim(), b.classList.contains("on")])'),
                 mix: $('(document.querySelector(".acc2-legend")||{}).textContent||""').replace(/\s+/g, ' ').trim()};
    snapshot(site, 'snpversity_step3', 'SNPVersity - step 3: selection summary and quick picks');
    // a quick pick replaces the selection and offers Undo
    $('applyPick("inAmes282")');
    Q.pick = {selected: $('[...S.selected].sort()'), on: $('[...document.querySelectorAll("#accPicks .qbtn.on")].map(b=>b.textContent.replace(/\\s+/g," ").trim())'),
              say: $('document.getElementById("accSay").textContent').replace(/\s+/g, ' ').trim()};
    $('document.querySelector("#accSay .ho-undo").click()');
    Q.undo = $('[...S.selected].sort()');
    // the drawer: SNPTrait on a draft of the selection; Apply writes it back
    $('openLineSelector()');
    Q.drawer = {open: $('!!document.querySelector("#traitDrawer .tdr")'), mode: $('TRAIT.mode'), draft: $('[...TRAIT.selected].sort()'),
                facetBlocks: $('document.querySelectorAll("#traitDrawer .facet-block").length')};
    $('traitToggleFacet("subpop","SS",true); traitToggleFacet("inAmes282","yes",true); traitToggleFacet("kernelType","Dent",true); traitSelectVisible(true)');
    Q.drawer.visible = $('traitVisible().map(r=>r.id).sort()');
    Q.drawer.foot = $('document.getElementById("tdrFoot").textContent').replace(/\s+/g, ' ').trim();
    snapshot(site, 'snpversity_drawer', 'SNPVersity - Browse & filter lines (SNPTrait drawer)');
    $('document.getElementById("tdrApply").click()');
    Q.applied = {selected: $('[...S.selected].sort()'), closed: $('!document.getElementById("traitDrawer")'), mode: $('TRAIT.mode'),
                 say: $('document.getElementById("accSay").textContent').replace(/\s+/g, ' ').trim()};
    // Escape cancels; the drawer keeps its own filters for the next opening
    $('openLineSelector()');
    Q.reopenFilters = $('[...(TRAIT.facets.subpop||[])]');
    $('traitSelectGroup("WiDiv",true); document.dispatchEvent(new window.KeyboardEvent("keydown",{key:"Escape"}))');
    Q.afterCancel = {selected: $('[...S.selected].sort()'), closed: $('!document.getElementById("traitDrawer")')};
    $('go("snptrait")');
    Q.pageAfter = {selected: $('[...TRAIT.selected]'), filtered: $('Object.values(TRAIT.facets).some(x=>x.size)')};
    // a pasted list: GRIN accession numbers and line names resolve
    $('go("snpversity")'); await site.wait(100);
    $('Handoff.setMode("replace"); applyAccList("PI 550473\\nMo17\\nnot-a-line", "pasted list")');
    Q.paste = {selected: $('[...S.selected].sort()'), unmatched: $('[...document.querySelectorAll("#uplReport li .mono")].map(e=>e.textContent)')};
    $('Handoff.setMode("add")');
  }

  /* ---------- the top bar's selection chip ----------
     SNPVersity's selection in every tool; it opens the same drawer. Outside SNPVersity, Apply
     sets the selection and a toast offers Undo. On the SNPTrait page the page's nodes wait
     outside the document while the drawer (same element ids) is open. */
  {
    const C = R.chip = {};
    $('S.selected=new Set(Data.defaultSelectionFor("zmgrin2026_imp")); go("snpgeo")');
    C.onGeo = $('document.getElementById("selChipN").textContent');
    $('document.getElementById("selChip").click()');
    C.open = {drawer: $('!!document.querySelector("#traitDrawer .tdr")'), title: $('document.getElementById("tdrTitle").textContent'), draft: $('TRAIT.selected.size')};
    $('traitClearFacets(); TRAIT.onlySel=false; traitSelectGroup("Other GRIN", true); document.getElementById("tdrApply").click()');
    C.applied = {n: $('S.selected.size'), chip: $('document.getElementById("selChipN").textContent'), tool: $('S.tool'),
                 toast: $('document.getElementById("snpToast").textContent').replace(/\s+/g, ' ').trim()};
    $('[...document.querySelectorAll("#snpToast button.act")].find(b=>b.textContent==="Undo").click()');
    C.undone = {n: $('S.selected.size'), chip: $('document.getElementById("selChipN").textContent')};
    $('go("snptrait")');
    $('document.getElementById("selChip").click()');
    C.onTrait = {grids: $('document.querySelectorAll("#traitGrid").length'), pageSetAside: $('!document.querySelector("#page #traitGrid")')};
    $('TraitDrawer.close(false)');
    C.traitBack = $('!!document.querySelector("#page #traitGrid")');
  }

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
  for (const ds of $('[...document.querySelectorAll("#gwxSendDataset option")].map(o=>o.value)')){   // the sets offered
    $(`(()=>{const s=document.getElementById("gwxSendDataset"); s.value="${ds}"; s.dispatchEvent(new Event("change"));})()`);
    W.perSet[ds] = {map: $('document.getElementById("gwxSendMap").textContent'),
                    missing: $('(document.querySelector("#gwxSendMap .gwx-ho-miss")||{}).textContent||""'),
                    replaceChecked: $('document.getElementById("gwxSendAccReplaceChk").checked')};
  }
  $('(()=>{const s=document.getElementById("gwxSendDataset"); s.value=s.options[0].value; s.dispatchEvent(new Event("change"));})()');
  // (a) default VCF set: the first one offered (the GRIN-linked set in this release; MaizeGDB 2026 HQ before)
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
      b73InDefault: $f('[...S.selected].some(id=>/^(B73_|ZmG_B73$)/.test(id))')};
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
               b73: $f('[...S.selected].filter(id=>/^(B73_|ZmG_B73$)/.test(id)).length'),
               lines: $f('new Set([...S.selected].map(id=>(Data.accessionById(id)||{}).founder)).size')};
    f.window.close();
  }

  /* ---------- SNPVersity annotation columns for MaizeGDB 2026 HQ / HC ----------
     No MaizeGDB 2026 HDF5 store exists locally, so the table is rendered from the real INFO
     of the MaizeGDB 2026 VCFs at the same chr2 window (sites only; fixtures/mgdb2026_*). The
     rows go through the same Data.parseVcf() + renderResults() path as a store query.
     Off while the initial release offers the GRIN-linked set only (uncomment with the sets). */
  if (false) for (const [ds, fx] of [['mgdb2026_hq', 'mgdb2026_hq_chr2_4491424_4499434.sites.vcf.gz'],
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
    /* MaizeGDB 2026's legacy layout (./distance/maizegdb_allchr_final_*): off while that set
       is not offered; an unknown dataset id now resolves to the GRIN-linked family.
    const mFocal = fs.readFileSync(path.join(dd, 'ids.txt'), 'utf8').split('\n')[3];
    C.mgdb = await cmpRun(s2, 'mgdb2026_hq', mFocal);
    */
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

  /* ---------- one gene's consequences at multi-gene sites (GRIN-linked) ----------
     SnpEff lists one consequence per gene; 1,464 of the 3,744 fixture sites list several. SNPFunction
     (Data.geneFunction), SNPFold (Data.queryFoldVariants) and SNPGeo's gene search must read the
     gene's own entry, with ESM only where it was computed for that entry; SNPVersity marks the
     other entries (+N) and lists every gene. check_gene_consequences.py recomputes all of it
     from the fixture VCFs and the ESM table. */
  {
    const P = R.pergene = {genes: {}};
    // + a chr2 gene when the full chr2 store (with Evo2 and ESM-C) is installed
    for (const gene of ['Zm00001eb374230', 'Zm00001eb404750', 'Zm00001eb374090', 'Zm00001eb056510'].concat((R.store || {}).chr2 > 100000 ? ['Zm00001eb067740'] : [])){
      const g = JSON.stringify(gene);
      P.genes[gene] = {
        interval: await $(`Data.lookupGene(${g})`),
        fn: await $(`Data.geneFunction(${g}, "zmgrin2026_imp").then(d=>({n:d.nVariants, nAcc:d.nAccessions, byClass:d.burden.byClass, damaging:d.damaging.length,
              v:d.variants.map(v=>({pos:v.pos, ref:v.ref, alt:v.alt, cls:v.consClass, sub:v.resi!=null?(v.aaRef||'')+v.resi+(v.aaAlt||''):null, esm:v.esm, esm2:v.esm2, esm3:v.esm3, evo2:v.evo2, esmc:v.esmc, hom:v.hom, het:v.het})),
              means:{evo2:d.burden.meanEvo2, esmc:d.burden.meanEsmc}}))`),
        fold: await $(`Data.queryFoldVariants(${g}, "zmgrin2026_imp").then(a=>a.map(v=>({pos:v.pos, resi:v.resi, variant:v.variant, cls:v.consClass, esm:v.esm, evo2:v.evo2, esmc:v.esmc})))`)};
    }
    $('S.dataset="zmgrin2026_imp"; S.geoInput=null; go("snpgeo"); document.getElementById("geoGeneInput").value="Zm00001eb374230"');
    await $('geoLookupGene()');
    // geoLookupGene() re-renders the page 200 ms later; let that land before leaving SNPGeo
    await until(site, 'GEO.rows && GEO.rows.length && document.querySelector("#geoMap svg")'); await site.wait(300);
    P.geo = $('JSON.parse(JSON.stringify(GEO.rows.map(r=>({pos:r.pos, gene:r.gene, effect:r.effect, sub:r.sub, esm1:r.esm1}))))');
    $('S.dataset="zmgrin2026_imp"; S.chr="chr9"; S.start=13118306; S.end=13124164; S.selected=new Set(Data.defaultSelectionFor("zmgrin2026_imp")); go("snpversity")');
    await site.wait(100);
    await $('runQuery()');
    P.versity = $(`(()=>{ const host=document.createElement('tbody'); host.innerHTML=S.results.rows.map(r=>rowHTML(r)).join('');
      return {window:[S.chr,S.start,S.end], rows:S.results.rows.length,
              markers:[...host.querySelectorAll('tr')].map(tr=>{ const m=tr.querySelector('.gene-more'); return m ? m.textContent : ''; }),
              genes:[...document.querySelectorAll('details .c-mono')].map(e=>e.textContent.trim())}; })()`);
    snapshot(site, 'snpversity_chr9_multigene', 'SNPVersity - chr9:13,118,306-13,124,164, sites with several consequences (+N)');
  }

  R.requests = site.log.filter(x => x.url && x.url.endsWith('.php')).map(x => ({url: x.url, ms: x.ms, n_genotypes: x.n_genotypes, reply: x.reply && x.reply.slice(0, 160)}));
  R.missingStatic = site.log.filter(x => x.status === 404).map(x => x.url);
  R.consoleErrors = site.errors;
  fs.writeFileSync(path.join(OUT, 'results.json'), JSON.stringify(R));
  const brief = JSON.parse(JSON.stringify(R)); delete brief.snpgeo_gene.allStats; delete brief.snpgeo_gene.accOrder;
  brief.gwas.grinSend.selected = brief.gwas.grinSend.selected.length; delete brief.gwas.defaultSend.ids;
  if (brief.pergene){ delete brief.pergene.geo; Object.values(brief.pergene.genes).forEach(g => { g.fn.v = g.fn.v.length; g.fold = g.fold.length; });
    brief.pergene.versity.markers = brief.pergene.versity.markers.filter(Boolean).length; }
  brief.snptrait.filterSS_Ames_Dent = brief.snptrait.filterSS_Ames_Dent.length; brief.snptrait.rangeKW = brief.snptrait.rangeKW.length; delete brief.snptrait.facetCounts; brief.snptrait.filterPlusIowa = brief.snptrait.filterPlusIowa.length;
  if (brief.compare) Object.values(brief.compare).forEach(v => { if (v && Array.isArray(v.rows)) { v.top3 = v.rows.slice().sort((x, y) => y[1] - x[1]).slice(0, 3); v.rows = v.rows.length; } });
  fs.writeFileSync(path.join(OUT, 'results_brief.json'), JSON.stringify(brief, null, 1));
  site.window.close();
  console.log('done');
})().catch(e => { console.error('HARNESS FAILURE', e); process.exit(1); });
