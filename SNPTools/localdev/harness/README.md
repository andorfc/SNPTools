# localdev/harness — headless end-to-end checks

`run_scenarios.js` loads `index.html` in jsdom with every `<script>` in page order.
Static files are served from the site root; each `*.php` request is executed by the real
PHP script through the PHP CLI (`php_shim.php`), so `processForm.php -> h5_to_vcf.py -> HDF5`
and `lookupGeneModel.php` run exactly as under `php -S`. No browser and no listening socket
are needed. d3 v5 is served from `node_modules`.

Scenarios (dataset `zmgrin2026_imp`, local test store from `../build_test_store.sh`):
1. SNPTrait: schema, facet counts, compound filter (SS x Ames282 x Dent), search, a
   numeric trait range (1000-kernel weight 250-300 g, then x Ames282), select-visible
   and hand-off to SNPVersity; neutral schema for `mgdb2026_hq`.
2. Help page lists SNPTrait (live) and SNPGeo.
3. SNPGeo gene search `Zm00001eb374090` (all 932 lines): table, map in North America and
   world views, three colour modes, country detail; per-site statistics dumped.
4. SNPGeo gene on a chromosome without a store (graceful error).
5. SNPVersity chr10:9,788,000-9,826,500 with the 25 NAM founders -> Send to SNPGeo
   (re-query of all 932 lines), then a forced partial hand-off (warning banner).
6. GWAS Explorer, Tassel Branch Number (intercept, NAM, Tibbs-Cortes et al. 2024): region
   chr2:4,491,424-4,499,434 (Zm00001eb067740 +-2 kb; 85 GWAS SNPs, 11 significant) -> Send to
   SNPVersity, (a) with the default VCF set (unchanged: 81 MaizeGDB 2026 NAM runs) and (b) with
   the GRIN-linked set (NAM names translated to 25 ZmG_* samples, CML103 reported as not
   available), then SNPVersity query and SNPGeo. The canvas is a no-op stub (nothing is drawn).

Outputs in `../out/`: `results.json`, `results_brief.json`, map PNG/SVG files (the app's own
`geoBuildExportSVG()` export, rasterised with resvg) and static HTML snapshots of each page.
`check_snptrait.py`, `check_snpgeo_counts.py` and `check_gwas_handoff.py` recompute the numbers independently from
the catalogue, trait side-file, region records and fixture VCF.

    cd SNPTools/localdev/harness
    npm install
    PHP_BIN=$(which php) PYTHON_PATH=/path/to/python-with-h5py npm test
