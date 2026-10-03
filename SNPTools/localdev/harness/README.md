# localdev/harness — headless end-to-end checks

`run_scenarios.js` loads `index.html` in jsdom with every `<script>` in page order.
Static files are served from the site root; each `*.php` request is executed by the real
PHP script through the PHP CLI (`php_shim.php`), so `processForm.php -> h5_to_vcf.py -> HDF5`
and `lookupGeneModel.php` run exactly as under `php -S`. No browser and no listening socket
are needed. d3 v5 is served from `node_modules`.

Scenarios (dataset `zmgrin2026_imp`, the only set offered in the initial release, local test store
from `../build_test_store.sh`; the MaizeGDB 2026 cases are kept in the code, switched off):
1. SNPTrait: schema, facet counts, compound filter (SS x Ames282 x Dent), search, a
   numeric trait range (1000-kernel weight 250-300 g, then x Ames282), select-visible
   and hand-off to SNPVersity (replace); while SS x Ames282 x Dent is ticked every section's counts
   are taken under the other sections' filters (Subpopulation still lists its values) and SS or NSS
   is the union; then the Send dialog: replace pre-ticked, add unticks
   it, and *add* keeps SNPVersity's 5 lines and adds the 14 (19). SNPTrait starts with nothing
   selected (not SNPVersity's selection). (The neutral-schema case for
   `mgdb2026_hq` is off with that set.)
   Then SNPVersity's step 3 (prototype): the quick pick Ames 282 (280 = inAmes282) with Undo, the
   SNPTrait drawer (draft = the 26 NAM lines; SS x Ames282 x Dent, select visible, Apply -> 39;
   Escape changes nothing; filters kept on reopening; the SNPTrait page's selection untouched), and
   a pasted list resolving `PI 550473` and `Mo17`.
   Then the top bar's selection chip: 26 on SNPGeo, opened there, + Other GRIN applied (56) with
   the toast, Undo (26), and opened on the SNPTrait page (the page's nodes set aside, one #traitGrid
   in the document, restored on close).
2. Help page lists SNPTrait (live) and SNPGeo.
3. SNPGeo gene search `Zm00001eb374090` (all 933 lines): table, map in North America and
   world views, three colour modes, country detail; per-site statistics dumped.
4. SNPGeo gene on a chromosome without a store (graceful error).
5. SNPVersity chr10:9,788,000-9,826,500 with the 26 NAM lines (B73 + 25 founders) -> Send to SNPGeo
   (re-query of all 933 lines), then a forced partial hand-off (warning banner).
6. GWAS Explorer, Tassel Branch Number (intercept, NAM, Tibbs-Cortes et al. 2024): region
   chr2:4,491,424-4,499,434 (Zm00001eb067740 +-2 kb; 85 GWAS SNPs, 11 significant) -> Send to
   SNPVersity, (a) with the default VCF set (the first offered: the GRIN-linked set, 26 ZmG_* NAM
   samples; it was MaizeGDB 2026 HQ, 81 runs, while that set was listed) and (b) with
   the GRIN-linked set (release v1.4: all 26 NAM names translated to 26 ZmG_* samples, none
   reported as not available), then SNPVersity query and SNPGeo. The canvas is a no-op stub (nothing is drawn).

7. SNPVersity annotation columns: the 14 columns (Gene model ... ESM-C, with Site QC after MAF when the
   chr2 store has its site-QC sidecar) for the GRIN-linked set
   (chr2 GWAS window from the store; MaizeGDB 2026 HQ and HC, from `../fixtures/mgdb2026_*`, only
   when those sets are offered); MQ and COMP are absent for the GRIN-linked set ('hidden') and
   every row has as many cells as the header; the Domain column is expected filled wherever `data/domains/`
   covers the site (it assumed the files absent and failed whenever they were installed);
   `check_annotation_columns.py` recomputes the filled-cell counts from the INFO, checks that
   unavailable/pending columns stay in the table empty, and checks the GRIN-linked INFO against
   its sources (SnpEff fields, MAF from the genotypes, ESM1/2/3 from the missense ESM table).

8. SNPCompare genome-wide scope per dataset family, with synthetic matrices written by
   `make_synthetic_distance.py` to a temp dir (`SNPTOOLS_DISTANCE_DIR`): zmgrin2026 (933 ids, all +
   SNP-only), the MaizeGDB 2026 legacy layout (only when that set is offered), an empty
   dir (scope disabled), and SNPTree's genome-wide Newick downloads; then, if installed, the real
   `distance/zmgrin2026/` files (`check_global_compare.py`, `check_ibs_files.py`: ids = catalogue,
   shapes, symmetry, diagonals, fractions, tree tips).
9. One gene's consequences at sites SnpEff annotates in several genes (1,464 of the 3,744 fixture sites):
   SNPFunction's burden and allele catalog (over usable sites: Site QC PASS or HET_ELEVATED, with the
   flagged and no-carrier damaging alleles counted apart) and SNPFold's coding variants for Zm00001eb374230,
   Zm00001eb404750, Zm00001eb374090 and Zm00001eb056510; SNPGeo's gene search for
   Zm00001eb374230 (site 13,120,567 reads A471G, not the neighbour's R65P); SNPVersity's "+N"
   markers and gene list on chr9:13,118,306-13,124,164. `check_gene_consequences.py` recomputes
   all of it from the fixture VCFs (each gene's own entry, carriers from the genotypes) and checks
   every ESM score shown against the ESM table row for that gene's substitution. With the full chr2
   store, Zm00001eb067740 is added: Evo2 (allele-level, every entry) and ESM-C (with the ESM entry)
   per variant, the burden means, and SNPFold's values; the annotation-column check counts the
   Evo2 / ESM-C cells against the INFO.
10. With a full chr2 store installed: the chr2 checks use its window
   (`../fixtures/chr2_store/`), and gene / 1-Mb query timings are recorded in `results.json`.
   With full stores for the other chromosomes, `check_gene_consequences.py` reads each fixture
   window (and each tested gene's whole interval) from the store through `store_windows.py`, which
   first confirms the store matches the release fixture there: same sites and REF/ALT, same
   genotypes, same value for every INFO key the fixture has. The scores the fixtures predate
   (PlantCAD, Evo2, ESM-C, maxR²) then come from the store; the ESM table check applies at fixture
   sites.
11. The code-review fixes (patch 0038), `check_review_fixes.py`: SNPTrait's Random % with no line
   visible; `randomSample()` drawing distinct ids flat over catalogue order; a SNPVersity result
   keeping its queried region (CHR column, JBrowse link, heading, Send to SNPTree) after the form
   moves to another chromosome; the newer of two Builds, two SNPFunction genes, two SNPGeo searches
   and two SNPFold genes winning when the older answers last (stubbed data calls resolved out of
   order, since PHP runs synchronously here); a SNPGeo search answered after leaving SNPGeo not
   drawing over SNPVersity; one variant query and one gene lookup per SNPFold gene. Then, in a
   temporary site root (the real `vcf/` is not pruned): `processForm.php` refusing a request over
   `SNPTOOLS_MAX_CELLS`, keeping string ids once each, filling a reserved name (INFO) with ./.,
   never reusing a file name, ignoring a path in `outName`, refusing `start=1e5`, pruning its own
   VCFs past the TTL; `ibsCompare.php` rebuilding a torn `.offidx`. Gene lookups (patch 0039):
   all 39,756 genes read back from `gff/genes_index.txt` exactly as from the store, misses stay
   misses; `lookupGeneModel.php` leaves a current index alone, rebuilds a stale, missing or
   truncated one, and reads the store itself when the index is stale and `gff/` is read-only.
12. Site QC counts in every VCF, `check_site_qc.py`. The page's query path (Data.queryVariants ->
   processForm.php -> h5_to_vcf.py) for y1 (chr6:91,643,082-91,646,759) with five lines and with all
   933, su1 (chr4:43,430,007-43,438,753) and chr10:96,075,873-96,278,559 with all 933 (two fixture
   windows when only the test stores are built), each once as is and once with `SNPTOOLS_SITEQC=0` for
   the PHP process. The checker recomputes NHET / NHOM from each row's 933 genotypes and SITEQC from the
   rule; the five-line INFO must equal the all-lines INFO (panel-wide counts); the `SNPTOOLS_SITEQC=0`
   VCF must equal the normal one minus the three fields and their ##INFO lines. In a temporary root
   that symlinks the store: no sidecar and a sidecar with a wrong `store_bytes` each give a `Note:` line
   and the `SNPTOOLS_SITEQC=0` VCF, a valid copy gives the normal one. The Python copies of the rule
   (`build_site_qc.py`, `h5_to_vcf.py`, `annotate_release_info.py`) must agree on every (het, hom) pair
   from 0 to 60. `annotate_release_info.py` is run on every fixture VCF with the inputs in
   `../fixtures/annotation/`: its counts must equal the fixture's genotypes and the sidecar at the same
   sites, every INFO key of the fixture keeps its value (`ESMC_score`, which the fixtures predate, is
   added), and a store built from its chr9 output gets no second set of the three keys from
   `h5_to_vcf.py`, with or without a sidecar.
13. Site QC in the pages, `check_site_qc_ui.py` (needs the full chr2/chr4/chr6/chr8/chr10 stores;
   skipped otherwise). SNPVersity on eight gene intervals (y1, su1, ZmWAK, Bx13, DGAT1-2, Htn1, tga1,
   crtRB1) with the 26 NAM lines: each row's class, the rows under the three filter choices, the
   summary line, the column after MAF and its cell text. SNPFunction on the same genes, read from the
   page: "<v> of <n> sites vary", the allele list and its two groups with their pills, knockout lines and
   those only in flagged calls, the burden's usable-site count, the banner; recomputed from the
   933-line VCF its query wrote. SNPImpact on su1 under each choice with the hidden-count note,
   SNPFold's "Usable alleles only" on Zm00001eb406050, SNPGeo's pills, and SNPTree / SNPMatrix /
   SNPCompare sites used on chr6:91,593,082-91,793,082 under both settings. The JS rule
   (`Data.siteQc`) against the rule on all 3,721 (het, hom) pairs 0-60. Then everything again with
   `SNPTOOLS_SITEQC=0`: no Site QC column, select, note, pill or site control, while SNPFunction (and
   SNPFold, built on it) give the same lists from their own carrier counts.

14. Line QC (`js/zmgrin.lineqc.js`), in `check_snptrait.py`: SNPTrait's Sample heterozygosity facet
   counts (Inbred 887, Elevated heterozygosity 26, Heterozygous sample 20) against the file, the 20
   heterozygous samples sent to SNPVersity, `sampleQC` and `hetShare` in both exports, the "het" marker
   (with its tooltip) on exactly those lines in SNPVersity's header cells and selection chips and in
   SNPFunction's carrier chips (ZmG_CH9 checked by name); then the same pages from a temporary root of
   symlinks without the file: no facet, no marker, no error. With the file loaded, scenario 1's facet
   checks include the added facet.

15. SNPCurate, `check_snpcurate.py` (needs the full stores): the table, the filters and the two sorts,
   every record as drawn (facts, the seven scores with their rank in the gene, subpopulation, country
   and trait rows, the reason a grey allele cannot be shown, references as DOI links) and the three
   buttons of SC0001. The checker recomputes every number from `data/curate/snpcurate.source.json`, the
   stores read through `h5_to_vcf.py`, the catalogue, the trait file and `gff/genes_index.txt`,
   independently of `tools/build_snpcurate.py`, and confirms the 7 non-site entries have no site.

16. SNPCurate marks in the other tools, also in `check_snpcurate.py`: SNPVersity on the tga1 interval with
   the 7 Z. parviglumis and 26 NAM lines (the gold badge at 46,648,374 opens SC0001; the interval line),
   SNPFunction's "Curated alleles" block for su1, Bx13 and DGAT1-2 (the Phe469 insertion listed though
   the allele list does not rank it), SNPImpact on tga1 under the default filter (N6K kept), SNPGeo's
   table and label. Test-only entries added to the loaded copy, never the source: su1 G627W (HET_ONLY)
   stays in SNPImpact's default view with its "Het only" pill and outline badge and is listed by
   SNPFunction; a coding variant of Zm00001eb406050 gets SNPFold's ring and table badge. Then the
   same pages without `js/snpcurate.data.js`: no badge, no block, no error.

Outputs in `../out/`: `results.json`, `results_brief.json`, map PNG/SVG files (the app's own
`geoBuildExportSVG()` export, rasterised with resvg) and static HTML snapshots of each page.
`check_snptrait.py`, `check_snpgeo_counts.py`, `check_gwas_handoff.py`, `check_gene_consequences.py`, `check_site_qc.py`, `check_site_qc_ui.py` and `check_snpcurate.py` recompute the numbers independently from
the catalogue, trait side-file, region records and fixture VCF.

    cd SNPTools/localdev/harness
    npm install
    PHP_BIN=$(which php) PYTHON_PATH=/path/to/python-with-h5py npm test
