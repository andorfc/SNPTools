# Running the maize SNPTools instance locally (SNPTrait + SNPGeo port)

Tree: `~/Documents/code/SNPTools_maize_port` (branch content `maize-snptrait-snpgeo` =
`main` @ `4bf370d8` + the 25 patches in `patches/`). Dataset with metadata and a local
test store: **MaizeGDB GRIN-linked 2026 · Imputed (Grzybowski 2023 sites)** (`zmgrin2026_imp`,
933 lines, release v1.4 = 926 imputed + 7 separately called, incl. NAM founder CML103) -- since
patch 0029 the only dataset offered (initial release; the MaizeGDB 2026 sets are commented out in
`js/data.js` with everything that serves them left in place). The server is a PHP built-in web server bound to 127.0.0.1 only; no admin rights needed.

## 1. One-time environment (no admin rights)

Either reuse the environment built for this port:

    ENV=~/.claude-science/conda/envs/snptools-local
    export PHP_BIN=$ENV/bin/php PYTHON_PATH=$ENV/bin/python PATH=$ENV/bin:$PATH

or create your own with any conda/mamba:

    conda create -y -n snptools-local -c conda-forge python=3.12 php nodejs h5py numpy
    conda activate snptools-local
    export PHP_BIN=$(which php) PYTHON_PATH=$(which python)

Versions used here: PHP 8.5.9 (CLI built-in server), Python 3.12 + h5py + numpy, Node 24 (tests only).

## 2. Test store (already built; rebuild if `hdf5/version3/` is empty)

    cd ~/Documents/code/SNPTools_maize_port/SNPTools
    make local-store PYTHON_PATH=$PYTHON_PATH      # 7 files hdf5/version3/zmgrin2026_chr{1,2,5,7,8,9,10}_impute.h5

The store is built from `localdev/fixtures/zmgrin2026_v1.4_chr*_testregions.vcf.gz` (3,495 sites x
933 samples; ZmG_CML103 is the last column). A store built before release v1.4 has 932 columns and
does not match the v1.4 catalogue: rerun `make local-store` once (a running server picks up the new
files; no restart needed).

**Full chromosome 2 (pilot of the per-chromosome build).** `hdf5/version3/zmgrin2026_chr2_impute.h5`
is the full Ceres build (5,179,690 sites x 933 samples, 2.0 GB; INFO with SnpEff fields, MAF, MAXR2
(4,665,215 sites), PlantCAD1/PlantCAD2 (4,865,809 SNPs) and Evo2 (665,539) rounded to 0.1, ESM1/2/3
and ESM-C (50,780 missense sites)), copied on 2026-09-30 from
`/90daydata/maizegdb/carson/grz2023/snptools_build/chr2/` (md5 9b8a52df8d08df7a875f3c84f4bf7d3b, as in
its `logs/MD5SUMS.txt`); the 3,495-site
demo cut is kept as `zmgrin2026_chr2_impute.demo.h5`. `make local-store` skips any store larger than
50 MB (`FORCE=1` replaces it with the demo cut). With the full chr2 store every chr2 region works,
not only the fixture windows. Measured through processForm.php -> h5_to_vcf.py (PHP CLI, this Mac):
gene Zm00001eb067740 (127 variants) 0.26 s for 26 NAM lines / 0.49 s for all 933; chr2:4-5 Mb
(16,583 variants) 0.43 s / 2.6 s including the browser-side parse (0.33 s / 0.85 s in PHP + Python);
the first query after installing the file is slower (0.41 s / 0.71 s for the gene) while the OS
cache warms.

**Genome-wide IBS for SNPCompare / SNPTree (GRIN-linked set).** Copy the Ceres matrices to
`distance/zmgrin2026/` (similarity.csv, missing_pct.csv, similarity_snp.csv, missing_pct_snp.csv,
ids.txt, tree_nj.nwk, tree_upgma.nwk; allele_distance.csv / co_called_sites.csv may sit there too).
They are data, not part of the patches. MaizeGDB 2026 keeps its files directly in `distance/`.

## 3. Start / stop

    cd ~/Documents/code/SNPTools_maize_port/SNPTools
    make start PHP_BIN=$PHP_BIN PYTHON_PATH=$PYTHON_PATH      # = localdev/start.sh
    open http://127.0.0.1:8877/index.html
    ...
    make stop                                                 # = localdev/stop.sh

`PORT=8899 make start` changes the port (8765 is used by Claude Science). Log: `localdev/server.log`; PID: `localdev/server.pid`.
If `stop` reports no PID but the port is busy: `lsof -iTCP:8877 -sTCP:LISTEN` and `kill <pid>`.
The page loads d3 v5 and fonts from public CDNs, so the browser needs internet access.

## 4. What to try

1. **SNPTrait** (sidebar, "Visualization & Search"): click the *MaizeGDB GRIN-linked 2026* pill.
   Lines are grouped by panel (NAM 26, Ames 282 254, WiDiv 623, Other GRIN 30). Try
   Subpopulation = SS, In Ames 282 = yes, Kernel type = Dent (14 lines), then search `iowa` (6).
   Under *GRIN trait ranges* pick *1000 Kernel Weight*, min 250, max 300, Add (299 lines).
   Within a filter section, ticked values combine with OR (SS and NSS: 182 + 117 = 299 lines);
   sections combine with AND. Each section is counted under the other sections' filters, so its
   other values stay listed with +n (patch 0035); before, ticking one value hid the rest.
   The same selector opens from SNPVersity (prototype, patch 0032): step 3 shows the selection
   (count, panel bar, subpopulation mix, chips), quick picks (NAM founders + B73 26, Ames 282 280,
   WiDiv 831, Other GRIN 30 -- panel membership, so a line can be in several -- random 2-25 %, All),
   each with Undo, and the paste/upload box (now also GRIN accessions such as `PI 550473` and line
   names). *Browse & filter lines…* opens SNPTrait as a drawer on a draft of the selection; the
   footer shows +added / -removed against SNPVersity, *Apply to SNPVersity* writes it back (Undo
   offered), Cancel / Escape / the backdrop drop it. The drawer keeps its own filters between
   openings and leaves the SNPTrait page's selection and filters alone.
   The top bar's **selection chip** (list icon + count, patch 0033) shows SNPVersity's selection in
   every tool and opens the same drawer; outside SNPVersity, Apply sets the selection SNPVersity
   queries next and a toast offers Undo / Open SNPVersity. On a 375 px screen the top bar now fits
   (16 px gutters; "SNPTools ›" and the button labels drop to icons below 560 px).
   *Send N lines to SNPVersity…* opens a dialog like GWAS Explorer's: *replace* the accessions
   selected in SNPVersity (pre-ticked) or *keep them and add* (it says how many lines are new).
   SNPTrait starts with nothing selected and keeps its own selection between visits.
2. **SNPGeo** ("Explore & Analyze"): pick the GRIN-linked dataset card; the gene box is
   pre-filled with `Zm00001eb374090` (chr9:12,838,008-12,843,999; 164 variants x 933 lines).
   Switch *North America (states / provinces)* / *World*, and the colour modes
   (carrier fraction, carriers among called, alternative allele frequency). Click a country
   for the state table and carrier list; arrow keys step through variants; PNG/SVG export.
   Other genes in the test store: Zm00001eb374230, Zm00001eb067740, Zm00001eb056510,
   Zm00001eb233650, Zm00001eb313510, Zm00001eb404740, Zm00001eb404760.
3. **SNPVersity -> SNPGeo**: choose the GRIN-linked dataset, region chr10:9,788,000-9,826,500,
   the default NAM selection, *Run*; then *Send to SNPGeo*. SNPGeo re-queries the region for
   all 933 lines (1,297 variants).

4. **GWAS Explorer -> SNPVersity with the GRIN-linked set** (NAM GWAS, Tibbs-Cortes et al. 2024):
   - *GWAS Explorer* (sidebar) -> trait dropdown -> **Tassel Branch Number** -> chip
     **Tassel Branch Number - Trait** (NAM · Tibbs-Cortes et al. 2024).
   - Type `chr2:4491424..4499434` in *Jump to SNP or region…* and press Enter; click
     **Select region** and drag across the visible window (the view is padded by ~8 %, so the
     panel shows roughly Chr2:4.490–4.500 Mb; the headless test selects exactly
     4,491,424–4,499,434: 85 SNPs, 11 above the threshold, lead S2_4498985, p = 3.6e-14).
   - **Send to SNPVersity →**. *Send this genomic region* and **Send accessions — replace ...** are
     pre-ticked. Each entry of **VCF set in SNPVersity** shows how many of the 26 NAM lines
     (B73 + 25 founders) it contains, and the text under the checkboxes names the missing ones:
     - (while MaizeGDB 2026 was offered: *High Quality* / *High Coverage*: "24 of 26 ... → 81
       samples (all runs); B73 = 4 samples" and "Not available in this set: CML52, NC358";)
     - *MaizeGDB GRIN-linked 2026 · Imputed (Grzybowski 2023 sites) (zmgrin2026_imp)*:
       "26 of 26 ... → 26 samples (ZmG_*); B73 = 1 sample", with no "Not available" line (release
       v1.4 adds ZmG_CML103; before v1.4 this set reported CML103 as not available).
     Choose the GRIN-linked set -> **Send**.
   - SNPVersity opens on the GRIN-linked set with chr2:4,491,424-4,499,434 and the 26 ZmG_* NAM
     lines (banner: "26 of 26 NAM lines incl. B73 (26 samples ...)"). **Build VCF & view** -> 249 variants.
   - **Send to SNPGeo** (re-queries all 933 lines) -> click the row `4494625 ...` (a GWAS SNP with
     p = 6.6e-9 that is also a release site) to map it.
   (Before patch 0029:) with *MaizeGDB 2026 · High Quality* SNPVersity received 81 accessions = all runs of the 24 NAM
   lines in that catalogue, including the 4 B73 reference runs (there is no local HDF5 for it, so
   the query itself fails locally). Note: SNPVersity's own default selection on page load is 23
   accessions (one run per tagged NAM founder, no B73); that is what remains if the accession
   checkbox is unticked.

Anything outside the fixture windows returns "No variants"; chromosomes 3, 4 and 6 have no
store ("HDF5 file not found"). The two `mgdb2026_*` datasets (commented out) have no HDF5 locally.

**Sites with several consequences (patch 0029).** SnpEff writes one consequence per gene: 1,464 of
the 3,744 GRIN-linked fixture sites list more than one (`GENEMODEL=Zm00001eb374100,Zm00001eb374090;
TYPE=downstream_gene_variant,intron_variant`). SNPVersity and SNPImpact show the first (most
severe) entry, and SNPVersity's Gene model column marks the rest with "+N" (hover). SNPFunction,
SNPFold and SNPGeo's gene search read the gene's own entry (`Data.rowForGene`); they had kept only
sites whose FIRST gene was theirs: Zm00001eb374230 kept 24 of its 63 sites (11 of 25 missense),
Zm00001eb374090 99 of 164, Zm00001eb404750 2 of 18. Try SNPFunction on `Zm00001eb374230`. A site
carries ONE ESM score, computed for the first missense entry of TYPE (`tools/annotate_release_info.py`),
so a gene whose substitution is not that entry shows no ESM there (13 of Zm00001eb374230's 25
missense sites) rather than its neighbour's score.

**SNPVersity annotation columns (GRIN-linked set).** 11 columns (Gene model ... ESM3): MQ and COMP
are left out for this set (patch 0030; status 'hidden' in `Data.annotationFields`), because the
Grzybowski et al. 2023 call set records no per-site MQ or coverage. Filled: Gene model, Effect, SNPEff
Impact (SnpEff 5.2a fields of the MaizeGDB Schnable scored VCFs), MAF (from the 933 release
genotypes), ESM1/ESM2/ESM3 (missense sites only; ESM2 = store ESM-2 650M, `esm2_store_score`).
Evo2 (INFO `evo2_score`: SNPs within 1 kb of a gene) and ESM-C (`ESMC_score`: missense) are
columns after PlantCAD2 and ESM3 (patch 0033); the full chr2 store fills them (chr2 GWAS window:
120 and 10 of 249 sites), the demo stores of the other chromosomes predate them, so they read
pending there. SNPImpact, SNPFunction (catalog, burden means, CSV), SNPFold and SNPGeo show them as
columns only where the region or gene has a score (Zm00001eb067740: Evo2 on 92 sites, ESM-C on 10).
Domain is the Pfam block covering the site in `data/domains/by_chr/<chr>.json` (canonical proteins'
domains mapped to the genome); "—" is the usual answer, since the blocks cover 0.8-1.0% of each
chromosome (46,930 of chr2's 5,179,690 sites fall in one). On the fixture sites the lookup agrees
with the protein-coordinate domains (`domains.by_gene.json`) on 232 of 237 coding sites; the other 5
lie 1-4 residues from a domain edge (residue numbering of another isoform, or the edge codon). Empty with a dashed header and a note above the
table: PlantCAD1/PlantCAD2 where a store has no scores yet.
maxR² is available for this set (highest PLINK 1.9 r² with any variant 400-5,000 bp away, from the
933 release genotypes, no filtering; blank = no partner variant in range or monomorphic); the
full chr2 store carries it (chr2 GWAS window: 233 of 249 sites); the demo stores of the other
chromosomes predate it, so their cells read NA until the rebuilt stores arrive. PlantCAD1/PlantCAD2
and Evo2 are rounded to 0.1 like the ESM scores. With the full chr2 store PlantCAD1/PlantCAD2 are
filled on chr2 (SNPs; indels read N/A); the other chromosomes still use the demo cut and show them
as pending.

5. **SNPCompare genome-wide (GRIN-linked set)**: open *SNPCompare*, choose **Dataset = MaizeGDB
   GRIN-linked 2026** (or arrive from SNPVersity with that set), focal `ZmG_B73`, scope
   **Genome-wide** -> **Table**: 933 lines; nearest ZmG_DJ7 (0.99823), ZmG_F42 (0.99575).
   *Genome-wide sites* switches between all 46,054,265 sites and SNPs only. `ZmG_MO17` -> nearest
   ZmG_SEAGULLSEVENTEEN (0.95391). Without `distance/zmgrin2026/` the scope stays disabled with a note.
6. **SNPTree**: with the GRIN-linked set, a *Genome-wide tree* card offers the precomputed
   Neighbour-Joining and UPGMA trees (933 tips) as Newick downloads; region trees are unchanged.
7. **SNPFunction Gene Ontology & pathways** (needs internet): pick the GRIN-linked dataset and
   analyze `Zm00001eb056510` (adh1). The *Gene Ontology & pathways* card is MaizeGDB's own view,
   read live from `https://claude.maizegdb.org/api/v1/records/gene/<id>?fields=function` (CORS
   open, ~0.4 s): 10 GO terms, 11 of 94 plant-slim categories lit, the ancestry graph per aspect;
   *Pathways* shows 6 E2P2 pathways (adh1's step in gold) and 4 KEGG maps. Set
   `window.SNPTOOLS_MAIZEGDB_BASE` before `js/snpfunction-ontology.js` loads to read another
   MaizeGDB host. Offline, or if MaizeGDB fails, the card says so and the old GO/KEGG lists from
   `data/function/annotations/` are shown instead.

## 5. Headless checks (no browser)

    cd ~/Documents/code/SNPTools_maize_port/SNPTools
    make test PHP_BIN=$PHP_BIN PYTHON_PATH=$PYTHON_PATH    # npm install + jsdom scenarios + independent checks
    ls localdev/out/                                       # results.json, map PNG/SVG, HTML snapshots

## 6. Turning the tree into the git branch

The sandbox that built this could not create `.git` directories, so the folder is a plain tree.
To get the branch with its 25 commits:

    git clone https://github.com/andorfc/SNPTools.git SNPTools_git && cd SNPTools_git
    git checkout -b maize-snptrait-snpgeo 4bf370d892615783aa52bcec979ad6340b30e4da
    git am --whitespace=nowarn /path/to/patches/*.patch      # author: "SNPTools maize port"
    # optional: git rebase -r --exec 'git commit --amend --no-edit --reset-author' 4bf370d8

`diff -r SNPTools_git/SNPTools ~/Documents/code/SNPTools_maize_port/SNPTools` should then
differ only in ignored local files (hdf5/version3/*.h5, vcf/*, localdev/out, node_modules).
