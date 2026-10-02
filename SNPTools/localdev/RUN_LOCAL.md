# Running the maize SNPTools instance locally (SNPTrait + SNPGeo port)

Tree: `~/Documents/code/SNPTools_maize_port` (branch content `maize-snptrait-snpgeo` =
`main` @ `4bf370d8` + the 21 patches in `patches/`). Dataset with metadata and a local
test store: **MaizeGDB GRIN-linked 2026 · Imputed (Grzybowski 2023 sites)** (`zmgrin2026_imp`,
933 lines, release v1.4 = 926 imputed + 7 separately called, incl. NAM founder CML103). The server is a PHP built-in web server bound to 127.0.0.1 only; no admin rights needed.

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
   *Send N lines to SNPVersity* hands the selection over.
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
     - *MaizeGDB 2026 · High Quality* / *High Coverage*: "24 of 26 ... → 81 samples (all runs);
       B73 = 4 samples" and "Not available in this set: CML52, NC358";
     - *MaizeGDB GRIN-linked 2026 · Imputed (Grzybowski 2023 sites) (zmgrin2026_imp)*:
       "26 of 26 ... → 26 samples (ZmG_*); B73 = 1 sample", with no "Not available" line (release
       v1.4 adds ZmG_CML103; before v1.4 this set reported CML103 as not available).
     Choose the GRIN-linked set -> **Send**.
   - SNPVersity opens on the GRIN-linked set with chr2:4,491,424-4,499,434 and the 26 ZmG_* NAM
     lines (banner: "26 of 26 NAM lines incl. B73 (26 samples ...)"). **Build VCF & view** -> 249 variants.
   - **Send to SNPGeo** (re-queries all 933 lines) -> click the row `4494625 ...` (a GWAS SNP with
     p = 6.6e-9 that is also a release site) to map it.
   With *MaizeGDB 2026 · High Quality* SNPVersity receives 81 accessions = all runs of the 24 NAM
   lines in that catalogue, including the 4 B73 reference runs (there is no local HDF5 for it, so
   the query itself fails locally). Note: SNPVersity's own default selection on page load is 23
   accessions (one run per tagged NAM founder, no B73); that is what remains if the accession
   checkbox is unticked.

Anything outside the fixture windows returns "No variants"; chromosomes 3, 4 and 6 have no
store ("HDF5 file not found"). Score/consequence columns read "—" because the release INFO
does not yet carry the annotation keys. The two `mgdb2026_*` datasets have no HDF5 locally.

## 5. Headless checks (no browser)

    cd ~/Documents/code/SNPTools_maize_port/SNPTools
    make test PHP_BIN=$PHP_BIN PYTHON_PATH=$PYTHON_PATH    # npm install + jsdom scenarios + independent checks
    ls localdev/out/                                       # results.json, map PNG/SVG, HTML snapshots

## 6. Turning the tree into the git branch

The sandbox that built this could not create `.git` directories, so the folder is a plain tree.
To get the branch with its 21 commits:

    git clone https://github.com/andorfc/SNPTools.git SNPTools_git && cd SNPTools_git
    git checkout -b maize-snptrait-snpgeo 4bf370d892615783aa52bcec979ad6340b30e4da
    git am --whitespace=nowarn /path/to/patches/*.patch      # author: "SNPTools maize port"
    # optional: git rebase -r --exec 'git commit --amend --no-edit --reset-author' 4bf370d8

`diff -r SNPTools_git/SNPTools ~/Documents/code/SNPTools_maize_port/SNPTools` should then
differ only in ignored local files (hdf5/version3/*.h5, vcf/*, localdev/out, node_modules).
