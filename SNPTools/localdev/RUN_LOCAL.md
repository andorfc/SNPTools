# Running the maize SNPTools instance locally (SNPTrait + SNPGeo port)

Tree: `~/Documents/code/SNPTools_maize_port` (branch content `maize-snptrait-snpgeo` =
`main` @ `4bf370d8` + the 17 patches in `patches/`). Dataset with metadata and a local
test store: **MaizeGDB GRIN-linked 2026 · Imputed (Grzybowski 2023 sites)** (`zmgrin2026_imp`,
932 lines). The server is a PHP built-in web server bound to 127.0.0.1 only; no admin rights needed.

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

## 3. Start / stop

    cd ~/Documents/code/SNPTools_maize_port/SNPTools
    make start PHP_BIN=$PHP_BIN PYTHON_PATH=$PYTHON_PATH      # = localdev/start.sh
    open http://127.0.0.1:8765/index.html
    ...
    make stop                                                 # = localdev/stop.sh

`PORT=8877 make start` changes the port. Log: `localdev/server.log`; PID: `localdev/server.pid`.
If `stop` reports no PID but the port is busy: `lsof -iTCP:8765 -sTCP:LISTEN` and `kill <pid>`.
The page loads d3 v5 and fonts from public CDNs, so the browser needs internet access.

## 4. What to try

1. **SNPTrait** (sidebar, "Visualization & Search"): click the *MaizeGDB GRIN-linked 2026* pill.
   Lines are grouped by panel (NAM 25, Ames 282 254, WiDiv 623, Other GRIN 30). Try
   Subpopulation = SS, In Ames 282 = yes, Kernel type = Dent (14 lines), then search `iowa` (6).
   Under *GRIN trait ranges* pick *1000 Kernel Weight*, min 250, max 300, Add (299 lines).
   *Send N lines to SNPVersity* hands the selection over.
2. **SNPGeo** ("Explore & Analyze"): pick the GRIN-linked dataset card; the gene box is
   pre-filled with `Zm00001eb374090` (chr9:12,838,008-12,843,999; 164 variants x 932 lines).
   Switch *North America (states / provinces)* / *World*, and the colour modes
   (carrier fraction, carriers among called, alternative allele frequency). Click a country
   for the state table and carrier list; arrow keys step through variants; PNG/SVG export.
   Other genes in the test store: Zm00001eb374230, Zm00001eb067740, Zm00001eb056510,
   Zm00001eb233650, Zm00001eb313510, Zm00001eb404740, Zm00001eb404760.
3. **SNPVersity -> SNPGeo**: choose the GRIN-linked dataset, region chr10:9,788,000-9,826,500,
   the default NAM selection, *Run*; then *Send to SNPGeo*. SNPGeo re-queries the region for
   all 932 lines (1,297 variants).

Anything outside the fixture windows returns "No variants"; chromosomes 3, 4 and 6 have no
store ("HDF5 file not found"). Score/consequence columns read "—" because the release INFO
does not yet carry the annotation keys. The two `mgdb2026_*` datasets have no HDF5 locally.

## 5. Headless checks (no browser)

    cd ~/Documents/code/SNPTools_maize_port/SNPTools
    make test PHP_BIN=$PHP_BIN PYTHON_PATH=$PYTHON_PATH    # npm install + jsdom scenarios + independent checks
    ls localdev/out/                                       # results.json, map PNG/SVG, HTML snapshots

## 6. Turning the tree into the git branch

The sandbox that built this could not create `.git` directories, so the folder is a plain tree.
To get the branch with its 17 commits:

    git clone https://github.com/andorfc/SNPTools.git SNPTools_git && cd SNPTools_git
    git checkout -b maize-snptrait-snpgeo 4bf370d892615783aa52bcec979ad6340b30e4da
    git am --whitespace=nowarn /path/to/patches/*.patch      # author: "SNPTools maize port"
    # optional: git rebase -r --exec 'git commit --amend --no-edit --reset-author' 4bf370d8

`diff -r SNPTools_git/SNPTools ~/Documents/code/SNPTools_maize_port/SNPTools` should then
differ only in ignored local files (hdf5/version3/*.h5, vcf/*, localdev/out, node_modules).
