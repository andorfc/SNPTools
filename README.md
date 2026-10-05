# SNPTools

> A browser-based suite of tools for exploring maize genomic variation, developed by
> [MaizeGDB](https://maizegdb.org) and funded by USDA-ARS.

SNPTools is a single-page web application for querying a maize variant dataset across hundreds of
lines, annotating each variant with its predicted effect and DNA / protein language-model scores,
and handing a selection to a family of connected tools: trait-based line selection, geography,
variant prioritization, gene function, protein structure, similarity and phylogeny. All
coordinates are on the **B73 v5** reference assembly (Zm-B73-REFERENCE-NAM-5.0).

This release offers one dataset, **MaizeGDB GRIN-linked 2026** (933 lines, 46,054,265 sites); the
earlier datasets are kept in the code, switched off (see [Datasets](#datasets)).

---

## Table of contents

- [How it works](#how-it-works)
- [Repository layout](#repository-layout)
- [The tools](#the-tools)
- [Datasets](#datasets)
- [Backend services](#backend-services)
- [Data files](#data-files)
- [Data builders](#data-builders)
- [Getting started](#getting-started)
- [Configuration](#configuration)
- [Testing](#testing)
- [Data provenance & citation](#data-provenance--citation)
- [Acknowledgments](#acknowledgments)
- [License](#license)

---

## How it works

The frontend is plain JavaScript with no build step: a small shell plus self-registering tool
modules, loaded in order by `index.html`.

- **`js/core.js`** loads first: icons, the tool registry, the shared state `S`, the router
  (`go(id)`), the left-rail navigation, tooltips, the top-bar selection chip and toasts. Each tool
  registers itself with `SNPTools.register(id, { render })`.
- **`js/data.js`** is the data layer. Everything that touches the variant store, the accession
  catalogue, gene models, domains or structures goes through the `Data` object.
- **Tool modules** (`js/snpversity.js`, `js/snptrait.js`, ...) each own one page and its state.
  A result handed from one tool to another is stored on `S` (`S.treeInput`, `S.impactInput`,
  `S.geoInput`, ...), which is how "Send selection to..." works. `js/snphandoff.js` holds the
  shared hand-off controls (replace or add to SNPVersity's selection).

```
Browser (core.js + tool modules)
   │ region + accessions                      │ gene id                    │ focal line
   ▼                                          ▼                            ▼
processForm.php ─► h5_to_vcf.py           lookupGeneModel.php         ibsCompare.php
   │                 │                        │ gff/genes_index.txt        │ distance/<family>/
   │                 ▼                        ▼                            ▼
   │   hdf5/grin2026/zmgrin2026_<chr>_impute.h5    {chromosome, start, end}   one IBS row
   ▼
vcf/snpv_*.vcf.gz ─► data.js parses it ─► rows ─► SNPVersity table and the other tools
                                   ├─ data/domains/by_chr/<chr>.json      Pfam domain by position
                                   ├─ data/genemodels/by_chr/<chr>.json   exon / CDS structure
                                   ├─ data/structures/<model>/structure-<gene>.js   (SNPFold)
                                   └─ MaizeGDB record API  /api/v1/records/gene/<id>?fields=function
                                      (SNPFunction's GO and pathway views, read live)
```

---

## Repository layout

```
.
├── README.md
└── SNPTools/                           # the web root
    ├── index.html                      # app shell + MaizeGDB masthead; loads every script
    ├── css/                            # main.css, paneffect-native.css, snpfunction-ontology.css
    ├── images/                         # MaizeGDB logos
    ├── js/
    │   ├── core.js  data.js  snphandoff.js  snphelp.js
    │   ├── snpversity.js  snpgwas.js  snptrait.js                    # Visualization & Search
    │   ├── snpimpact.js  snpfunction.js  snpfunction-ontology.js
    │   ├── snpfold.js  snpgeo.js  snpgeo.regions.js                  # Explore & Analyze
    │   ├── snpcompare.js  snptree.js  snpmatrix.js                   # Compare & Relate
    │   ├── snppaneffect.js  paneffect-*.js  pe/                      # PanEffect (native port)
    │   ├── zmgrin.catalog.js           # GRIN-linked lines + SNPTrait schema (generated)
    │   ├── zmgrin.lineqc.js            # per-line heterozygosity (generated, tools/build_site_qc.py)
    │   ├── snpcurate.js  snpcurate.data.js   # SNPCurate page + its registry (generated, tools/build_snpcurate.py)
    │   ├── accessions.catalog.js       # catalogue of the earlier datasets (generated)
    │   ├── accessions.real.js          # legacy flat accession list (fallback)
    │   └── mgdb-genome-groups.js       # NAM founder heterotic-group colours (from MaizeGDB)
    ├── processForm.php                 # region + accessions -> VCF (runs h5_to_vcf.py)
    ├── h5_to_vcf.py                    # HDF5 store -> VCF extractor
    ├── lookupGeneModel.php             # gene model -> coordinates
    ├── gene_index_lib.php              # sorted gene index used by lookupGeneModel.php
    ├── ibsCompare.php                  # precomputed genome-wide IBS rows and trees
    ├── gff/                            # genes_data.serialized (source) + genes_index.txt
    ├── hdf5/grin2026/                  # variant stores (not in git)
    ├── distance/<family>/              # genome-wide IBS matrices and trees (not in git)
    ├── data/                           # static annotation files (only data/gwas is in git)
    ├── paneffect/                      # PanEffect score files (not in git)
    ├── vcf/                            # VCFs written by processForm.php (not in git)
    ├── tools/                          # data builders (see below) + build_gene_index.php
    ├── build_catalog.py                # builds accessions.catalog.js (earlier datasets)
    ├── verify_ibs_ids.py               # checks distance/ids.txt against that catalogue
    ├── Makefile                        # builders, local server, tests (make help)
    └── localdev/                       # local instance, fixtures and headless test harness
```

---

## The tools

| Tool | What it does |
| --- | --- |
| **SNPVersity** | The front door. Choose a dataset, an interval or a B73 v5 gene model, and the lines (quick picks, a pasted list of IDs / line names / GRIN accessions, or *Browse & filter*, which opens SNPTrait over the page). Returns a colour-coded genotype table and a downloadable VCF with predicted effects, Pfam domain, MAF, maxR², Site QC, and PlantCAD1/2, Evo2, ESM1/2/3 and ESM-C scores. A site with consequences in several genes shows its most severe one and "+N" for the others. Every site is shown with its Site QC class; a filter hides flagged and no-carrier sites or keeps passing ones only. |
| **GWAS Explorer** | Manhattan plots of published NAM GWAS results (`data/gwas/`); drag-select a peak and send the region and/or the lines to SNPVersity. |
| **SNPTrait** | Select lines by GRIN passport and evaluation data: panel, subpopulation, origin, improvement status, kernel type, binned and numeric trait ranges, and sample heterozygosity (Inbred / Elevated heterozygosity / Heterozygous sample, from `js/zmgrin.lineqc.js`). Several values in one filter section are ORed, sections are ANDed. Send the lines to SNPVersity (replace or add). |
| **SNPImpact** | Ranks a region's variants by a combined PlantCAD + ESM score with the predicted consequence and Pfam domain; filters and a shortlist. Flagged and no-carrier sites are hidden by default (Site QC filter; the note above the table says how many). Opened empty, it offers SNPVersity's current result ("Load data from SNPVersity"), as SNPTree and SNPCompare do. |
| **SNPFunction** | A gene dossier: domains, MaizeGDB's Gene Ontology and pathway views (read live from the MaizeGDB record API; the local annotation file is the fallback), and the gene's variant burden and damaging alleles across the whole panel, with their carriers. It classes every site from its own carrier counts: the allele list, knockout lines and burden use usable sites, flagged calls and alleles with no carrier sit in two collapsed groups, and a banner warns when most protein-changing sites of the gene are flagged. An allele's carrier panel opens its homozygous carriers, heterozygous carriers, all carriers or the reference lines (homozygous for the reference allele) in SNPVersity; SNPFold's carrier panel does the same. |
| **SNPCurate** | A registry of published maize alleles (21 in release v1.4: validated causal changes, published markers and tagging sites), each with its carriers by subpopulation and country, Site QC class, scores, priority, rank in its gene and GRIN trait means in this release, or the reason the release cannot show it. Gold, outline and grey marks; hand-offs of the carriers to SNPVersity, the gene to SNPFunction and the site to SNPGeo. The same marks, opening the record, appear on curated alleles in SNPVersity (beside the position, and a line naming the curated alleles of the interval), SNPImpact (never hidden by the Site QC or score filters), SNPFunction (a "Curated alleles" block for the gene, grey ones included), SNPFold (a ring on the track and the table) and SNPGeo. |
| **SNPFold** | Coding variants on the predicted protein structure (AlphaFold2, Boltz2 or ESMFold model): a linear protein browser (domains, secondary structure, pLDDT, InterProScan sites, disorder) and a 3D view. "Usable alleles only" (on by default) leaves flagged and no-carrier variants out of the track, the table and the 3D view. |
| **SNPGeo** | Where each allele is found: countries, U.S. states, Canadian provinces and Mexican states coloured by reference/alternative composition, carrier fraction or allele frequency, from the GRIN origin of each line. Each variant shows its Site QC class. |
| **SNPCompare** | Lines ranked by identity-by-state to a focal line: genome-wide (precomputed), in the region, or both with the difference. Opens SNPMatrix for the region's pairwise matrix. The table lists each line's PI number (GRIN accession). The region scope uses usable sites by default ("Region sites"); the genome-wide matrices use all sites. |
| **SNPTree** | Neighbour-joining / UPGMA trees from the region's genotypes (Newick, MEGA, PHYLIP); the precomputed genome-wide trees as downloads. Usable sites by default, with a switch to all sites (SNPMatrix likewise). |
| **PanEffect** | Every amino-acid substitution of a protein scored by ESM, for B73 and across the pan-genome, with domains and secondary structure (native port of MaizeGDB's PanEffect). |

SNPImpute, SNPDensity and SNPGermplasm are roadmap placeholders and are not in the navigation.

---

## Datasets

Datasets are defined in `js/data.js` (`DATASETS`); `processForm.php` maps each id to its store.

| Dataset id | Name | Lines | Sites | Notes |
| --- | --- | --- | --- | --- |
| `zmgrin2026_imp` | MaizeGDB GRIN-linked 2026 (release v1.4) | 933 | 46,054,265 | Grzybowski et al. (2023) sites, Beagle 5 imputation; 926 imputed lines + 7 lines called separately at the same sites (including NAM founder CML103); every line linked to its USDA GRIN accession |

Per site the stores carry SnpEff 5.2a effects (TYPE, EFFECT, GENEMODEL, SUB), MAF and maxR² from
the 933 genotypes, PlantCAD1, PlantCAD2 and Evo2 (every SNP: 43,296,332 sites), ESM1b, ESM2, ESM3
and ESM-C (missense sites) and the site-QC fields NHET, NHOM and SITEQC (every site), on all ten
chromosomes. Mapping quality and coverage are not recorded for this call set, so SNPVersity leaves
those columns out (`Data.annotationFields`).

The stores are the v1.4.2 annotation build of the release (2026-10-04): the lines, sites and
genotypes of v1.4, unchanged, with Evo2 extended from the SNPs within 1 kb of a gene to every SNP
and the three site-QC fields written into each row's INFO. The app and the VCFs it writes keep the
label "release v1.4".

**Site QC.** The lines are inbreds, yet 27.5% of non-reference calls are heterozygous, mostly at
sites where heterozygous carriers outnumber homozygous ones (reads from another copy of the sequence
mapping there), and 10.9% of sites have no carrier among the 933 lines (the site list was built on a
larger panel). Each site is classed from its heterozygous (NHET) and homozygous (NHOM) carriers among
all 933 lines, first match wins: `NO_CARRIER`, `HET_ONLY` (no homozygous carrier), `HET_EXCESS`
(NHET > NHOM), `HET_ELEVATED` (heterozygous carriers a quarter or more of all carriers), `PASS`. Of
46,054,265 sites: 18,177,811 PASS, 3,756,376 HET_ELEVATED, 14,354,590 HET_EXCESS, 4,752,396 HET_ONLY,
5,013,092 NO_CARRIER. *Usable* means PASS or HET_ELEVATED; *flagged* means HET_ONLY or HET_EXCESS.
Twenty lines are heterozygous as sequenced (more than 5% of the 16,463,157 clean sites, PASS with at
least 3 homozygous carriers; the median line is at 0.78%): SNPTrait filters them (*Sample
heterozygosity*), and SNPVersity and SNPFunction mark them "het".
Flagged sites are kept in the data and in every download; each view's default
(`Data.SITE_QC_DEFAULTS`) is: SNPVersity and SNPGeo show every site with its class; SNPImpact,
SNPFunction and SNPFold hide flagged and no-carrier sites; SNPTree, SNPMatrix and SNPCompare (region
scope) use usable sites. The pages read the counts from each VCF's INFO: the v1.4.2 stores carry
them, and for a store that does not `h5_to_vcf.py` adds them from the sidecars
`tools/build_site_qc.py` writes. Without them (`SNPTOOLS_SITEQC=0`, or such a store without its
sidecar) the views behave as before, except SNPFunction (and SNPFold, which reads SNPFunction's
full-panel query), which counts carriers itself.

The MaizeGDB 2026 (High Quality / High Coverage), MaizeGDB 2024, Schnable 2023 and NAM 2021 sets are
commented out of `DATASETS`; their backend handling (`processForm.php`, `ibsCompare.php`, the
accession catalogue) is kept, so uncommenting them offers them again once their stores are installed.

### Accession catalogue

`js/zmgrin.catalog.js` defines `window.SNP_CATALOG.families.zmgrin2026` (the 933 lines with their
GRIN accession, panel memberships, subpopulation and origin) and `window.SNPTRAIT_SCHEMA.zmgrin2026`
(SNPTrait's facets); `data/traits/zmgrin2026.traits.json` holds the GRIN trait summaries. Sample ids
(`ZmG_<genotype>`) are the column names in the stores, so a selection maps directly to HDF5 columns.

---

## Backend services

| Endpoint | Called by | Purpose |
| --- | --- | --- |
| `processForm.php` | `Data.queryVariants` | POST `chr`, `start`, `end` (whole numbers), `dataSet`, `genotypes` (JSON list of ids), `outName`. Runs `h5_to_vcf.py`, writes `vcf/snpv_*.vcf.gz` and answers `{status: success, outFile, variants}`, `{status: empty}` or `{status: error, message}`. Refuses requests over `SNPTOOLS_MAX_CELLS`, never overwrites an existing VCF, deletes its VCFs after `SNPTOOLS_VCF_TTL_HOURS`. Replies name no server path (details go to the PHP error log). |
| `lookupGeneModel.php` | `Data.lookupGene` | `?geneModelId=<id>` -> `{chromosome, start, end, ID}` by a binary search of `gff/genes_index.txt`, the sorted index of `gff/genes_data.serialized` (`gene_index_lib.php`). The index rebuilds itself when the source changes and `gff/` is writable (or run `php tools/build_gene_index.php`); without a current index the source answers directly. Unknown ids answer `{"id":"empty", ...}`. |
| `ibsCompare.php` | SNPCompare, SNPTree | `?focal=<id>&dataset=<family>[&sites=snp]` -> `{focal, rows:[{id, similarity, missing}]}` from `distance/<family>/`; `?probe=1&dataset=<family>` says what is installed; `?tree=nj\|upgma&dataset=<family>` serves a precomputed Newick. Without `dataset` it serves `zmgrin2026`. |
| MaizeGDB record API (external) | SNPFunction | `GET <base>/api/v1/records/gene/<id>?fields=function`: GO terms with evidence, plant-slim categories, E2P2 pathways and KEGG maps. CORS-open. The base is `window.SNPTOOLS_MAIZEGDB_BASE` (default `https://claude.maizegdb.org`). |

`h5_to_vcf.py <store.h5> <out.vcf[.gz]> <start> <end> <accessions.json>` reads the store in
chunk-aligned blocks (peak memory about 0.3 GB for any interval), fills an id the store lacks with
`./.`, and writes gzip when the output name ends in `.gz`. Store files are named
`hdf5/grin2026/<family>_<chr>_<tier>.h5` (`zmgrin2026_chr1_impute.h5`); each holds the datasets
`CHROM, POS, REF, ALT, QUAL, INFO` and one int8 genotype column per sample (0 = 0/0, 1 = 0/1,
2 = 1/1, 3 = missing), written by `tools/vcf_to_h5.py`. Each row's INFO ends with `NHET`, `NHOM`
(heterozygous and homozygous-alternate carriers among all 933 lines, whatever lines were asked for)
and `SITEQC` (`NO_CARRIER`, `HET_ONLY`, `HET_EXCESS`, `HET_ELEVATED` or `PASS`, derived from the two
counts). A store built from a release VCF that has the three (v1.4.2) carries them in its own rows,
which are written as stored. To the rows of a store that does not, `h5_to_vcf.py` adds them from
the store's current site-QC sidecar (`<store>.siteqc.h5`, `tools/build_site_qc.py`). A missing or
stale sidecar (site count or store size changed) puts a `Note:` line in the extractor's output, and
a store without the three in its own rows is then written without them. With `SNPTOOLS_SITEQC=0`
every VCF is written without them, whatever their source.

INFO keys the browser reads: `GENEMODEL`, `TYPE`, `EFFECT`, `SUB` (parallel comma lists, one entry
per affected gene, most severe first), `MAF`, `MAXR2`, `plantcad1_score`, `plantcad2_score`,
`evo2_score`, `ESM1_score`, `ESM2_score`, `ESM3_score`, `ESMC_score`, `NHET`, `NHOM`, `SITEQC`; and for
the earlier datasets `MQ`, `CVP`, `DNA_SCORE`, `AA_SCORE`.

---

## Data files

Static files under `SNPTools/`, loaded on demand and cached by `js/data.js` and the tools:

| Path | Used by | Content |
| --- | --- | --- |
| `hdf5/grin2026/zmgrin2026_<chr>_impute.h5` | processForm.php | the variant stores (1.4-2.8 GB per chromosome) |
| `hdf5/grin2026/zmgrin2026_<chr>_impute.siteqc.h5` | build_snpcurate.py, build_site_qc.py `--summary`; h5_to_vcf.py for a store without the fields in its INFO | site-QC sidecar per store: heterozygous, homozygous-alternate and missing calls per site over all 933 lines (`NHET`, `NHOM`, `NMISS`; 13-26 MB each; not in git) |
| `data/qc/zmgrin2026.siteqc.summary.json`, `zmgrin2026.lineqc.tsv` | reports, help text | site classes per chromosome and in total; per-line heterozygosity (generated) |
| `data/curate/snpcurate.source.json` | build_snpcurate.py | SNPCurate's registry of published alleles, edited by hand |
| `js/snpcurate.data.js` | SNPCurate (optional) | the registry with each site's carriers, class, scores, priority and trait means in this release (`window.SNP_CURATE`, generated) |
| `js/zmgrin.lineqc.js` | SNPTrait, SNPVersity, SNPFunction (optional) | per-line heterozygous share at clean sites and its class (`window.SNP_LINE_QC`, generated); `Data.accessionsFor` adds `hetShare` and `sampleQC` to each line |
| `data/domains/by_chr/<chr>.json` (fallback `domains.by_chr.json`) | SNPVersity, SNPImpact | Pfam domain blocks by genomic position |
| `data/domains/domains.by_gene.json`, `domains.by_protein.json` | SNPImpact, SNPFold, SNPFunction | each gene's canonical-protein domains |
| `data/genemodels/by_chr/<chr>.json` | SNPImpact, SNPFunction, SNPFold | canonical exon / CDS structure |
| `data/structures/{alphafold,boltz,esmfold}/structure-<gene>.js` | SNPFold | per-gene model with pLDDT, secondary structure and domains (`window.SNPFOLD_STRUCT` / `SNPFOLD_PDB`) |
| `data/results.sites.tsv`, `data/B73_iupred2a_data_*` | SNPFold | InterProScan residue sites; IUPred2A disorder |
| `data/function/annotations/<gene>.json` | SNPFunction | local functional annotation (fallback for the GO / pathway views; its GO is computational, see `data/function/README.md`) |
| `data/traits/zmgrin2026.traits.json` | SNPTrait | GRIN trait summaries per line |
| `data/geo/countries.geo.json`, `admin1_na.geo.json` | SNPGeo | Natural Earth base maps |
| `data/gwas/` | GWAS Explorer | published NAM GWAS results + `manifest.json` (in git) |
| `distance/zmgrin2026/` | SNPCompare, SNPTree | genome-wide IBS matrices (`similarity*.csv`, `missing_pct*.csv`, `ids.txt`) and trees |
| `paneffect/` | PanEffect | ESM score files per gene |
| `gff/genes_data.serialized`, `gff/genes_index.txt` | lookupGeneModel.php | gene-model coordinates (in git) |

---

## Data builders

The builders for the GRIN-linked dataset are in `SNPTools/tools/` (details in
`tools/README_builders.md`; run order in `make help`):

| Script | Output |
| --- | --- |
| `build_maize_samples.py` | `data/zmgrin2026_samples.tsv` (the 933 sample ids and their GRIN links) |
| `build_maize_trait_catalog.py` | `js/zmgrin.catalog.js`, `data/traits/zmgrin2026.traits.json` |
| `build_maize_snpgeo_data.py`, `build_geo_layers.py` | `js/snpgeo.regions.js`, `data/geo/*.geo.json` |
| `annotate_release_info.py` | the SnpEff, MAF, maxR², site QC (`NHET`, `NHOM`, `SITEQC`) and language-model INFO fields of a release VCF |
| `maxr2_chr.sh`, `check_r2.py` | the maxR² table (PLINK 1.9) and its independent check |
| `vcf_to_h5.py` | `hdf5/grin2026/zmgrin2026_<chr>_impute.h5` from an annotated release VCF |
| `build_snpcurate.py` | `js/snpcurate.data.js` from `data/curate/snpcurate.source.json` (`make curate`); `--table` writes the benchmark table |
| `build_site_qc.py` | the `.siteqc.h5` sidecar beside each store; with `--summary`, `data/qc/` and `js/zmgrin.lineqc.js` (`make site-qc`) |
| `build_gene_index.php` | `gff/genes_index.txt` |

The domain, gene-model and structure files (`data/domains/`, `data/genemodels/`,
`data/structures/`) and the PanEffect score files are produced by MaizeGDB pipelines outside this
repository.

---

## Getting started

### Requirements

- A web server that serves static files and runs PHP (Apache or nginx + PHP-FPM; locally the PHP
  built-in server). Tested with PHP 8.5.
- Python 3 with `h5py` and `numpy` for `h5_to_vcf.py` (tested with Python 3.12, h5py 3.16,
  numpy 2.5).
- The variant stores in `hdf5/grin2026/` with their site-QC sidecars (`make site-qc`), and the data
  files above.

### Local instance

`SNPTools/localdev/RUN_LOCAL.md` has the full recipe (environment, stores, what to try):

```bash
cd SNPTools
make start PHP_BIN=/path/to/php PYTHON_PATH=/path/to/python    # http://127.0.0.1:8877/index.html
make stop
```

`make local-store` builds small test stores from `localdev/fixtures/` when the full stores are not
installed (it never replaces a store larger than 50 MB unless `FORCE=1`).

### Deployment

1. Put `SNPTools/` under the web root.
2. Install the stores in `hdf5/grin2026/`, the matrices in `distance/zmgrin2026/` and the data files.
   Run `make site-qc` where the stores are (or copy their `.siteqc.h5` sidecars with them): it writes
   the sidecars (about 1.5 minutes with `JOBS=4`) and regenerates `data/qc/` and `js/zmgrin.lineqc.js`,
   which are also in git. The v1.4.2 stores carry the site-QC fields themselves, so the pages show
   Site QC with or without the sidecars; `make curate` and the summary need them, and so do the pages
   for a store without the fields. After replacing a store, run `make site-qc` and then `make curate`
   (`js/snpcurate.data.js`, also in git) again.
3. Make `vcf/` writable by the web server (and `gff/`, if the gene index should rebuild itself).
4. Set the environment variables below for the PHP process, and point SNPFunction at a MaizeGDB
   host that serves `/api/v1` if not the default.

---

## Configuration

Browser side, in `js/data.js` (`CFG`):

| Key | Default | Purpose |
| --- | --- | --- |
| `endpoint` | `processForm.php` | variant query endpoint |
| `vcfDir` | `vcf/` | web path of the VCFs (must match `processForm.php`) |
| `geneEndpoint` | `lookupGeneModel.php` | gene model -> coordinates |
| `structDir` | `data/structures/` | SNPFold structure files |
| `domainsDir` / `domainsUrl` / `domainsGeneUrl` | `data/domains/...` | Pfam domain files |
| `geneModelsDir` | `data/genemodels/by_chr/` | exon / CDS structure |
| `tableMaxSpan` | `20_000_000` | wider intervals return a VCF download instead of a table |
| `buildCellsMax` | `2e9` | mirror of the server's build limit, for the run-bar warning only |

`window.SNPTOOLS_MAIZEGDB_BASE`, set before `js/snpfunction-ontology.js` loads, changes the MaizeGDB
host SNPFunction reads. SNPCompare has its own `CFG` (`ibsCompare.php` endpoint, default dataset, a
demo toggle).

Server side, environment variables of the PHP process (`SetEnv` in Apache, `ENV` in Docker):

| Variable | Default | Purpose |
| --- | --- | --- |
| `PYTHON_PATH` | `python3` | interpreter with h5py + numpy that runs `h5_to_vcf.py` |
| `SNPTOOLS_MAX_CELLS` | `2e9` | largest request built, in genotype calls (variants x lines; about 100 Mb for all 933 lines); larger ones are refused with a message |
| `SNPTOOLS_VCF_TTL_HOURS` | `24` | built VCFs are deleted after this long (checked at most every 10 minutes); `0` keeps them |
| `SNPTOOLS_DEBUG` | unset | `1` adds the extractor's command line and output to `processForm.php` replies (local debugging only: they name server paths) |
| `SNPTOOLS_DISTANCE_DIR` | `./distance/` | where `ibsCompare.php` looks for the matrices |
| `SNPTOOLS_SITEQC` | unset | `0` builds every VCF without `NHET`, `NHOM` and `SITEQC`: `h5_to_vcf.py` adds none from the site-QC sidecars and takes out those a store carries itself, so the pages show no Site QC |

---

## Testing

`SNPTools/localdev/harness/` runs the app in jsdom against the real PHP endpoints and stores, then
recomputes what each tool shows independently from the fixture VCFs and data files
(`localdev/harness/README.md`):

```bash
cd SNPTools
make test PHP_BIN=/path/to/php PYTHON_PATH=/path/to/python
```

---

## Data provenance & citation

- Variants and genotypes: Grzybowski MW, Mural RV, Xu G, Turkus J, Yang J, Schnable JC (2023).
  A common resequencing-based genetic marker data set for global maize diversity. *Plant J*
  113(6):1109-1121. https://doi.org/10.1111/tpj.16123 — imputed and linked to USDA GRIN accessions
  by MaizeGDB (release v1.4).
- Passport, origin and trait data: USDA-ARS GRIN-Global.
- Effects: SnpEff 5.2a. DNA language-model scores: PlantCAD1, PlantCAD2, Evo2 7B. Protein
  language-model scores: ESM1b, ESM-2, ESM3, ESM C.
- Protein structures: AlphaFold2, Boltz-2, ESMFold; secondary structure from DSSP; domains from
  Pfam / InterProScan; disorder from IUPred2A.
- Gene Ontology and pathways: MaizeGDB (GO with evidence codes, E2P2 pathways, EnTAP KEGG maps).
- GWAS: Tibbs-Cortes et al. (2024), https://doi.org/10.1101/gr.279027.124, and the studies listed in
  `data/gwas/manifest.json`.
- Base maps: Natural Earth (public domain).

<!-- TODO: add the preferred citation for SNPTools itself. -->

---

## Acknowledgments

Developed by [MaizeGDB](https://maizegdb.org). Funded by USDA-ARS. Uses d3, 3Dmol.js and jsdom (tests).

<!-- TODO: add contributors and award numbers. -->

---

## License

<!-- TODO: add a license. If this is a U.S. government work, state the applicable public terms;
     otherwise add a LICENSE file and reference it here. -->
