# Maize SNPTrait / SNPGeo data builders

All builders are Python 3.8+ standard library unless noted. Outputs mirror `SNPTools/`.

| Script | Output | Inputs |
|---|---|---|
| `build_maize_samples.py` | `data/zmgrin2026_samples.tsv` | T0 population members, `grz2023_grin_crosswalk.tsv`, optional `--release-manifest` (vcf_name -> sample_id as written to the release VCFs) and `--extra-samples` (separately called lines) |
| `build_maize_trait_catalog.py` | `js/zmgrin.catalog.js` (`SNP_CATALOG.families.zmgrin2026`, `SNPTRAIT_SCHEMA.zmgrin2026`), `data/traits/zmgrin2026.traits.json` | samples TSV, T5 passport/geo, T4b trait summary, T4c trait dictionary, crosswalk |
| `build_maize_snpgeo_data.py` | `js/snpgeo.regions.js` (and `data/geo/countries.geo.json`) | Natural Earth admin-0, samples TSV, T5, crosswalk |
| `build_geo_layers.py` | `data/geo/countries.geo.json`, `data/geo/admin1_na.geo.json` | Natural Earth v5.1.2 admin-0 110m, admin-1 10m (see `data/geo/PROVENANCE.md`) |
| `vcf_to_h5.py` (h5py, numpy) | `hdf5/grin2026/zmgrin2026_<chr>_impute.h5` | one release VCF per chromosome |
| `build_snpcurate.py` (h5py, numpy) | `js/snpcurate.data.js`; with `--table`, the benchmark table (.tsv + .md) | `data/curate/snpcurate.source.json`, the stores and sidecars, `js/zmgrin.catalog.js`, `data/traits/zmgrin2026.traits.json`, `gff/genes_index.txt` |
| `build_site_qc.py` (h5py, numpy) | `hdf5/grin2026/zmgrin2026_<chr>_impute.siteqc.h5`; with `--summary`, `data/qc/zmgrin2026.siteqc.summary.json`, `data/qc/zmgrin2026.lineqc.tsv`, `js/zmgrin.lineqc.js` | the stores (read-only) |

The files committed in this branch were built for release v1.4 (933 samples = 926 imputed
+ 7 companion direct-call lines; v1.4 adds NAM founder CML103, `ZmG_CML103`, PI 690319). The
catalogue and trait side-file were rebuilt with `build_maize_trait_catalog.py --bins-from
data/traits/zmgrin2026.traits.json` so the v1.3 tertile labels stay unchanged; the CML103 trait
values come from the v1.4 genotype-phenotype linkage table (see its `source` field). The region
record for CML103 was taken from `build_maize_snpgeo_data.py` and added to the v1.3 file (the
builder in this tree also resolves ISO3 for 27 class-C lines and writes "United States of America";
regenerating the whole file would change those 776 records, which is left as a follow-up). Run order and pinned sources are in
`../Makefile` (`make help`). The prototype validators from the design phase
(`validate_maize_port.py`, `smoke_test_snptools.js`) ran the unmodified fusarium modules
with proposed patches; they are superseded by `localdev/harness/` now that the changes are
applied. `validate_maize_port.py` is kept for its structural checks of the catalogue shape.

## build_site_qc.py (site QC sidecars and line QC)

`make site-qc` (= `build_site_qc.py --jobs 4`, then `--summary`). For each store
`<name>.h5` it writes `<name>.siteqc.h5` in the same folder, never opening the store for writing.
Datasets, in the store's site order, chunked (65,536,) with gzip level 4 like the store: `POS` int64
(a copy, for alignment checks), `NHET`, `NHOM`, `NMISS` int16 (heterozygous, homozygous-alternate and
missing calls over all the store's samples). Attributes: `format` = `snptools-siteqc-1`, `store`,
`store_bytes`, `n_sites`, `n_samples`, `built`. The file is written under a temporary name and
renamed. A sidecar is valid while its `n_sites`, `store_bytes` and `POS` match its store;
`--check` (`make site-qc-check`) verifies that and exits 1 on any mismatch. `h5_to_vcf.py` uses a
sidecar while its `n_sites` and `store_bytes` match (positions are not re-read per request) and
appends `NHET`, `NHOM` and `SITEQC` to every row's INFO that does not have them already (the release
v1.4.2 and later stores carry them in their own rows, see `annotate_release_info.py` below); otherwise it
prints a `Note:` line and writes the rows as stored. `SNPTOOLS_SITEQC=0` writes every VCF without
the three, whatever their source: nothing is merged, the fields a store carries itself are taken out
of its rows, and the header has no `##INFO` line for them.

Only counts are stored; the class is derived by one rule (first match wins), written once per
language: `site_qc_codes` here and in `h5_to_vcf.py`, `site_qc` in `annotate_release_info.py`:

| Code | Rule | Group |
|---|---|---|
| `NO_CARRIER` | NHET + NHOM = 0 | no carrier |
| `HET_ONLY` | NHOM = 0 | flagged |
| `HET_EXCESS` | NHET > NHOM | flagged |
| `HET_ELEVATED` | 4 x NHET >= NHET + NHOM | caution (usable) |
| `PASS` | everything else | usable |

`--summary` reads the stores and sidecars again and writes the class counts per chromosome and in
total, and the per-line table. A clean site is `PASS` with NHOM >= 3 (16,463,157 sites in release
v1.4). A line's share is its heterozygous calls at clean sites over the clean sites: class H
(heterozygous sample) above 0.05, E (elevated heterozygosity) above 0.02, I (inbred) otherwise.
Release v1.4: 20 H, 26 E, 887 I.

## build_snpcurate.py (SNPCurate)

`make curate`. The registry is `data/curate/snpcurate.source.json`, edited by hand (one entry per
published allele: gene and symbol, label, kind `causal` / `published_marker` / `tag`, variant type,
status `site` / `not_a_site` / `not_locatable` / `structural` / `unreliable` / `not_annotated`, the site
(chr, pos, ref, alt) or the position or interval where known, `site_gene` when the site is annotated in
another gene, evidence `validated` / `associated` / `tag`, the trait and its GRIN descriptor, references
with DOIs, a note, `needs_review`, `tagged_by`, and an `expect` block of het / hom / class for tests).
For each site entry the builder finds exactly one site (else it stops) and computes the carriers
(per subpopulation and country of origin, with each group's total), the Site QC class, the allele
frequency, MAF and maxR², the consequence and substitution of `site_gene` or `gene` read as the app
reads them, the seven scores, the combined score and priority (as `Data.impactPriority`), the trait by
genotype, and for a missense entry its rank by each score among the gene's missense variants with a
carrier. A non-site entry is checked to have no site at its position or codon. It stops when a count
differs from the `expect` block. Marks: gold (causal, validated, a site, nothing to review), outline
(any other site), grey (not a site). Release v1.4: 21 entries, 14 sites; gold 2, outline 12, grey 7.
`--table <out.tsv>` also writes the benchmark table, one row per entry, and the same as Markdown.

## annotate_release_info.py (SNPVersity annotation INFO for zmgrin2026_imp)

Adds TYPE/EFFECT/GENEMODEL/SUB (from the MaizeGDB Schnable scored VCFs, SnpEff 5.2a), MAF (from
the release genotypes) and ESM1_score/ESM2_score/ESM3_score (missense ESM table; ESM2 = store
ESM-2 650M = esm2_store_score) to a release VCF before `vcf_to_h5.py`. Streams both sorted VCFs
(one chromosome per run); run it per chromosome on the full release to annotate the 46M-site
store. MQ/CVC/CVP are deliberately not written (not available for the Grzybowski et al. 2023 call
set; see the docstring). `--dna-scores <tsv>` (chr,pos,ref,alt,plantcad1_score,plantcad2_score,
plantcad2_onepass_score,evo2_score; position-sorted) adds the Atlas PlantCAD1, PlantCAD2 (masked,
512 bp), PlantCAD2 one-pass and Evo2 scores (all SNPs; until release v1.4.2 the Evo2 table held the
genic +/-1 kb SNPs only), rounded to 1 decimal like ESM (`--pc-decimals`, default 1; full precision
stays in the score tables). The one-pass score (release v1.4.3) is ln P(alt)/P(ref) at the unmasked
base with 8,192-bp windows; it is not on the scale of plantcad2_score (its median is about -3.4
against -0.7). Any of the four score columns may be present: only those are written and declared,
so a v1.4.2 table without the one-pass column still works. ESMC_score (llr_esmc) is written with the other ESM scores. `--maxr2 <tsv>`
(chr,pos,ref,alt,MAXR2; position-sorted) adds MAXR2: the highest PLINK 1.9 --r2 of the site with
any variant 400-5,000 bp away, from the 933 release genotypes with no MAF/missingness/r2 filtering,
written with up to 6 decimals (1.0, 0.509182) like the MaizeGDB 2026 stores. NHET, NHOM and SITEQC
(site QC, the rule of `build_site_qc.py`) are counted from the file's own genotypes and written last;
a store built from such a VCF carries them itself, and `h5_to_vcf.py` then adds nothing to its rows. The rebuilt
full-chromosome stores on Ceres use both options. The installed stores are the release v1.4.3 builds
(2026-10-07, `grz2023/snptools_build_v1.4.3/chr<N>/` on Ceres; the annotated VCFs are in
`grz2023/release_v1.4.3/annotated_933/`): an annotation-only update of v1.4 with the same lines,
sites and genotypes (v1.4.3 adds plantcad2_onepass_score to v1.4.2, every other field
byte-identical), in which 43,296,332 SNPs carry PlantCAD1, PlantCAD2, the PlantCAD2 one-pass score
and Evo2 and all
46,054,265 sites carry NHET, NHOM and SITEQC. Each store was built from the annotator's own output,
so its INFO reads `MAF=0.0000` or a score of `-1.0` where the bgzipped release VCF, rewritten by
bcftools, reads `0` and `-1`: the same values. The MAXR2 table comes from `maxr2_chr.sh <chr>`
(Ceres: bcftools + PLINK 1.9 --r2 on the 933-sample merge, max over partners 400-5,000 bp away);
`check_r2.py` recomputes 300 random pairs from the PLINK .bed as an independent check.

Natural Earth: https://github.com/nvkelso/natural-earth-vector tag v5.1.2
(f1890d9f152c896d250a77557a5751a93d494776), public domain.
