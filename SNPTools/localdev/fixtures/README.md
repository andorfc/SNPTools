# Local test fixtures (dataset `zmgrin2026_imp`, release v1.4)

Cut on Ceres on 2026-09-28 with bcftools 1.20 from
`/90daydata/maizegdb/carson/grz2023/release_v1.3/all/grz2023_release_v1.3_all_<chr>.vcf.gz`
(926 Beagle-imputed samples; records and genotypes are identical in release v1.4); the 6
default companion lines of v1.3 (ZmG_CM174, ZmG_CML144, ZmG_CML312, ZmG_CML451, ZmG_PHT69,
ZmG_R109B) were appended from `companion/grz2023_release_v1.3_direct_calls_8lines.vcf.gz`
(GT only, matched on CHROM/POS/REF/ALT; all 3,495 sites matched). On 2026-09-29 the NAM founder
**ZmG_CML103** (GRIN PI 690319) was appended as the last column from
`release_v1.4/companion/grz2023_release_v1.4_direct_calls_8lines.vcf.gz` (bcftools 1.20 on
Ceres; GT only, same matching; 3,495/3,495 sites matched: 2,833 0/0, 67 0/1, 377 1/1, 218 ./.).
The per-sample FT field is not applied: its 41 GRZHARD genotypes are already ./. and its 14
LowQual genotypes are 0/0 reference calls, which are kept. The demo regions do not overlap the
chr10 30-130 Mb segment flagged in the release's CML103 identity caveat. Sample columns = the
933 `sample_id`s of `data/zmgrin2026_samples.tsv`, in that order (CML103 last). INFO carries the SNPVersity annotation keys written by `tools/annotate_release_info.py` (2026-09-29):
TYPE/EFFECT/GENEMODEL/SUB copied from the MaizeGDB Schnable scored VCFs (SnpEff 5.2a; the
3,495/3,495 sites matched on CHROM/POS/REF/ALT), MAF from the 933 genotypes of the file, and
ESM1_score/ESM2_score/ESM3_score for the 133 missense sites (grz2023_missense_esm llr_esm1b /
llr_esm2 = store ESM-2 650M, i.e. esm2_store_score / llr_esm3, 1 decimal). No MQ/CVC/CVP (not
available for the Grzybowski call set), no MAXR2 and no PlantCAD scores: the full-chromosome
builds add those with `--maxr2` and `--dna-scores` (PlantCAD1/2 and Evo2 rounded to 0.1,
`--pc-decimals 1`). Inputs for re-running it are in `annotation/`:

    for f in zmgrin2026_v1.4_chr*_testregions.vcf.gz; do
      python3 ../../tools/annotate_release_info.py --vcf $f \
        --snpeff annotation/schnable_scored_testregions.sites.vcf.gz \
        --esm annotation/grz2023_missense_esm_testregions.tsv.gz --out /tmp/a.vcf.gz && mv /tmp/a.vcf.gz $f
    done

`chr2_store/zmgrin2026_v1.4.2_chr2_4491424_4499434.annotated.vcf.gz`: the GWAS window of the
annotated chr2 release VCF, used by the checks when the full chr2 store is installed. Cut on Ceres
on 2026-10-05 with bcftools 1.20 (`bcftools view -r chr2:4491424-4499434 -Oz`) from the release
v1.4.2 file
`/90daydata/maizegdb/carson/grz2023/release_v1.4.2/annotated_933/zmgrin2026_v1.4.2_chr2_933.annotated.vcf.gz`
(md5 6b83f45f72661a57c6c7f9796e497a94, as in the build's `logs/MD5SUMS.txt`): the release VCF, not
the store, so the checks compare what the app shows with an independent source. Same 249 sites and
genotypes as the demo fixture (the 7 companion lines stand in another column order). Its INFO adds
plantcad1/plantcad2 and evo2 (the window's 158 SNPs), ESMC_score (10), MAXR2 (233) and NHET, NHOM
and SITEQC (all 249); its header carries the release's 933 `##SAMPLE` lines. bcftools writes whole
numbers without decimals (MAF 0, a score of -1) where the store, built from the annotator's own
output, reads 0.0000 and -1.0: the same values (41 MAF, 44 score and 2 MAXR2 values of the window
differ in text only).
It replaces `zmgrin2026_v1.4_chr2_4491424_4499434.annotated.vcf.gz`, cut on 2026-09-30 from the
release v1.4 build (`grz2023/snptools_build/chr2/zmgrin2026_v1.4_chr2_933.annotated.vcf.gz`, store
md5 9b8a52df8d08df7a875f3c84f4bf7d3b), in which evo2 covered the 120 SNPs within 1 kb of a gene
and INFO had no site-QC fields; every key of that file has the same value in this one, and the 38
sites that gained evo2 are upstream (25) and downstream (13) gene variants.

`annotation/schnable_scored_testregions.sites.vcf.gz`: columns 1-8 of Atlas
`/90daydata/maizegdb/carson/grz2023_scoring/protein/vcf/chr<N>_schnable_scored.vcf.gz` in the
fixture windows. `annotation/grz2023_missense_esm_testregions.tsv.gz`: the 147 rows of
grz2023_missense_esm.tsv.gz at fixture sites. `mgdb2026_{hq,hc}_chr2_4491424_4499434.sites.vcf.gz`:
columns 1-8 of Atlas `MaizeGDB2026_AI_scores/VCF/chr2_high_{quality,coverage}.vcf` in the chr2
window (414 / 942 sites, no genotypes), used by the harness to render the annotation columns of
the MaizeGDB 2026 sets without a local store.

| Region | Why | Sites |
|---|---|---|
| chr9:12,836,008-12,846,000 | Zm00001eb374090 ±2 kb (SNPGeo default gene) | part of 451 on chr9 |
| chr9:13,118,306-13,124,164 | Zm00001eb374230 ±2 kb | part of 451 on chr9 |
| chr2:4,491,424-4,499,434 | Zm00001eb067740 ±2 kb | 249 |
| chr1:280,588,357-280,596,023 | Zm00001eb056510 ±2 kb | 288 |
| chr5:90,719,578-90,729,950 | Zm00001eb233650 ±2 kb | 433 |
| chr7:123,683,735-123,693,964 | Zm00001eb313510 ±2 kb | 289 |
| chr10:129,631-131,683 / 216,406-222,251 | Zm00001eb404740, Zm00001eb404760 | part of 1,517 on chr10 |
| chr10:9,788,000-9,826,500 | SNPVersity default window (no annotated gene) | part of 1,517 on chr10 |
| chr8:163,450,112-163,454,880 | old GENE_MODELS interval for Zm00001eb374090 (no annotated gene) | 268 |

Gene intervals are from `gff/genes_data.serialized` (what lookupGeneModel.php returns).
With only the test stores built from these files (`make local-store`), queries outside these
windows return "No variants" and chromosomes 3, 4 and 6 have no store. Since 2026-10-02 the full
stores of all ten chromosomes are installed locally (since 2026-10-04 the release v1.4.2 builds);
the checks then read each window from the
store (`harness/store_windows.py`), after confirming it matches these fixtures site by site.
