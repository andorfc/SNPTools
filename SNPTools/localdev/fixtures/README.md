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
llr_esm2 = store ESM-2 650M, i.e. esm2_store_score / llr_esm3, 1 decimal). No MQ/CVC/CVP/MAXR2
and no PlantCAD scores (pending). Inputs for re-running it are in `annotation/`:

    for f in zmgrin2026_v1.4_chr*_testregions.vcf.gz; do
      python3 ../../tools/annotate_release_info.py --vcf $f \
        --snpeff annotation/schnable_scored_testregions.sites.vcf.gz \
        --esm annotation/grz2023_missense_esm_testregions.tsv.gz --out /tmp/a.vcf.gz && mv /tmp/a.vcf.gz $f
    done

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
Queries outside these windows return "No variants"; chromosomes 3, 4 and 6 have no store.
