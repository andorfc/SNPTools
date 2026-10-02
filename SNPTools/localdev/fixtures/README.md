# Local test fixtures (dataset `zmgrin2026_imp`, release v1.3)

Cut on Ceres on 2026-09-28 with bcftools 1.20 from
`/90daydata/maizegdb/carson/grz2023/release_v1.3/all/grz2023_release_v1.3_all_<chr>.vcf.gz`
(926 Beagle-imputed samples); the 6 default companion lines (ZmG_CM174, ZmG_CML144,
ZmG_CML312, ZmG_CML451, ZmG_PHT69, ZmG_R109B) were appended from
`companion/grz2023_release_v1.3_direct_calls_8lines.vcf.gz` (GT only, matched on
CHROM/POS/REF/ALT; all 3,495 sites matched). Sample columns = the 932 `sample_id`s of
`data/zmgrin2026_samples.tsv`. INFO is the release INFO (AC/AN/AF/DR2/IMP); the
annotation keys SNPTools shows (GENEMODEL, TYPE, EFFECT, plantcad*/ESM* scores) are not
merged yet, so those columns read "—".

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
