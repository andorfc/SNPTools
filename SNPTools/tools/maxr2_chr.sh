#!/bin/bash
# maxR2 for one chromosome of release v1.4 (933 default samples), as in the MaizeGDB 2026 stores:
# highest pairwise LD r2 (PLINK 1.9 --r2, genotypic allele-count correlation) with any variant 400-5,000 bp away.
# No MAF / r2 / missingness filtering. Output: chrN.maxr2.tsv.gz  (chr pos ref alt MAXR2 [plink value as printed])
set -eo pipefail
C=$1; SD=$2
G=/90daydata/maizegdb/carson/grz2023; R=$G/release_v1.4; B=$G/snptools_build/$C; O=$G/snptools_build/maxr2; mkdir -p $B $O
T=/tmp/mx_$C; mkdir -p $T
M=$B/zmgrin2026_v1.4_${C}_933.gt.vcf.gz
if [ ! -s $M.tbi ]; then   # same 933-sample merge as stageB_ceres.sh (reused by the SNPTools build)
  bcftools annotate --threads 2 -x ^FORMAT/GT $R/all/grz2023_release_v1.4_all_$C.vcf.gz -Oz -o $T/imp.vcf.gz; bcftools index -t $T/imp.vcf.gz
  bcftools view -s ^ZmG_KI9 -r $C $R/companion/grz2023_release_v1.4_direct_calls_8lines.vcf.gz -Ou | \
    bcftools annotate -x INFO,FILTER,^FORMAT/GT -Oz -o $T/comp.vcf.gz; bcftools index -t $T/comp.vcf.gz
  bcftools merge --threads 2 -m none $T/imp.vcf.gz $T/comp.vcf.gz -Oz -o $M.tmp.vcf.gz
  mv $M.tmp.vcf.gz $M; bcftools index -t $M
fi
N=$(bcftools index -n $M)
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' $M -Oz -o $T/id.vcf.gz
plink --vcf $T/id.vcf.gz --double-id --allow-extra-chr --chr-set 10 no-xy no-mt --vcf-half-call m \
      --keep-allele-order --make-bed --out $T/g --threads ${THREADS:-4} > $T/plink_bed.log 2>&1
mkfifo $T/r2.ld
# reducer: max r2 per variant over partners at 400..5000 bp; also keep the first 300,000 pair lines for validation
awk -v S=$T/pairs_sample.txt 'NR==1{next} { if (NR<=300001) print > S
       d=$5-$2; if (d<0) d=-d; if (d<400 || $7=="nan" || $7=="-nan") next
       if (!($3 in m) || $7+0>m[$3]+0) m[$3]=$7; if (!($6 in m) || $7+0>m[$6]+0) m[$6]=$7 }
     END{for (k in m) print k"\t"m[k]}' $T/r2.ld > $T/max.txt &
AWK=$!
plink --bfile $T/g --allow-extra-chr --chr-set 10 no-xy no-mt --r2 --ld-window 999999 --ld-window-kb 5 --ld-window-r2 0 \
      --out $T/r2 --threads ${THREADS:-4} > $T/plink_r2.log 2>&1 || { timeout 10 bash -c "echo > $T/r2.ld" || true; kill $AWK 2>/dev/null || true; cp $T/plink_r2.log $O/$C.plink_r2.log; echo "plink --r2 failed"; exit 7; }
wait $AWK
# emit in VCF order; variants without any partner in 400-5000 bp (or monomorphic) get no MAXR2
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $M | awk -F'\t' 'BEGIN{OFS="\t"; print "chr","pos","ref","alt","MAXR2"}
   NR==FNR{m[$1]=$2; next} {k=$1":"$2":"$3":"$4; print $1,$2,$3,$4,(k in m? m[k] : "")}' $T/max.txt - | gzip -5 > $O/$C.maxr2.tsv.gz
NV=$(gzip -dc $O/$C.maxr2.tsv.gz | tail -n +2 | awk -F'\t' '$5!=""' | wc -l)
NP=$(grep -oE '[0-9]+ variants and [0-9]+ (people|samples) pass' $T/plink_bed.log | head -1 || true)
echo -e "$C\t$N\t$NV\t$NP" > $O/$C.maxr2.summary.tsv
cp $T/plink_bed.log $O/$C.plink_bed.log; cp $T/plink_r2.log $O/$C.plink_r2.log; cp $T/pairs_sample.txt $O/$C.pairs_sample.txt
[ "$C" = chr10 ] && cp $T/g.bed $T/g.bim $T/g.fam $O/ 2>/dev/null || true
rm -rf $T
