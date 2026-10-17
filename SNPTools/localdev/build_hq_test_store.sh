#!/usr/bin/env bash
# Build small local test stores for dataset zmgrin2026_hq (the High Quality set of release v1.4.3)
# when the full stores (1.0-1.9 GB each, 12.9 GB for the ten) do not fit: 500-kb windows around the
# example genes, the SNPVersity default window and the SNPCurate genes, cut from the HQ release VCFs
# on Ceres with bcftools and stored with tools/vcf_to_h5.py (162,764 sites x 933 samples, ~100 MB).
# Output: ../hdf5/grin2026/zmgrin2026_<chr>_HQ.h5 (ignored by hdf5/.gitignore). The HQ site-QC
# sidecars beside them describe the full stores, so h5_to_vcf.py notes the mismatch and reads NHET,
# NHOM and SITEQC from each row's INFO, which the HQ stores carry. Queries outside the windows
# return "No variants". A store larger than 50 MB (a full one) is never replaced.
# Needs ssh access to Ceres and a python with h5py + numpy:
#   PYTHON_PATH=/path/to/python ./build_hq_test_store.sh
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"; root="$(dirname "$here")"
PY="${PYTHON_PATH:-python3}"
HOST="${CERES_HOST:-ceres}"
SRC="${HQ_VCF_DIR:-/90daydata/maizegdb/carson/grz2023/release_v1.4.3_HQ/vcf}"
SSH=(ssh -n -o PKCS11Provider=none -o BatchMode=yes "$HOST")
tmp="$(mktemp -d "${TMPDIR:-/tmp}/hqcut.XXXXXX")"
# anchor +-250 kb, overlapping windows merged
while read -r chr list; do
  out="$root/hdf5/grin2026/zmgrin2026_${chr}_HQ.h5"
  if [ -f "$out" ] && [ "$(wc -c < "$out")" -gt 50000000 ]; then
    echo "skip $chr: $(basename "$out") is a full store"; continue
  fi
  vcf="$tmp/zmgrin2026_v1.4.3_${chr}_933.HQ.testwindows.vcf.gz"
  "${SSH[@]}" "module load bcftools/1.20 >/dev/null 2>&1; bcftools view -r $list $SRC/zmgrin2026_v1.4.3_${chr}_933.HQ.vcf.gz -Oz" > "$vcf"
  rm -f "$out"
  "$PY" "$root/tools/vcf_to_h5.py" "$vcf" "$out"
  rm -f "$vcf"
done <<'EOF'
chr1 chr1:47674442-48198490,chr1:280340357-280840357
chr2 chr2:4243424-4743424,chr2:237250176-237750176
chr3 chr3:151271664-151771664
chr4 chr4:43180007-43680007,chr4:46397932-46897932,chr4:184898556-185398556,chr4:205164194-205664194
chr5 chr5:6712520-7212520,chr5:90471578-90971578
chr6 chr6:91393082-91893082,chr6:114933677-115433677
chr7 chr7:5757514-6257514,chr7:123435735-123935735
chr8 chr8:126428681-126928681,chr8:142337452-142837452,chr8:155802116-156302116,chr8:163200112-163700112
chr9 chr9:12588008-13370306,chr9:158240315-158740315
chr10 chr10:1-470251,chr10:9538000-10076500,chr10:95925873-96425873,chr10:138437404-138937404
EOF
rmdir "$tmp"
