#!/usr/bin/env bash
# Build the small local HDF5 test store for dataset zmgrin2026_imp from the
# fixture VCFs in fixtures/ (release v1.4: 926 imputed + 7 companion lines = 933
# sample ids; 8 example genes and 2 test windows, 3,495 sites; see fixtures/README.md).
# Output: ../hdf5/grin2026/zmgrin2026_<chr>_impute.h5 and its site-QC sidecar
# zmgrin2026_<chr>_impute.siteqc.h5 (tools/build_site_qc.py; both ignored by hdf5/.gitignore).
# Needs a python with h5py + numpy:   PYTHON_PATH=/path/to/python ./build_test_store.sh
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"; root="$(dirname "$here")"
PY="${PYTHON_PATH:-python3}"
mkdir -p "$root/hdf5/grin2026"
for f in "$here"/fixtures/zmgrin2026_v1.4_chr*_testregions.vcf.gz; do
  chr="$(basename "$f" | sed -E 's/^zmgrin2026_v1\.4_(chr[0-9]+)_testregions\.vcf\.gz$/\1/')"
  out="$root/hdf5/grin2026/zmgrin2026_${chr}_impute.h5"
  # never overwrite a full-chromosome store (e.g. the Ceres chr2 build) with the demo cut
  if [ -f "$out" ] && [ "$(wc -c < "$out")" -gt 50000000 ] && [ "${FORCE:-0}" != 1 ]; then
    echo "skip $chr: $(basename "$out") is a full store ($(( $(wc -c < "$out") / 1000000 )) MB); FORCE=1 to replace it with the demo cut"
    continue
  fi
  "$PY" "$root/tools/vcf_to_h5.py" "$f" "$out"
  "$PY" "$root/tools/build_site_qc.py" "$out"
done
