#!/usr/bin/env bash
# Build the small local HDF5 test store for dataset zmgrin2026_imp from the
# fixture VCFs in fixtures/ (release v1.4: 926 imputed + 7 companion lines = 933
# sample ids; 8 example genes and 2 test windows, 3,495 sites; see fixtures/README.md).
# Output: ../hdf5/version3/zmgrin2026_<chr>_impute.h5 (ignored by hdf5/.gitignore).
# Needs a python with h5py + numpy:   PYTHON_PATH=/path/to/python ./build_test_store.sh
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"; root="$(dirname "$here")"
PY="${PYTHON_PATH:-python3}"
mkdir -p "$root/hdf5/version3"
for f in "$here"/fixtures/zmgrin2026_v1.4_chr*_testregions.vcf.gz; do
  chr="$(basename "$f" | sed -E 's/^zmgrin2026_v1\.4_(chr[0-9]+)_testregions\.vcf\.gz$/\1/')"
  "$PY" "$root/tools/vcf_to_h5.py" "$f" "$root/hdf5/version3/zmgrin2026_${chr}_impute.h5"
done
