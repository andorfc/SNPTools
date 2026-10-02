#!/usr/bin/env bash
# Start a local SNPTools instance: PHP built-in web server on 127.0.0.1:${PORT:-8765},
# document root = SNPTools/ (so processForm.php and lookupGeneModel.php run for real).
#   PHP_BIN      php executable (default: php on PATH)
#   PYTHON_PATH  python with h5py + numpy, used by processForm.php -> h5_to_vcf.py
#   PORT         default 8765
# Logs to localdev/server.log; PID in localdev/server.pid.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"; root="$(dirname "$here")"
PHP="${PHP_BIN:-php}"; PORT="${PORT:-8765}"
export PYTHON_PATH="${PYTHON_PATH:-python3}"
if [ -f "$here/server.pid" ] && kill -0 "$(cat "$here/server.pid")" 2>/dev/null; then
  echo "already running (pid $(cat "$here/server.pid")) on port $(cat "$here/server.port" 2>/dev/null || echo "$PORT")"; exit 0
fi
"$PYTHON_PATH" -c 'import h5py, numpy' || { echo "PYTHON_PATH=$PYTHON_PATH lacks h5py/numpy" >&2; exit 1; }
ls "$root"/hdf5/version3/zmgrin2026_chr*_impute.h5 >/dev/null 2>&1 || echo "note: no zmgrin2026 HDF5 store yet - run localdev/build_test_store.sh" >&2
mkdir -p "$root/vcf"
cd "$root"
nohup "$PHP" -d max_execution_time=0 -S "127.0.0.1:$PORT" -t "$root" > "$here/server.log" 2>&1 &
echo $! > "$here/server.pid"; echo "$PORT" > "$here/server.port"
sleep 1
if kill -0 "$(cat "$here/server.pid")" 2>/dev/null; then
  echo "SNPTools running at http://127.0.0.1:$PORT/index.html (pid $(cat "$here/server.pid"))"
else
  echo "server failed to start; see $here/server.log" >&2; cat "$here/server.log" >&2; exit 1
fi
