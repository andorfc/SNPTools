#!/usr/bin/env bash
# Stop the local SNPTools instance started by start.sh.
here="$(cd "$(dirname "$0")" && pwd)"
if [ -f "$here/server.pid" ]; then
  pid="$(cat "$here/server.pid")"
  if kill -0 "$pid" 2>/dev/null; then kill "$pid" && echo "stopped pid $pid"; else echo "pid $pid not running"; fi
  rm -f "$here/server.pid" "$here/server.port"
else
  echo "no server.pid - not running (check: lsof -iTCP:${PORT:-8877} -sTCP:LISTEN)"
fi
