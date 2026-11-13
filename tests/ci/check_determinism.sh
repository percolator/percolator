#!/usr/bin/env bash
# Checks that the results of percolator do not depend on the number of
# threads: runs option_matrix.sh with --num-threads 1 and with
# --num-threads 4 and compares every output file byte for byte (stderr,
# which holds timings, is not compared).
#
# usage: tests/ci/check_determinism.sh <build-dir> <out-dir>
#   MATRIX_PSMS sets the number of PSMs of the synthetic inputs (default
#   100000, so that the parallel sorts and loops split the work).
# Exits with status 1 if a run fails or an output differs.
set -u
if [ $# -ne 2 ]; then
  sed -n '2,10p' "$0" >&2
  exit 2
fi
BUILD=$1
OUT=$2
HERE=$(cd "$(dirname "$0")" && pwd)
export MATRIX_PSMS=${MATRIX_PSMS:-100000}

status=0
bash "$HERE/option_matrix.sh" "$BUILD" "$OUT/threads1" --num-threads 1 || status=1
bash "$HERE/option_matrix.sh" "$BUILD" "$OUT/threads4" --num-threads 4 || status=1

compared=0
differing=0
for f in $(cd "$OUT/threads1" && find . -type f ! -name stderr.txt | sort); do
  compared=$((compared + 1))
  if ! cmp -s "$OUT/threads1/$f" "$OUT/threads4/$f"; then
    differing=$((differing + 1))
    echo "DIFF ${f#./}"
  fi
done
echo "$compared files compared, $differing differ between 1 and 4 threads"
[ "$differing" -eq 0 ] || status=1
exit $status
