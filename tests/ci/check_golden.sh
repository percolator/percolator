#!/usr/bin/env bash
# Checks that the results of percolator and qvality are byte for byte what
# they were: runs option_matrix.sh with synthetic inputs of 8,000 PSMs and
# compares the SHA-256 hash of every output file (stderr, which holds
# timings and paths, is left out) with the manifest
# tests/ci/option_matrix.sha256.
#
# usage: tests/ci/check_golden.sh [--update | --report | --strict]
#                                 <build-dir> <out-dir>
#   --update  writes the manifest from this run instead of comparing. Use it
#             in a commit that changes results on purpose, or adds runs to
#             option_matrix.sh, with a Release build at default flags on
#             x86-64 Linux, and say in the commit message why the hashes
#             change.
#   --report  lists the differences but exits with status 0.
#   --strict  exits with status 1 if a file differs.
#   Without an option the check is strict on x86-64 Linux and a report
#   elsewhere: the manifest was made there, and another processor
#   architecture or math library may round the last digit differently.
#   The hashes of this run are written to <out-dir>/hashes.sha256.
# Exits with status 1 if a run fails or, when strict, a file differs.
set -u
mode=auto
case "${1:-}" in
  --update) mode=update; shift ;;
  --report) mode=report; shift ;;
  --strict) mode=strict; shift ;;
esac
if [ $# -ne 2 ]; then
  sed -n '2,21p' "$0" >&2
  exit 2
fi
BUILD=$1
OUT=$2
HERE=$(cd "$(dirname "$0")" && pwd)
MANIFEST=$HERE/option_matrix.sha256
if [ "$mode" = auto ]; then
  if [ "$(uname -s)" = Linux ] && [ "$(uname -m)" = x86_64 ]; then
    mode=strict
  else
    mode=report
  fi
fi
if command -v sha256sum > /dev/null 2>&1; then
  sha256() { sha256sum "$@"; }
elif command -v shasum > /dev/null 2>&1; then
  sha256() { shasum -a 256 "$@"; }
else
  echo "error: neither sha256sum nor shasum found" >&2
  exit 2
fi

# the manifest holds for this input size only
export MATRIX_PSMS=8000
rm -rf "$OUT"
mkdir -p "$OUT"
bash "$HERE/option_matrix.sh" "$BUILD" "$OUT" > "$OUT/matrix.log" 2>&1
matrix_status=$?
if [ "$matrix_status" -eq 0 ]; then
  tail -n 1 "$OUT/matrix.log"
else
  grep -v '^ok ' "$OUT/matrix.log"
  echo "error: option_matrix.sh failed; the hashes are not compared" >&2
  exit 1
fi

HASHES=$OUT/hashes.sha256
(cd "$OUT" && find . -type f ! -name stderr.txt ! -name matrix.log \
    ! -name hashes.sha256 | sed 's|^\./||' | LC_ALL=C sort |
  while read -r f; do sha256 "$f"; done) |
  awk '{ sub(/^\*/, "", $2); print $1 "  " $2 }' > "$HASHES"
total=$(wc -l < "$HASHES" | tr -d ' ')

if [ "$mode" = update ]; then
  cp "$HASHES" "$MANIFEST"
  echo "wrote $MANIFEST ($total files)"
  exit 0
fi
if [ ! -f "$MANIFEST" ]; then
  echo "error: $MANIFEST not found" >&2
  exit 2
fi

differences=$(awk '
  NR == FNR { want[$2] = $1; next }
  { got[$2] = $1 }
  END {
    for (f in want) {
      if (!(f in got)) print "MISSING " f
      else if (got[f] != want[f]) print "DIFF " f
    }
    for (f in got) if (!(f in want)) print "NEW " f
  }' "$MANIFEST" "$HASHES" | LC_ALL=C sort -k 2)
if [ -z "$differences" ]; then
  echo "$total files, all as in the manifest"
  exit 0
fi
printf '%s\n' "$differences"
n=$(printf '%s\n' "$differences" | wc -l | tr -d ' ')
echo "$total files, $n not as in the manifest (DIFF: other content, MISSING:" \
  "not written, NEW: not in the manifest)"
if printf '%s\n' "$differences" | grep -q ' inputs/'; then
  echo "The synthetic inputs differ: tests/benchmark/make_synthetic_pin.py or" \
    "awk gives other files here, so the outputs cannot be compared."
fi
if [ "$mode" = report ]; then
  echo "Reported only ($(uname -s) $(uname -m)); not counted as a failure."
  exit 0
fi
echo "If the change of results is intended, write the manifest again with:" \
  "bash tests/ci/check_golden.sh --update <build-dir> <out-dir>"
exit 1
