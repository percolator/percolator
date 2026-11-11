#!/usr/bin/env bash
# Checks that percolator and qvality end with exit status 1 and an error
# message, and not with status 0 or a crash, on input they cannot use and
# on output they cannot write: unknown and invalid options, missing, empty
# and malformed input, inputs without targets or decoys, output files that
# are the input, lie in a missing directory or are on a full disk
# (/dev/full, where it exists).
#
# usage: tests/ci/exit_status.sh <build-dir> <out-dir>
#   PYTHON names the Python 3 interpreter.
# Exits with status 1 if a case ends with another status than expected,
# without the expected message, or with a sanitizer report.
set -u
if [ $# -ne 2 ]; then
  sed -n '2,12p' "$0" >&2
  exit 2
fi
BUILD=$1
OUT=$2
SRC=$(cd "$(dirname "$0")/../.." && pwd)
PERC=$(cd "$BUILD" && pwd)/src/percolator
QVAL=$(cd "$BUILD" && pwd)/src/qvality/qvality
QD=$SRC/data/qvality
PY=${PYTHON:-python3}

for exe in "$PERC" "$QVAL"; do
  if [ ! -x "$exe" ]; then
    echo "error: $exe not found" >&2
    exit 2
  fi
done

mkdir -p "$OUT"
# the cases run inside the output directory, with relative file names
cd "$OUT" || exit 2
"$PY" "$SRC/tests/benchmark/make_synthetic_pin.py" -n 600 -f 10 \
  --mode concatenated --seed 1 -o good.pin || exit 2
: > empty.pin
head -n 1 good.pin > header_only.pin
awk -F '\t' 'NR == 1 || $2 == 1' good.pin > no_decoys.pin
awk -F '\t' 'NR == 1 || $2 == -1' good.pin > no_targets.pin
awk -F '\t' 'BEGIN { OFS = "\t" } NR == 3 { $6 = "abc" } { print }' \
  good.pin > text_feature.pin
awk -F '\t' 'BEGIN { OFS = "\t" } NR == 3 { $6 = "nan" } { print }' \
  good.pin > nan_feature.pin
awk -F '\t' 'BEGIN { OFS = "\t" } NR == 3 { print $1, $2, $3; next } { print }' \
  good.pin > short_row.pin
cp good.pin same.pin

cases=0
failures=0

# runs a command and checks its exit status and its stderr
expect() {  # name, expected status, pattern expected in stderr, command
  local name=$1 want=$2 pattern=$3
  shift 3
  "$@" > "$name.stdout.txt" 2> "$name.stderr.txt"
  local got=$?
  cases=$((cases + 1))
  local why=""
  if [ "$got" -ne "$want" ]; then
    why="exit status $got, expected $want"
  elif grep -q -E \
      "runtime error:|AddressSanitizer|LeakSanitizer|UndefinedBehaviorSanitizer" \
      "$name.stderr.txt"; then
    why="sanitizer report"
  elif [ -n "$pattern" ] && ! grep -q -e "$pattern" "$name.stderr.txt"; then
    why="no \"$pattern\" in stderr"
  fi
  if [ -n "$why" ]; then
    failures=$((failures + 1))
    echo "FAIL $name ($why)"
    tail -n 5 "$name.stderr.txt"
  else
    echo "ok   $name"
  fi
}

# a run that must fail: exit status 1 and the message
fails() {  # name, pattern, command
  local name=$1 pattern=$2
  shift 2
  expect "$name" 1 "$pattern" "$@"
}

# the control: the input of the cases below is good
expect good_input 0 "" "$PERC" -S 1 -m good_psms.txt good.pin
expect qv_good_input 0 "" "$QVAL" "$QD/target.xcorr" "$QD/null.xcorr"

# options
fails unknown_option "is invalid" "$PERC" --no-such-option good.pin
fails no_input "too few arguments" "$PERC"
fails num_folds_0 "num-folds option requires" "$PERC" --num-folds 0 good.pin
fails num_folds_text "num-folds option requires" \
  "$PERC" --num-folds abc good.pin
fails num_threads_0 "num-threads option requires" \
  "$PERC" --num-threads 0 good.pin
fails test_fdr_2 "testFDR option requires" "$PERC" -t 2 good.pin
fails subset_negative "subset-max-train option requires" \
  "$PERC" -N -5 good.pin
fails merge_normalization_unknown "merge-normalization has to be" \
  "$PERC" --merge-normalization nope good.pin
fails direction_check_unknown "direction-check has to be" \
  "$PERC" --direction-check nope good.pin
fails pooled_with_test "cannot be combined" \
  "$PERC" --direction-check pooled --merge-normalization test good.pin
fails search_input_unknown "search-input option has to be" \
  "$PERC" -I nope good.pin
fails default_direction_unknown "not found" "$PERC" -V nosuchfeature good.pin
fails parameter_file_missing "could not open" \
  "$PERC" --parameter-file no_such_file.txt good.pin
fails fasta_missing "Could not process the fasta database" \
  "$PERC" -f no_such_file.fasta -P decoy_ good.pin

# input
fails input_missing "Unable to read" "$PERC" no_such_file.pin
fails input_empty "Unable to read" "$PERC" empty.pin
fails input_header_only "too few features" "$PERC" header_only.pin
fails input_no_decoys "no decoy PSMs" "$PERC" no_decoys.pin
fails input_no_targets "no target PSMs" "$PERC" no_targets.pin
fails input_text_feature "error reading in feature vector" \
  "$PERC" text_feature.pin
fails input_nan_feature "strange feature" "$PERC" nan_feature.pin
fails input_short_row "error reading" "$PERC" short_row.pin
# mix-max with unequal numbers of target and decoy PSMs: a concatenated
# search taken for separate searches; --no-terminate runs it anyway
fails mix_max_unequal "needs one target and one decoy PSM per spectrum" \
  "$PERC" -I separate good.pin
fails mix_max_unequal_y "needs one target and one decoy PSM per spectrum" \
  "$PERC" -y good.pin
# RESET does not use mix-max and runs such an input
expect reset_unequal 0 "Running the Percolator-RESET algorithm" \
  "$PERC" -S 1 -I separate --reset-algorithm -r reset_unequal.txt good.pin
expect mix_max_unequal_no_terminate 0 "No-terminate flag set: ignoring error" \
  "$PERC" -S 1 -I separate --no-terminate -m unequal_psms.txt good.pin

# output
fails output_is_input "which is the input file" "$PERC" -m same.pin same.pin
if ! cmp -s good.pin same.pin; then
  failures=$((failures + 1))
  echo "FAIL output_is_input (the input file was changed)"
fi
fails output_twice "would both be written" \
  "$PERC" -m twice.txt -r twice.txt good.pin
fails output_missing_directory "Could not open the file" \
  "$PERC" -m no_such_directory/psms.txt good.pin
fails weights_missing_directory "Could not open the file" \
  "$PERC" -w no_such_directory/weights.txt good.pin

# qvality
fails qv_unknown_option "is invalid" \
  "$QVAL" --no-such-option "$QD/target.xcorr" "$QD/null.xcorr"
fails qv_input_missing "Could not open" "$QVAL" no_such_file.txt "$QD/null.xcorr"
fails qv_input_empty "contains no values" "$QVAL" empty.pin "$QD/null.xcorr"

# a full disk
if [ -c /dev/full ]; then
  for option in m M r B w J; do
    fails "full_disk_$option" "Could not write the file /dev/full" \
      "$PERC" -S 1 "-$option" /dev/full good.pin
  done
  fails full_disk_proteins "Could not write the file /dev/full" \
    "$PERC" -S 1 -f auto -P decoy_ -l /dev/full good.pin
  fails full_disk_stdout "Could not write the results to standard output" \
    bash -c '"$0" -S 1 good.pin > /dev/full' "$PERC"
  fails qv_full_disk_o "Could not write the file /dev/full" \
    "$QVAL" -o /dev/full "$QD/target.xcorr" "$QD/null.xcorr"
  fails qv_full_disk_stdout "Could not write the results to standard output" \
    bash -c '"$0" "$1" "$2" > /dev/full' "$QVAL" "$QD/target.xcorr" \
    "$QD/null.xcorr"
else
  echo "skip full_disk_* (no /dev/full)"
fi

echo "$cases cases, $failures failed"
[ "$failures" -eq 0 ]
