#!/usr/bin/env bash
# Runs percolator and qvality from a build tree over a matrix of option
# combinations: the example input in data/, synthetic inputs made with
# tests/benchmark/make_synthetic_pin.py, and the qvality example.
#
# The sanitizer CI job runs it on an ASan+UBSan build, where every run must
# exit with status 0 and print no sanitizer report. check_determinism.sh
# runs it twice with different thread counts and compares the outputs.
#
# usage: tests/ci/option_matrix.sh <build-dir> <out-dir> [percolator options]
#   The percolator options (e.g. --num-threads 4) are added to every
#   percolator run. MATRIX_PSMS sets the number of PSMs of the synthetic
#   inputs (default 20000). PYTHON names the Python 3 interpreter.
#   The binaries are those of a single-configuration build tree
#   (<build-dir>/src/percolator) or of the Release configuration of a Visual
#   Studio one (<build-dir>/src/Release/percolator.exe; run the script with
#   Git Bash there). MATRIX_PERCOLATOR and MATRIX_QVALITY name other
#   programs to run, e.g. the commands of an installed wheel.
#   tests/ci/matrix_hashes.py lists the hashes of the outputs, for comparing
#   platforms.
# Exits with status 1 if any run fails.
set -u
if [ $# -lt 2 ]; then
  sed -n '2,22p' "$0" >&2
  exit 2
fi
BUILD=$1
OUT=$2
shift 2
EXTRA=("$@")
SRC=$(cd "$(dirname "$0")/../.." && pwd)
PERC=$BUILD/src/percolator
QVAL=$BUILD/src/qvality/qvality
if [ ! -e "$PERC" ] && [ -e "$BUILD/src/Release/percolator.exe" ]; then
  PERC=$BUILD/src/Release/percolator.exe
  QVAL=$BUILD/src/qvality/Release/qvality.exe
fi
PERC=${MATRIX_PERCOLATOR:-$PERC}
QVAL=${MATRIX_QVALITY:-$QVAL}
EX=$SRC/data/percolator/tab/percolatorTab
QD=$SRC/data/qvality
NPSMS=${MATRIX_PSMS:-20000}
PY=${PYTHON:-python3}

for exe in "$PERC" "$QVAL"; do
  if ! command -v "$exe" > /dev/null; then
    echo "error: $exe not found" >&2
    exit 2
  fi
done

IN=$OUT/inputs
mkdir -p "$IN"
"$PY" "$SRC/tests/benchmark/make_synthetic_pin.py" -n "$NPSMS" -f 20 \
  --mode concatenated --seed 1 -o "$IN/concatenated.pin" || exit 2
"$PY" "$SRC/tests/benchmark/make_synthetic_pin.py" -n "$NPSMS" -f 20 \
  --mode separate --seed 2 -o "$IN/separate.pin" || exit 2
# the concatenated input split into two files, for the multi-file path
head -n 1 "$IN/concatenated.pin" > "$IN/part1.pin"
head -n 1 "$IN/concatenated.pin" > "$IN/part2.pin"
awk 'NR > 1 && NR % 2 == 0' "$IN/concatenated.pin" >> "$IN/part1.pin"
awk 'NR > 1 && NR % 2 == 1' "$IN/concatenated.pin" >> "$IN/part2.pin"
# a FASTA file with the proteins of the concatenated input
awk -f "$SRC/tests/ci/proteins_fasta.awk" "$IN/concatenated.pin" \
  > "$IN/concatenated.fasta" || exit 2

runs=0
failures=0

# checks the exit status and stderr of a finished run
check() {  # name, exit status, stderr file
  runs=$((runs + 1))
  if [ "$2" -ne 0 ] || grep -q -E \
      "runtime error:|AddressSanitizer|LeakSanitizer|UndefinedBehaviorSanitizer" \
      "$3"; then
    failures=$((failures + 1))
    echo "FAIL $1 (exit status $2)"
    tail -n 30 "$3"
  else
    echo "ok   $1"
  fi
}

# runs percolator with all outputs written into $OUT/<name>/
perc() {  # name, then options and input files
  local name=$1
  shift
  local o=$OUT/$name
  mkdir -p "$o"
  "$PERC" -S 1 ${EXTRA[@]+"${EXTRA[@]}"} -m "$o/psms.txt" -M "$o/decoy_psms.txt" \
    -r "$o/peptides.txt" -B "$o/decoy_peptides.txt" -w "$o/weights.txt" \
    "$@" > "$o/stdout.txt" 2> "$o/stderr.txt"
  check "$name" $? "$o/stderr.txt"
}

qval() {  # name, then options and input files
  local name=$1
  shift
  local o=$OUT/$name
  mkdir -p "$o"
  "$QVAL" "$@" > "$o/stdout.txt" 2> "$o/stderr.txt"
  check "$name" $? "$o/stderr.txt"
}

perc ex_default "$EX"
perc ex_only_psms -U "$EX"
perc ex_tdc -Y "$EX"
perc ex_tdc_pep --tdc-pep "$EX"
perc ex_classic --classic-training "$EX"
perc ex_quick_nested -x --nested-xval-bins 2 "$EX"
perc ex_nested_best_positive --nested-xval-bins 3 --train-best-positive "$EX"
perc ex_subset -N 5000 "$EX"
perc ex_reset --reset-algorithm "$EX"
perc ex_reset_tdc --reset-algorithm -Y "$EX"
perc ex_composition --reset-algorithm --composition-match "$EX"
perc ex_irls_pep --irls-pep "$EX"
perc ex_pava_pep --pava-pep "$EX"
perc ex_rank_pep --rank-pep "$EX"
perc ex_unitnorm -u "$EX"
perc ex_each_iteration -R "$EX"
perc ex_direction -V -deltCn "$EX"
perc ex_tab_out -J "$OUT/ex_tab_out/features.tab" "$EX"
perc ex_init_weights -W "$OUT/ex_default/weights.txt" "$EX"
perc ex_static_weights -W "$OUT/ex_default/weights.txt" --static "$EX"
perc ex_protein -f auto -P decoy_ -l "$OUT/ex_protein/proteins.txt" \
  -L "$OUT/ex_protein/decoy_proteins.txt" "$EX"
perc ex_protein_cg -f auto -P decoy_ -c -g \
  -l "$OUT/ex_protein_cg/proteins.txt" "$EX"
perc ex_spectral -f auto -P decoy_ --spectral-counting-fdr 0.01 \
  -l "$OUT/ex_spectral/proteins.txt" "$EX"
perc syn_concatenated "$IN/concatenated.pin"
perc syn_concatenated_subset -N $((NPSMS / 4)) "$IN/concatenated.pin"
perc syn_concatenated_tdc_pep --tdc-pep "$IN/concatenated.pin"
perc syn_separate "$IN/separate.pin"
perc syn_separate_tdc -Y "$IN/separate.pin"
perc syn_separate_tdc_pava --tdc-pep --pava-pep "$IN/separate.pin"
perc syn_multi_file "$IN/part1.pin" "$IN/part2.pin"
perc syn_folds1 --num-folds 1 "$IN/concatenated.pin"
perc syn_folds5 --num-folds 5 "$IN/concatenated.pin"
perc syn_folds10_separate --num-folds 10 "$IN/separate.pin"
perc syn_in_place --in-place-training "$IN/concatenated.pin"
perc syn_in_place_subset --in-place-training -N $((NPSMS / 4)) \
  "$IN/separate.pin"
perc syn_merge_train --merge-normalization train "$IN/concatenated.pin"
perc syn_merge_train_decoys --merge-normalization train-decoys \
  "$IN/separate.pin"
perc syn_merge_none --merge-normalization none "$IN/concatenated.pin"
perc syn_direction_pooled --direction-check pooled "$IN/concatenated.pin"
perc syn_direction_second_chance --direction-check second-chance \
  "$IN/separate.pin"
perc syn_search_input_separate -I separate "$IN/concatenated.pin"
perc syn_search_input_concatenated -I concatenated "$IN/separate.pin"
perc syn_train_fdr_initial --train-fdr-initial 0.05 "$IN/concatenated.pin"
perc syn_fasta -f "$IN/concatenated.fasta" -z trypsin -P decoy_ \
  -l "$OUT/syn_fasta/proteins.txt" -L "$OUT/syn_fasta/decoy_proteins.txt" \
  "$IN/concatenated.pin"
perc syn_stdin --stdinput-tab < "$IN/concatenated.pin"
perc syn_tab_in -j "$IN/concatenated.pin"
qval qv_default "$QD/target.xcorr" "$QD/null.xcorr"
qval qv_irls_pep --irls-pep "$QD/target.xcorr" "$QD/null.xcorr"
qval qv_pava_pep --pava-pep "$QD/target.xcorr" "$QD/null.xcorr"
qval qv_tdc_pep --tdc-pep -d "$QD/target.xcorr" "$QD/null.xcorr"
qval qv_tdc -Y "$QD/target.xcorr" "$QD/null.xcorr"

echo "$runs runs, $failures failed"
[ "$failures" -eq 0 ]
