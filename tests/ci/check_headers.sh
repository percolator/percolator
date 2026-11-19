#!/usr/bin/env bash
# Checks that every header in src/ compiles on its own, so that it includes
# what it uses, and twice in one unit, so that it has an include guard.
# Crux and OpenMS compile these headers in their own order.
#
# usage: tests/ci/check_headers.sh <build-dir> [compiler and flags]
#   The build directory holds the headers that CMake generates (src/).
#   The compiler defaults to "c++ -fopenmp".
# Exits with status 1 if a header fails.
set -u
if [ $# -lt 1 ]; then
  sed -n '2,9p' "$0" >&2
  exit 2
fi
BUILD=$1
shift
CXX=("$@")
[ ${#CXX[@]} -eq 0 ] && CXX=(c++ -fopenmp)
SRC=$(cd "$(dirname "$0")/../.." && pwd)/src
GEN=$(cd "$BUILD" && pwd)/src

compile() {  # the lines of the unit
  printf '%s\n' "$@" | "${CXX[@]}" -std=c++17 -fsyntax-only -x c++ \
    -I"$SRC" -I"$SRC/picked_protein" -I"$GEN" - 2>&1
}

headers=0
failures=0
for path in "$SRC"/*.h; do
  h=${path##*/}
  headers=$((headers + 1))
  if ! out=$(compile "#include \"$h\""); then
    failures=$((failures + 1))
    echo "FAIL $h (alone)"
    printf '%s\n' "$out" | grep -m 4 'error'
  elif ! out=$(compile "#include \"$h\"" "#include \"$h\""); then
    failures=$((failures + 1))
    echo "FAIL $h (included twice)"
    printf '%s\n' "$out" | grep -m 4 'error'
  fi
done
echo "$headers headers, $failures failed"
[ "$failures" -eq 0 ]
