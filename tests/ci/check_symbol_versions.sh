#!/usr/bin/env bash
# Checks that a Linux binary needs no newer versioned symbols from the
# shared libraries than a given system provides, e.g. RHEL 8: glibc 2.28,
# libstdc++ of GCC 8 (GLIBCXX_3.4.25, CXXABI_1.3.11), libgomp of GCC 8
# (GOMP_4.5).
#
# usage: tests/ci/check_symbol_versions.sh <binary> [<PREFIX_version> ...]
#   e.g. check_symbol_versions.sh percolator GLIBC_2.28 GLIBCXX_3.4.25
# Prints the libraries the binary needs and the newest version it requires
# of each prefix; exits with status 1 if one is newer than its limit.
set -u
if [ $# -lt 1 ]; then
  sed -n '2,10p' "$0" >&2
  exit 2
fi
BIN=$1
shift
objdump -p "$BIN" | grep NEEDED
versions=$(objdump -T "$BIN" | grep -oE '\b[A-Z]+_[0-9]+(\.[0-9]+)*\b' | sort -u)
status=0
for limit in "$@"; do
  prefix=${limit%_*}
  newest=$(printf '%s\n' "$versions" | grep -E "^${prefix}_[0-9]" | sort -V | tail -n 1)
  if [ -z "$newest" ]; then
    echo "$prefix: not needed"
    continue
  fi
  highest=$(printf '%s\n%s\n' "$newest" "$limit" | sort -V | tail -n 1)
  if [ "$highest" = "$limit" ]; then
    echo "$prefix: needs $newest (limit $limit)"
  else
    echo "$prefix: needs $newest, newer than the limit $limit"
    status=1
  fi
done
exit $status
