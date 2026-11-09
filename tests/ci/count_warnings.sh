#!/usr/bin/env bash
# Counts the distinct compiler warnings in a build log (a warning in a header
# that several files include counts once) and compares the count with a
# budget, which should only go down.
#
# usage: tests/ci/count_warnings.sh <build-log> <budget>
# Exits with status 1 if there are more warnings than the budget.
set -u
if [ $# -ne 2 ]; then
  sed -n '2,7p' "$0" >&2
  exit 2
fi
LOG=$1
BUDGET=$2
warnings=$(grep -E '^[^ ]+:[0-9]+:[0-9]+: warning: ' "$LOG" | sort -u)
if [ -z "$warnings" ]; then
  n=0
else
  n=$(printf '%s\n' "$warnings" | wc -l | tr -d ' ')
  echo "Warnings by option:"
  printf '%s\n' "$warnings" | grep -oE '\[-W[^]]+\]$' | sort | uniq -c | sort -rn
fi
echo "$n distinct compiler warnings (budget $BUDGET)"
if [ "$n" -gt "$BUDGET" ]; then
  echo "More warnings than the budget; the new ones are among:"
  printf '%s\n' "$warnings"
  exit 1
fi
if [ "$n" -lt "$BUDGET" ]; then
  echo "Fewer warnings than the budget: lower WARNING_BUDGET in .github/workflows/build.yml to $n."
fi
