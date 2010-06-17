#!/bin/bash

cd $1

# We are at percolator 'bin' folder now.

# Performance tests.
# Using the "time list_of_executed_commands" syntax
# to measure overall execution time.
# Once the execution time is measured it should be
# parsed using regex and decision about performance
# gain/regression should be made.

time_before=$(date '+%s')
run_command=$(qvality 2>&1)
time_after=$(date '+%s')

# Checking time values.
if [ $((time_after - time_before)) -lt 5 ]; then
  echo "$0 finished successfully!"
  exit 0
else
  echo "ERROR: $0 exited with an error code."
  exit 1
fi


