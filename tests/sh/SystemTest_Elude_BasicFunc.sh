#!/bin/bash
if [ "$1" != "" ]; then
cd $1
fi

# We are at percolator 'bin' folder now.

run_command=$(elude -t ../data/train_data/train_2.txt -e ../data/train_data/test_2.txt predictions.txt -s model.txt -i in_source_fragments.txt -g retention_index.txt -u)

if [ $? -eq 0 ]; then
  echo "$0 finished successfully!"
  exit 0
else
  echo "ERROR: $0 exited with an error code."
  exit 1
fi
