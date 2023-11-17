#!/bin/bash

# Ensure the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: ./script.sh <FASTA_FILE> <WINDOW_LENGTHS>"
    exit 1
fi

FASTA_FILE=$1
WINDOW_LENGTHS=$2  # This should be a space-separated list of window lengths

for WINDOW_LENGTH in $WINDOW_LENGTHS; do
  echo "calculating window length $WINDOW_LENGTH" 1>&2
  python mutation.py $FASTA_FILE -l $WINDOW_LENGTH
  echo ""
done
