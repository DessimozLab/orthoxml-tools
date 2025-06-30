#!/bin/bash

set -e  # Exit immediately if a command fails
set -u  # Treat unset variables as errors

# Define test file paths
EXAMPLES_DIR="examples/data"
INFILE="$EXAMPLES_DIR/ex3-int-taxon.orthoxml"
OUT_GENE_STATS="tests_output/gene_stats.json"
OUT_EXPORT_PAIRS="tests_output/export_pairs.tsv"
OUT_FILTERED="tests_output/filtered.orthoxml"

echo "Running orthoxml CLI tests..."

echo -e "\n[1] Test: stats"
orthoxml stats --infile "$INFILE"

echo -e "\n[2] Test: gene-stats"
orthoxml gene-stats --infile "$INFILE"

echo -e "\n[3] Test: gene-stats with --outfile"
orthoxml gene-stats --infile "$INFILE" --outfile "$OUT_GENE_STATS"
cat "$OUT_GENE_STATS"

echo -e "\n[4] Test: taxonomy"
orthoxml taxonomy --infile "$INFILE"

echo -e "\n[5] Test: export pairs"
orthoxml export pairs --infile "$INFILE"

echo -e "\n[6] Test: export pairs with --outfile"
orthoxml export pairs --infile "$INFILE" --outfile "$OUT_EXPORT_PAIRS"
cat "$OUT_EXPORT_PAIRS"

echo -e "\n[7] Test: split"
orthoxml split --infile "$INFILE"

echo -e "\n[8] Test: filter"
orthoxml filter --infile "$INFILE" --score-name CompletenessScore --threshold 0.5 --strategy topdown

echo -e "\n[9] Test: filter with --outfile"
orthoxml filter --infile "$INFILE" --score-name CompletenessScore --threshold 0.5 --strategy topdown --outfile "$OUT_FILTERED"
cat "$OUT_FILTERED"

echo -e "\n[10] Test: help commands"
orthoxml -h
orthoxml stats -h

echo -e "\nAll tests completed successfully."
