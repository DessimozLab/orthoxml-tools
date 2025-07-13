#!/bin/bash

set -e  # Exit immediately if a command fails
set -u  # Treat unset variables as errors

# Define test file paths
EXAMPLES_DIR="examples/data"
INFILE="$EXAMPLES_DIR/ex3-int-taxon.orthoxml"
FILTER_INFILE="$EXAMPLES_DIR/sample-for-filter.orthoxml"
VALIDATE_INFILE="$EXAMPLES_DIR/ex3.orthoxml"
MULTIPLE_RHOGS_INFILE="$EXAMPLES_DIR/ex4-int-taxon-multiple-rhogs.orthoxml"
OUT_GENE_STATS="tests_output/gene_stats.json"
OUT_EXPORT_PAIRS="tests_output/export_pairs.tsv"
OUT_FILTERED="tests_output/filtered.orthoxml"
OUT_EXPORT_OGS="tests_output/ogs.tsv"

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

echo -e "\n[5.1] Test: export ortho pairs with --outfile (default chunk & buffer sizes)"
orthoxml export-pairs ortho --infile "$INFILE" --outfile "$OUT_EXPORT_PAIRS"
cat "$OUT_EXPORT_PAIRS"

echo -e "\n[5.2] Test: export para pairs with --outfile (default chunk & buffer sizes)"
orthoxml export-pairs para --infile "$INFILE" --outfile "$OUT_EXPORT_PAIRS"
cat "$OUT_EXPORT_PAIRS"

echo -e "\n[5.3] Test: export pairs with --outfile (default chunk & buffer sizes) custom id"
orthoxml export-pairs ortho --infile "$INFILE" --outfile "$OUT_EXPORT_PAIRS" --id geneId
cat "$OUT_EXPORT_PAIRS"

echo -e "\n[5.4] Test: export pairs with custom --chunk-size and --buffer-size"
orthoxml export-pairs \
    ortho \
    --infile "$INFILE" \
    --outfile "$OUT_EXPORT_PAIRS" \
    --chunk-size 5000 \
    --buffer-size 1048576
cat "$OUT_EXPORT_PAIRS"

echo -e "\n[6] Test: export ortho groups with --outfile"
orthoxml export-ogs --infile "$INFILE" --outfile "$OUT_EXPORT_OGS"
cat "$OUT_EXPORT_OGS"

echo -e "\n[7] Test: split"
orthoxml split --infile "$MULTIPLE_RHOGS_INFILE" --outdir "tests_output/splits"

echo -e "\n[8.1] Test: filter cascade-remove"
orthoxml filter \
    --infile "$FILTER_INFILE" \
    --strategy cascade-remove \
    --threshold 0.24 \
    --outfile "$OUT_FILTERED"

#echo -e "\n[8.2] Test: filter reparent"
#orthoxml filter \
#    --infile "$FILTER_INFILE" \
#    --strategy  reparent \
#    --threshold 0.24 \
#    --outfile "$OUT_FILTERED"
#cat "$OUT_FILTERED"

echo -e "\n[8.3] Test: filter extract"
orthoxml filter \
    --infile "$FILTER_INFILE" \
    --strategy extract \
    --threshold 0.24 \
    --outfile "$OUT_FILTERED"
cat "$OUT_FILTERED"

echo -e "\n[9] Test: OrthoXML to NHX conversion"
orthoxml to-nhx \
    --infile "$MULTIPLE_RHOGS_INFILE" \
    --outdir "tests_output/trees" \
    --xref-tag geneId

echo -e "\n[10] Test: Newick (NHX) to OrthoXML conversion"
orthoxml from-nhx --infile "$EXAMPLES_DIR/sample.nhx" --outfile "tests_output/from_nhx.orthoxml"
orthoxml from-nhx --infile "$EXAMPLES_DIR/sample2.nhx" "$EXAMPLES_DIR/sample.nhx" --outfile "tests_output/from_nhx21.orthoxml"

echo -e "\n[11] Test: Orthofinder CSV to OrthoXML conversion"
orthoxml from-orthofinder --infile examples/data/OrthofinderOrthogroups.csv --outfile tests_output/orthofinder.orthoxml

echo -e "\n[12] Test: help commands"
orthoxml -h
orthoxml stats -h

echo -e "\n[13] Test: version"
orthoxml --version
orthoxml -v

echo -e "\n[14] Test: validation"
orthoxml validate --infile "$VALIDATE_INFILE"

echo -e "\nAll tests completed successfully."
