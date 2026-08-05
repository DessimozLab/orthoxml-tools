# orthoxml-tools

Tools for working with OrthoXML files.

## Table of contents

- [What is OrthoXML?](#what-is-orthoxml)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Common workflows](#common-workflows)
  - [Inspect an OrthoXML file](#inspect-an-orthoxml-file)
  - [Export data for downstream analysis](#export-data-for-downstream-analysis)
  - [Transform or subset data](#transform-or-subset-data)
  - [Convert between formats](#convert-between-formats)
- [Command reference](#command-reference)
- [Getting help](#getting-help)
- [Legacy API](#legacy-api)
- [Development and testing](#development-and-testing)

## What is OrthoXML?

> OrthoXML is a standard for sharing and exchanging orthology predictions. It provides a structured way to describe orthology relationships while preserving database-specific annotations. More details are available at [OrthoXML](https://github.com/qfo/orthoxml/tree/main).

## Installation

Install the package from PyPI:

```bash
pip install orthoxml-tools
```

Input OrthoXML files can be plain text or compressed with gzip (.gz) or bzip2 (.bz2).

## Quick start

The CLI follows this pattern:

```bash
orthoxml-tools [options] <subcommand> [options]
```

A few common first steps:

```bash
# Validate a file against the schema declared inside it
orthoxml-tools validate --infile examples/data/ex1.orthoxml

# Summarize the tree structure
orthoxml-tools stats --infile examples/data/ex3-int-taxon.orthoxml

# Print a human-readable taxonomy tree
orthoxml-tools taxonomy --infile examples/data/ex3-int-taxon.orthoxml
```

## Common workflows

### Inspect an OrthoXML file

Use these commands to understand the content of an OrthoXML file before processing it further.

```bash
# Basic statistics
orthoxml-tools stats --infile path/to/file.orthoxml

# Gene counts per taxon
orthoxml-tools gene-stats --infile path/to/file.orthoxml --outfile gene_stats.json

# Show the taxonomy tree in a human-readable format
orthoxml-tools taxonomy --infile path/to/file.orthoxml

# Write the taxonomy tree to NHX format with annotated internal nodes
orthoxml-tools taxonomy --infile path/to/file.orthoxml --outfile species-tree.nhx

# Export orthologous pairs
orthoxml-tools export-pairs ortho --infile path/to/file.orthoxml --outfile orthos.tsv
```

### Export data for downstream analysis

Export pairwise relationships or ortholog groups into tabular formats that are easy to analyze with other tools.

```bash
# Export orthologous pairs
orthoxml-tools export-pairs ortho \
  --infile examples/data/ex1-int-taxon.orthoxml \
  --outfile orthos.tsv

# Export paralogous pairs
orthoxml-tools export-pairs para \
  --infile examples/data/ex1-int-taxon.orthoxml \
  --outfile paras.tsv

# Export ortholog groups as a simple two-column table
orthoxml-tools export-ogs \
  --infile examples/data/sample-for-og.orthoxml \
  --outfile ogs.tsv \
  --id protId
```

### Transform or subset data

Use these commands to focus on a subset of the tree or to remove incomplete groups.

```bash
# Keep only selected species
orthoxml-tools subset \
  --infile examples/data/sample-for-subset.orthoxml \
  --outfile mammals.orthoxml \
  --species "Homo sapiens" "Mus musculus"

# Extract one or more HOGs as standalone root groups
orthoxml-tools subset \
  --infile examples/data/sample-for-subset.orthoxml \
  --outfile opistokonta.orthoxml \
  --hog-ids HOG_Opistokonta

# Filter by completeness score
orthoxml-tools filter \
  --infile examples/data/sample-for-filter.orthoxml \
  --threshold 0.24 \
  --strategy cascade-remove \
  --outfile filtered.orthoxml
```

### Convert between formats

The package also supports conversions to and from common phylogenetic and orthology formats.

```bash
# Convert OrthoXML to NHX trees
orthoxml-tools to-nhx \
  --infile examples/data/sample-for-nhx.orthoxml \
  --outdir ./tests_output/trees \
  --xref-tag protId \
  --encode-levels

# Convert NHX back to OrthoXML
orthoxml-tools from-nhx \
  --infile examples/data/sample.nhx \
  --outfile ./tests_output/from_nhx.orthoxml

# Convert OrthoFinder-style CSV to OrthoXML
orthoxml-tools from-csv \
  --infile examples/data/InputOrthogroups.csv \
  --outfile ./tests_output/orthofinder.orthoxml

# Convert OrthoXML back to OrthoFinder-style CSV
orthoxml-tools to-csv \
  --infile  ./tests_output/orthofinder.orthoxml \
  --outfile examples/data/InputOrthogroups.csv
```

## Command reference

### validate
Validate an OrthoXML file against the schema version declared in the file.

```bash
orthoxml-tools validate --infile path/to/file.orthoxml
```

Options:
- `--infile <file>`: Input OrthoXML file (required).

### stats
Display basic tree statistics.

```bash
orthoxml-tools stats --infile path/to/file.orthoxml
```

### gene-stats
Display gene counts per taxon.

```bash
orthoxml-tools gene-stats --infile path/to/file.orthoxml [--outfile <file>]
```

Options:
- `--infile <file>`: Input OrthoXML file (required).
- `--outfile <file>`: Write counts to a JSON file when provided.

### taxonomy
Print a taxonomy tree for the provided OrthoXML file.

```bash
orthoxml-tools taxonomy --infile path/to/file.orthoxml
```

Options:
- `--infile <file>`: Input OrthoXML file (required).
- `--outfile <file>`: When provided, write the taxonomy tree to this file in NHX format with internal nodes annotated.

### export-pairs
Export ortholog or paralog pairs as tab-separated output.

```bash
orthoxml-tools export-pairs <ortho|para> \
  --infile <file> \
  --outfile <file> \
  [--id <tag>] \
  [--chunk-size <number>] \
  [--buffer-size <bytes>]
```

Options:
- `--infile <file>`: Input OrthoXML file (required).
- `--outfile <file>`: Output file (required).
- `--id <tag>`: Identifier to use in the output (`id`, `geneId`, or `protId`).
- `--chunk-size <number>`: Number of pairs to buffer per write (default: 20,000).
- `--buffer-size <bytes>`: I/O buffer size in bytes (default: 4 MiB).

### export-ogs
Export orthologous groups as a simple TSV file.

```bash
orthoxml-tools export-ogs --infile path/to/file.orthoxml --outfile path/to/output.tsv [--id <tag>]
```

### subset
Extract a subset of an OrthoXML file by species and/or HOG IDs.

```bash
orthoxml-tools subset --infile path/to/file.orthoxml --outfile path/to/output.orthoxml \
  [--species SPECIES [SPECIES ...]] \
  [--species-file FILE] \
  [--hog-ids HOG_ID [HOG_ID ...]] \
  [--hog-ids-file FILE]
```

Options:
- `--infile <file>`: Input OrthoXML file (required).
- `--outfile <file>`: Output OrthoXML file (required).
- `--species <name> [<name> ...]`: One or more species names to keep.
- `--species-file <file>`: Plain-text file with one species name per line.
- `--hog-ids <id> [<id> ...]`: One or more HOG IDs to extract as new root groups.
- `--hog-ids-file <file>`: Plain-text file with one HOG ID per line.

Notes:
- HOG IDs can refer to any nesting level; the matched subtree is promoted to a root HOG in the output.
- If both a parent and child HOG ID are supplied, only the parent is extracted.

### split
Split the tree into multiple trees based on root HOGs.

```bash
orthoxml-tools split --infile path/to/file.orthoxml --outdir path/to/output_folder
```

### filter
Filter the tree by completeness score using a chosen strategy.

```bash
orthoxml-tools filter \
  --infile path/to/file.orthoxml \
  --threshold <value> \
  --strategy <cascade-remove|extract|reparent> \
  --outfile path/to/output.orthoxml
```

### to-nhx
Convert OrthoXML to Newick/NHX format.

```bash
orthoxml-tools to-nhx --infile path/to/file.orthoxml --outdir path/to/output_folder --xref-tag protId
```

Options:
- `--infile <file>`: Input OrthoXML file (required).
- `--outdir <folder>`: Output folder for generated files (required).
- `--xref-tag <tag>`: Gene attribute to use as the leaf label (default: `protId`).
- `--encode-levels`: Include group-level information as NHX comments.

### from-nhx
Convert Newick/NHX files back to OrthoXML.

```bash
orthoxml-tools from-nhx --infile path/to/file.nhx --outfile path/to/file.orthoxml [--species-encode nhx|underscore]
```

### from-csv
Convert an OrthoFinder-style CSV file to OrthoXML.

```bash
orthoxml-tools from-csv --infile path/to/file.csv --outfile path/to/file.orthoxml
```

> Note: Because the CSV format does not preserve full hierarchical structure, the resulting OrthoXML is reported at the root level and should be considered an exploratory conversion.

### to-csv
Export an OrthoXML file to the OrthoFinder-style TSV format.

```bash
orthoxml-tools to-csv --infile path/to/file.orthoxml --outfile path/to/file.tsv [--id <attribute>]
```

## Getting help

To see help for any command:

```bash
orthoxml-tools --help
orthoxml-tools -h
orthoxml-tools stats --help
orthoxml-tools stats -h
```

## Legacy API

The older object-oriented interface is deprecated and will be removed in v1.0.0. The legacy documentation remains in [LEGACY-README.md](LEGACY-README.md).

## Development and testing

```bash
uv install ".[test]"
pytest -vv

# CLI smoke test
tests/test_cli.sh
```
