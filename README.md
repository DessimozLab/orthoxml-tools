# orthoxml-tools

Tools for working with OrthoXML files.

## What is OrthoXML Format?

> OrthoXML is a standard for sharing and exchanging orthology predictions. OrthoXML is designed broadly to allow the storage and comparison of orthology data from any ortholog database. It establishes a structure for describing orthology relationships while still allowing flexibility for database-specific information to be encapsulated in the same format.  
> [OrthoXML](https://github.com/qfo/orthoxml/tree/main)

# Installation

```
pip install orthoxml-tools
```

# Usage

```bash
orthoxml-tools [options] <subcommand> [options]
```

> Note: Input OrthoXML files can be in plain text or compressed format. Both gzip (.gz) and bzip2 (.bz2) compression are supported.

## Subcommands

### 🛠️ **validate**
Validate an OrthoXML file against the schema version specified in the file itself.

```bash
orthoxml-tools validate --infile path/to/file.orthoxml
```

**Options:**
- `--infile <file>`: Specify the input file (required).

**Example:**
```bash
orthoxml-tools validate --infile examples/data/ex1.orthoxml
```

### 🛠️ **stats**
Display basic statistics.

```bash
orthoxml-tools stats --infile path/to/file.orthoxml [--outfile <file>] 
```

**Options:**
- `--infile <file>`: Specify the input file (required).

**Example:**
```bash
orthoxml-tools stats --infile examples/data/ex3-int-taxon.orthoxml
```

### 🛠️ **gene-stats**
Display statistics for gene count per taxon.

```bash
orthoxml-tools gene-stats --infile path/to/file.orthoxml [--outfile <file>]
```

**Options:**
- `--infile <file>`: Specify the input file (required).
- `--outfile <file>`: Write stats to a txt file.

**Example:**
```bash
orthoxml-tools gene-stats --infile examples/data/ex3-int-taxon.orthoxml --outfile gene_stats.txt
```

### 🛠️ **filter**
Filter orthology groups based on CompletenessScore score and a threshold and strategy.

```bash
orthoxml-tools filter --infile path/to/file.orthoxml --threshold <value> --strategy <cascade-remove|extract|reparent> --outfile <file>
```

**Options:**
- `--infile <file>`: Specify the input file. (required)
- `--threshold <value>`: Set the threshold for filtering. value below this will be removed. (required)
- `--strategy <cascade-remove|extract|reparent>`: Choose the filtering strategy (default is `cascade-remove`).
- `--outfile <file>`: Save output to a file. if not specified, the output will be printed to stdout. (required)


**Examples:**
```bash
 orthoxml-tools filter --infile examples/data/sample-for-filter.orthoxml --score-name CompletenessScore --strategy top-down --threshold 0.24 --outfile tests_output/filtered_stream.orthoxml
```

### 🛠️ **taxonomy**
Print a human-readable taxonomy tree from the OrthoXML file.

```bash
orthoxml-tools taxonomy --infile path/to/file.orthoxml
```

**Example:**
```bash
>>> orthoxml-tools taxonomy --infile examples/data/ex3-int-taxon.orthoxml
Root
├── Mus musculus
└── Primates
    ├── Homo sapiens
    └── Pan troglodytes
```

### 🛠️ **export-pairs**
Export pairs (orthologs or paralogs) in TSV form, with configurable chunking and buffering.

```bash
orthoxml-tools export-pairs <ortho|para> \
    --infile <file> \
    --outfile <file> \
    [--id <tag>] \
    [--chunk-size <number>] \
    [--buffer-size <bytes>]
```

**Positional arguments:**
<ortho|para>
Choose which pair type to export:
- `ortho`: orthologous pairs
- `para`: paralogous pairs

**Options:**
- `--infile <file>`: Input OrthoXML file (required).
- `--outfile <file>`: Write output CSV to this file (required).
- `--id <tag>`: Gene attribute to use as identifier (default: id).
- `--chunk-size <number>`: Number of pairs to process per chunk (default: 20_000).
- `--buffer-size <bytes>`: I/O buffer size in bytes (default: 4194304).

**Examples:**

```bash
# [5.1] Export ortholog pairs with default chunk & buffer sizes
orthoxml-tools export-pairs ortho \
    --infile examples/data/ex1-int-taxon.orthoxml \
    --outfile orthos.csv

# [5.2] Export paralog pairs with default chunk & buffer sizes
orthoxml-tools export-pairs para \
    --infile examples/data/ex1-int-taxon.orthoxml \
    --outfile paras.csv

# [5.3] Export ortholog pairs using `geneId` as the identifier column
orthoxml-tools export-pairs ortho \
    --infile examples/data/ex1-int-taxon.orthoxml \
    --outfile orthos_geneid.csv \
    --id geneId

# [5.4] Export ortholog pairs with custom chunk and buffer sizes
orthoxml-tools export-pairs ortho \
    --infile examples/data/ex1-int-taxon.orthoxml \
    --outfile orthos_custom.csv \
    --chunk-size 5000 \
    --buffer-size 1048576
```


### 🛠️ **export-ogs**
Export Orthologous Groups as TSV file.

```bash
orthoxml-tools export-ogs --infile path/to/file.orthoxml --outfile path/to/output.tsv [--id <tag>]
```

**Options:**
- `--infile <file>`: Input OrthoXML file (required).
- `--outfile <file>`: Write output CSV to this file (required).
- `--id <tag>`: Gene attribute to use as identifier (default: id).

**Examples:**
```bash
orthoxml-tools export-ogs --infile examples/data/sample-for-og.orthoxml --outfile tests_output/ogs.tsv --id protId
```

### 🛠️ **split**
Split the tree into multiple trees based on rootHOGs.

```bash
orthoxml-tools split --infile path/to/file.orthoxml --outdir path/to/output_folder
```

**Options:**
- `--infile <file>`: Specify the input OrthoXML file (required).
- `--outdir <folder>`: Specify the output folder where the trees will be saved.
- 
**Examples:**
```bash
orthoxml-tools split --infile examples/data/ex4-int-taxon-multiple-rhogs.orthoxml --outdir tests_output/splits
```

## File Conversions

### 🛠️ **OrthoXML to Newick Tree (NHX)**
Convert OrthoXML to Newick (NHX) format.

```bash
orthoxml-tools to-nhx --infile path/to/file.orthoxml --outdir path/to/output_folder --xref-tag [geneId,protId,...]    
```

**Options:**
- `--infile <file>`: Specify the input OrthoXML file (required).
- `--outdir <folder>`: Specify the output folder where the NHX files will be saved (required).
- `--xref-tag <tag>`: Specify the attribute of the `<gene>` element to use as the label for the leaves. Default is `protId`.
- `--encode-levels`: If set, encode group levels as NHX comments in the output tree. This is useful for visualizing the hierarchy of orthologous groups.
  
**Example:**
```bash
orthoxml-tools to-nhx --infile examples/data/sample-for-nhx.orthoxml --outdir ./tests_output/trees --xref-tag protId --encode-levels
```

### 🛠️ **Newick Tree (NHX) to OrthoXML**
Convert Newick (NHX) format to OrthoXML.

```bash
orthoxml-tools from-nhx --infile path/to/file.nhx --outfile path/to/file.orthoxml
```

**Options:**
- `--infile <file>`: Specify the input nhx file or files. (at least one file is required).
  - You can specify multiple files by providing them as a space-separated list.
  - If you provide multiple files, they will be combined into a single OrthoXML output.
- `--outfile <folder>`: Specify the output OrthoXML file (required).

**Example:**
```bash
orthoxml-tools from-nhx --infile examples/data/sample.nhx --outfile ./tests_output/from_nhx.orthoxml
orthoxml-tools from-nhx --infile examples/data/sample2.nhx examples/data/sample.nhx --outfile ./tests_output/from_nhx21.orthoxml 
```

### 🛠️ CSV to OrthoXML (exploratory feature)
Convert a CSV file to OrthoXML. The CSV file is structured such that each row represents an orthogroup (OG), each column corresponds to a species, and each cell contains a gene name. This format is generated by OrthoFinder e.g. `examples/data/InputOrthogroups.csv`.

> [!WARNING]
> Note that since the CSV does not contain the full information required to represent the hierarchical structure of HOGs, the output OrthoXML file is reported at the root level. It should not be considered a full-fledged OrthoXML file.

```bash
orthoxml-tools from-csv --infile path/to/file.csv --outfile path/to/file.orthoxml
```

**Options:**
- `--infile <file>`: Specify the input orthogroups.csv file (required).
- `--outfile <folder>`: Specify the output OrthoXML file (required).

**Example:**
```bash
orthoxml-tools from-csv --infile examples/data/InputOrthogroups.csv --outfile tests_output/orthofinder.orthoxml
```


### 🛠️ **filter**
Filter the OrthoXML tree by a completeness score. 

- `--score-name <str>`: Name of the field for completeness score annotation (e.g. 'CompletenessScore') 
- `--threshold <float>`: Threshold value for the completeness score
- `--strategy <bottomup|topdown>`: Filtering strategy. Bottom-up will keep complete subHOGs even if they parents are incomplete.
- `--outfile <file>`: If provided, write the filtered OrthoXML to this file; otherwise, print to stdout

```bash
orthoxml-tools tests/test-data/case_filtering.orthoxml filter --score-name CompletenessScore \
                                                        --threshold 0.75 \
                                                        --strategy bottomup \
                                                        --outfile output-oxml.orthoxml 
```

### **Help**
To see help for any command:

```bash
orthoxml-tools --help
orthoxml-tools -h
orthoxml-tools stats --help
orthoxml-tools stats -h
```

## Legacy API

The `orthoxml-tools` package used to provides a object oriented interface for working with OrthoXML files. This API is deprecated and will be removed in v1.0.0. Please use the new streaming CLI method. The documentation on it can be found [here](LEGACY-README.md).

## Testing

```
uv install `.[test]`
pytest -vv

# test cli
tests/test_cli.sh
```
