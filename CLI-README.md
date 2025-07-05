# Usage from CLI

The `orthoxml-tools` package also provides a command-line interface for working with OrthoXML files. After installation, you can access the CLI via:

```bash
orthoxml [options] <subcommand> [options]
```

## Subcommands

### **validate**
Validate an OrthoXML file against the schema version specified in the file itself.

```bash
orthoxml validate --infile path/to/file.orthoxml
```

**Options:**
- `--infile <file>`: Specify the input file (required).

**Example:**
```bash
orthoxml validate --infile examples/data/ex1.orthoxml
```

### **stats**
Display basic statistics.

```bash
orthoxml stats --infile path/to/file.orthoxml [--outfile <file>] 
```

**Options:**
- `--infile <file>`: Specify the input file (required).

**Example:**
```bash
orthoxml stats --infile examples/data/ex1.orthoxml
```

### **gene-stats**
Display statistics for gene count per taxon.

```bash
orthoxml gene-stats --infile path/to/file.orthoxml [--outfile <file>]
```

**Options:**
- `--infile <file>`: Specify the input file (required).
- `--outfile <file>`: Write stats to a CSV file.

**Example:**
```bash
orthoxml gene-stats --infile examples/data/ex1.orthoxml --outfile gene_stats.csv
```

### **filter**
Filter orthology groups based on a specified score and threshold.

```bash
orthoxml filter --infile path/to/file.orthoxml --score-name <name> --threshold <value> --strategy <topdown|bottomup>
```

**Options:**
- `--infile <file>`: Specify the input file (required).
- `--score-name <name>`: Specify the score to filter by (e.g., CompletenessScore).
- `--threshold <value>`: Set the threshold for filtering.
- `--strategy <topdown|bottomup>`: Choose the filtering strategy (default is `topdown`).
- `--outfile <file>`: Save output to a file. if not specified, the output will be printed to stdout.


**Examples:**
```bash
orthoxml filter --infile data/test_case_2.orthoxml --score-name CompletenessScore --threshold 0.9 --strategy topdown 
```

with file output:
```bash
orthoxml filter --infile data/test_case_2.orthoxml --score-name CompletenessScore --threshold 0.9 --strategy topdown --outfile filtered.orthoxml
```

### **taxonomy**
Print a human-readable taxonomy tree from the OrthoXML file.

```bash
orthoxml taxonomy --infile path/to/file.orthoxml
```

**Example:**
```bash
>>> orthoxml taxonomy --infile examples/data/ex3-int-taxon.orthoxml
Root
├── Mus musculus
└── Primates
    ├── Homo sapiens
    └── Pan troglodytes
```

### **export**
Export orthology data as pairs or groups.

```bash
orthoxml export <pairs|groups> --infile path/to/file.orthoxml [--outfile <file>]
```

**Options:**
- `--infile <file>`: Specify the input file (required).
- `--outfile <file>`: Save output to a file.

**Examples:**
```bash
orthoxml export pairs --infile examples/data/ex1-int-taxon.orthoxml  --outfile pairs.csv
orthoxml export groups --infile examples/data/ex1-int-taxon.orthoxml
```

### **split**
Split the tree into multiple trees based on rootHOGs.

```bash
orthoxml split --infile path/to/file.orthoxml --outdir path/to/output_folder
```

**Options:**
- `--infile <file>`: Specify the input OrthoXML file (required).
- `--outdir <folder>`: Specify the output folder where the trees will be saved.
- 
**Examples:**
```bash
orthoxml split --infile examples/data/ex4-int-taxon-multiple-rhogs.orthoxml --outdir tests_output/splits
```

## File Conversions

### OrthoXML to Newick Tree (NHX)
Convert OrthoXML to Newick (NHX) format.

```bash
orthoxml to-nhx --infile path/to/file.orthoxml --outdir path/to/output_folder --xref-tag [geneId,protId,...]    
```

**Options:**
- `--infile <file>`: Specify the input OrthoXML file (required).
- `--outdir <folder>`: Specify the output folder where the NHX files will be saved (required).
- `--xref-tag <tag>`: Specify the attribute of the `<gene>` element to use as the label for the leaves. Default is `protId`.

**Example:**
```bash
orthoxml to-nhx --infile examples/data/ex4-int-taxon-multiple-rhogs.orthoxml --outdir ./tests_output/trees --xref-tag geneId
```

### Newick Tree (NHX) to OrthoXML
Convert Newick (NHX) format to OrthoXML.

```bash
orthoxml from-nhx --infile path/to/file.nhx --outfile path/to/file.orthoxml
```

**Options:**
- `--infile <file>`: Specify the input nhx file or files. (at least one file is required).
  - You can specify multiple files by providing them as a space-separated list.
  - If you provide multiple files, they will be combined into a single OrthoXML output.
- `--outfile <folder>`: Specify the output OrthoXML file (required).

**Example:**
```bash
orthoxml from-nhx --infile examples/data/sample.nhx --outfile ./tests_output/from_nhx.orthoxml
orthoxml from-nhx --infile examples/data/sample2.nhx examples/data/sample.nhx --outfile ./tests_output/from_nhx21.orthoxml 
```



### **Help**
To see help for any command:

```bash
orthoxml --help
orthoxml -h
orthoxml stats --help
orthoxml stats -h
```
