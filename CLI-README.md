# Usage from CLI

The `orthoxml-tools` package also provides a command-line interface for working with OrthoXML files. After installation, you can access the CLI via:

```bash
orthoxml [options] <subcommand> [options]
```

**Global options:**
- `--validate`: Validate the OrthoXML file.

## Subcommands

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
orthoxml taxonomy --infile examples/data/ex1-int-taxon.orthoxml
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
orthoxml split --infile path/to/file.orthoxml
```


### **Help**
To see help for any command:

```bash
orthoxml --help
orthoxml -h
orthoxml stats --help
orthoxml stats -h
```
