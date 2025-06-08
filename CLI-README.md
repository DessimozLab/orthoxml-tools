# Usage from CLI

The `orthoxml-tools` package also provides a command-line interface for working with OrthoXML files. After installation, you can access the CLI via:

```bash
orthoxml FILE [options] <subcommand> [options]
```

**Global options:**
- `--validate`: Validate the OrthoXML file.

## Subcommands

### **stats**
Display basic statistics.

```bash
orthoxml path/to/file.xml stats 
```

**Options:**
- `--outfile <file>`: Write stats to a CSV file.

**Example:**
```bash
orthoxml examples/data/ex1.orthoxml --validate stats --outfile stats.csv
```

### **gene-stats**
Display statistics for gene count per taxon.

```bash
orthoxml path/to/file.xml gene-stats
```

**Options:**
- `--outfile <file>`: Write stats to a CSV file.

**Example:**
```bash
orthoxml examples/data/ex1.orthoxml gene-stats --outfile gene_stats.csv
```

### **filter**
Filter orthology groups based on a specified score and threshold.

```bash
orthoxml path/to/file.xml filter --score-name <name> --threshold <value> --strategy <topdown|bottomup>
```

**Options:**
- `--score-name <name>`: Specify the score to filter by (e.g., CompletenessScore).
- `--threshold <value>`: Set the threshold for filtering.
- `--strategy <topdown|bottomup>`: Choose the filtering strategy (default is `topdown`).
- `--outfile <file>`: Save output to a file. if not specified, the output will be printed to stdout.


**Examples:**
```bash
orthoxml data/test_case_2.orthoxml filter --score-name CompletenessScore --threshold 0.9 --strategy topdown 
```

with file output:
```bash
orthoxml data/test_case_2.orthoxml filter --score-name CompletenessScore --threshold 0.9 --strategy topdown --outfile filtered.orthoxml
```

### **taxonomy**
Print a human-readable taxonomy tree from the OrthoXML file.

```bash
orthoxml path/to/file.xml taxonomy
```

**Example:**
```bash
orthoxml examples/data/ex1-int-taxon.orthoxml --validate taxonomy
```

### **export**
Export orthology data as pairs or groups.

```bash
orthoxml path/to/file.xml export <pairs|groups> 
```

**Options:**
- `--outfile <file>`: Save output to a file.

**Examples:**
```bash
orthoxml examples/data/ex1-int-taxon.orthoxml export pairs --outfile pairs.csv
orthoxml examples/data/ex1-int-taxon.orthoxml --validate export groups
```

### **split**
Split the tree into multiple trees based on rootHOGs.

```bash
orthoxml split path/to/file.xml
```


### **Help**
To see help for any command:

```bash
orthoxml --help
orthoxml stats --help
```
