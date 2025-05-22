# Usage from CLI

The `orthoxml-tools` package also provides a command-line interface for working with OrthoXML files. After installation, you can access the CLI via:

```bash
orthoxml FILE [options] <subcommand> [options]
```

**Global options:**
- `--validate`: Validate the OrthoXML file.
- `--completeness <threshold>`: Filter entries by CompletenessScore.
  
## Subcommands

### **stats**
Display basic statistics and gene count per taxon.

```bash
orthoxml path/to/file.xml stats 
```

**Options:**
- `--outfile <file>`: Write stats to a CSV file.

**Example:**
```bash
orthoxml examples/data/ex1.orthoxml --validate stats --outfile stats.csv
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
