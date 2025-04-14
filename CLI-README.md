# Usage from CLI

The `orthoxml-tools` package also provides a command-line interface for working with OrthoXML files. After installation, you can access the CLI via:

```bash
orthoxml <subcommand> [options]
```

## Subcommands

### **stats**
Display basic statistics and gene count per taxon.

```bash
orthoxml stats path/to/file.xml
```

**Options:**
- `--outfile <file>`: Write stats to a CSV file.
- `--validate`: Validate the OrthoXML file.
- `--completeness <threshold>`: Filter entries by CompletenessScore.

**Example:**
```bash
orthoxml stats data/orthoxml.xml --outfile stats.csv --validate
```

### **taxonomy**
Print a human-readable taxonomy tree from the OrthoXML file.

```bash
orthoxml taxonomy path/to/file.xml
```

**Options:**
- `--validate`
- `--completeness <threshold>`

**Example:**
```bash
orthoxml taxonomy data/orthoxml.xml --validate
```

### **export**
Export orthology data as pairs or groups.

```bash
orthoxml export <pairs|groups> path/to/file.xml
```

**Options:**
- `--outfile <file>`: Save output to a file.
- `--validate`
- `--completeness <threshold>`

**Examples:**
```bash
orthoxml export pairs data/orthoxml.xml --outfile pairs.csv
orthoxml export groups data/orthoxml.xml --validate
```

### **split**
Split the tree into multiple trees based on rootHOGs.

```bash
orthoxml split path/to/file.xml
```

**Options:**
- `--validate`
- `--completeness <threshold>`

**Example:**
```bash
orthoxml split data/orthoxml.xml
```

### **Help**
To see help for any command:

```bash
orthoxml --help
orthoxml stats --help
```
