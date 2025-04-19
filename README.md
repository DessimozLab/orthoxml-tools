# orthoxml-tools

Tools for working with OrthoXML files.

## What is OrthoXML Format?

> OrthoXML is a standard for sharing and exchaning orthology predictions. OrthoXML is designed broadly to allow the storage and comparison of orthology data from any ortholog database. It establishes a structure for describing orthology relationships while still allowing flexibility for database-specific information to be encapsulated in the same format.  
> [OrthoXML](https://github.com/qfo/orthoxml/tree/main)

# Installation

```
pip install orthoxml
```

# Usage

```python
>>> from orthoxml import OrthoXMLTree
>>> otree = OrthoXMLTree.from_file("data/sample.orthoxml", validate=True)
>>> otree
2025-02-11 11:43:17 - loaders - INFO - OrthoXML file is valid for version 0.5
OrthoXMLTree(genes=[5 genes], species=[3 species], groups=[0 groups], taxonomy=[0 taxons], orthoxml_version=0.5)
```

### Filter Based on CompletenessScore at Loading
```python
>>> from orthoxml import OrthoXMLTree
>>> otree = OrthoXMLTree.from_file("data/sample.orthoxml", CompletenessScore_threshold=0.95, validate=True)
>>> otree
2025-02-11 11:43:17 - loaders - INFO - OrthoXML file is valid for version 0.5
OrthoXMLTree(genes=[5 genes], species=[3 species], groups=[0 groups], taxonomy=[0 taxons], orthoxml_version=0.5)
```

### Accessing Specific Data

*   **Groups**

```python
>>> otree.groups
OrthologGroup(taxonId=5, geneRefs=['5'], orthologGroups=[OrthologGroup(taxonId=4, geneRefs=['4'], orthologGroups=[], paralogGroups=[ParalogGroup(taxonId=None, geneRefs=['1', '2', '3'], orthologGroups=[], paralogGroups=[])])], paralogGroups=[])
```

*   **Genes**

```python
>>> otree.genes
defaultdict(orthoxml.models.Gene,
            {'1': Gene(id=1, geneId=hsa1, protId=None),
             '2': Gene(id=2, geneId=hsa2, protId=None),
             '3': Gene(id=3, geneId=hsa3, protId=None),
             '4': Gene(id=4, geneId=ptr1, protId=None),
             '5': Gene(id=5, geneId=mmu1, protId=None)})
```

*   **Taxonomy**

```python
>>> otree.taxonomy
Taxon(id=5, name=Root, children=[Taxon(id=3, name=Mus musculus, children=[]), Taxon(id=4, name=Primates, children=[Taxon(id=1, name=Homo sapiens, children=[]), Taxon(id=2, name=Pan troglodytes, children=[])])])
```

For a more human-readable tree structure:

```python
>>> print(otree.taxonomy.to_str())
Root
├── Mus musculus
└── Primates
    ├── Homo sapiens
    └── Pan troglodytes
```

*   **Species**

```python
>>> otree.species
[Species(name=Homo sapiens, NCBITaxId=9606, genes=[Gene(id=1, geneId=hsa1), Gene(id=2, geneId=hsa2), Gene(id=3, geneId=hsa3)]),
 Species(name=Pan troglodytes, NCBITaxId=9598, genes=[Gene(id=4, geneId=ptr1)]),
 Species(name=Mus musculus, NCBITaxId=10090, genes=[Gene(id=5, geneId=mmu1)])]
```

### Statistics of the OrthoXML tree

*   **Basic Stats**
```python
>>> otree.base_stats()
{'genes': 10,
 'species': 3,
 'groups': 3,
 'taxonomy': 0,
 'orthoxml_version': '0.5'}
```

*   **Gene Number per Taxonomic Level Stats**
```python
>>> otree.gene_stats()
{'5': 4, '3': 3, '4': 3, '2': 6, '1': 10}
>>> otree.gene_stats(filepath="out.csv", sep=",") # to also writes the stats to file with two columns: taxonId and gene_count
{'5': 4, '3': 3, '4': 3, '2': 6, '1': 10}
```

### Manipulate the Tree

* **Split an instance of OrthoXML Tree to separate OrthoXML Trees based on rootHOGs**
```python
>>> otrees = otree.split_by_rootHOGs()
>>> otrees[0].groups
OrthologGroup(taxonId=1, geneRefs=['1000000002'], orthologGroups=[OrthologGroup(taxonId=2, geneRefs=['1001000001', '1002000001'], orthologGroups=[], paralogGroups=[])], paralogGroups=[])
```

### Export Options

*   **Orthologous Pairs**

```python
>>> otree.to_ortho_pairs()
[('1', '2'), ('1', '3')]
>>> otree.to_ortho_pairs(filepath="out.csv") # to also writes the pairs to file
[('1', '2'), ('1', '3')]
```

*   **Get Orthologous Pairs of an Specific Gene**

```python
>>> otree.to_ortho_pairs_of_gene("1001000001")
[('1001000001', '1002000001'), ('1000000002', '1001000001')]
>>> otree.to_ortho_pairs_of_gene("1001000001", filepath="out.csv") # to also writes the pairs to file
[('1001000001', '1002000001'), ('1000000002', '1001000001')]
```

*   **Orthologous Groups**

```python
>>> otree.to_ogs()
[['1000000002', '1001000001', '1002000001'],
 ['1000000003', '1001000002', '1002000002'],
 ['1000000004', '1001000003', '1002000003']]
>>> otree.to_ogs(filepath="out.csv") # to also writes the groups to file
[['1000000002', '1001000001', '1002000001'],
 ['1000000003', '1001000002', '1002000002'],
 ['1000000004', '1001000003', '1002000003']]
```

# Usage from CLI
The `orthoxml-tools` package also provides a command-line interface for working with OrthoXML files. For more information on it checkout [CLI-README.md](./CLI-README.md)

## Testing

```
uv install `.[test]`
pytest -vv
```