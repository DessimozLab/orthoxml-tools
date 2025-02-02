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
from orthoxml import OrthoXMLTree

tree = OrthoXMLTree.from_file("data/orthoxml.xml")

print(tree)
```