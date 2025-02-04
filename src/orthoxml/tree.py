# tree.py

from .loaders import load_orthoxml_file, parse_orthoxml
from .exceptions import OrthoXMLParsingError

class OrthoXMLTree:
    def __init__(self, genes, species, groups, xml_tree):
        self.genes = genes
        self.species = species
        self.groups = groups
        self.xml_tree = xml_tree

    def __repr__(self):
        return f"<OrthoXMLTree: {len(self.genes)} genes, {len(self.species)} species, {len(self.groups)} groups>"
    
    def __str__(self):
        return f"OrthoXMLTree: {len(self.species)} species, {len(self.groups)} groups"
    
    @classmethod
    def from_file(cls, filepath) -> "OrthoXMLTree":
        """
        Factory method to load an OrthoXMLDocument from a file.
        
        :param filepath: Path to the OrthoXML file.
        :return: An instance of OrthoXMLTree.
        """
        try:
            # Use the loader function to get the XML tree from a file.
            xml_tree = load_orthoxml_file(filepath)
        except Exception as e:
            raise OrthoXMLParsingError(f"Failed to load OrthoXML file: {e}")
        
        # Parse the XML tree into model objects.
        genes, species, groups = parse_orthoxml(xml_tree)

        return cls(genes, species, groups, xml_tree)

    @classmethod
    def from_string(cls, xml_str):
        pass

    def to_orthoxml(self, filepath=None, pretty=True):
        pass