# tree.py

from .loaders import load_orthoxml_file, parse_orthoxml
from .exceptions import OrthoXMLParsingError
from lxml import etree
from .models import Gene, Species, OrthologyGroup

class OrthoXMLTree:
    def __init__(
        self,
        genes: dict[str, Gene],
        species: dict[str, Species],
        groups: list[OrthologyGroup],
        xml_tree: etree.ElementTree,
        orthoxml_version: str = None
    ):
        self.genes = genes
        self.species = species
        self.groups = groups
        self.xml_tree = xml_tree
        self.orthoxml_version = orthoxml_version

    def __repr__(self):
        return f"<OrthoXMLTree: {len(self.genes)} genes, {len(self.species)} species, {len(self.groups)} groups>"
    
    def __str__(self):
        return f"OrthoXMLTree: {len(self.genes)} genes, {len(self.species)} species, {len(self.groups)} groups"
    
    @classmethod
    def from_file(
        cls, 
        filepath: str, 
        validate: bool = False,
    ) -> "OrthoXMLTree":
        """
        Create an OrthoXMLTree instance from an OrthoXML file.

        Args:
            filepath: Path to the OrthoXML file
            orthoxml_version: OrthoXML schema version to use (default: None)

        Returns:
            OrthoXMLTree: Initialized OrthoXMLTree instance

        Raises:
            OrthoXMLParsingError: If there's an error loading or parsing the file
        """
        try:
            # Load XML document and validate against schema
            xml_tree = load_orthoxml_file(filepath, validate)
            
            # Parse XML elements into domain models
            genes, species, groups, orthoxml_version = parse_orthoxml(xml_tree)

            return cls(
                genes=genes,
                species=species,
                groups=groups,
                xml_tree=xml_tree,
                orthoxml_version=orthoxml_version
            )

        except etree.XMLSyntaxError as e:
            raise OrthoXMLParsingError(f"Invalid XML syntax: {str(e)}") from e
        except Exception as e:
            raise OrthoXMLParsingError(f"Error parsing OrthoXML: {str(e)}") from e

    @classmethod
    def from_string(cls, xml_str):
        pass

    def to_orthoxml(self, filepath=None, pretty=True):
        pass