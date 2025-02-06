# tree.py

from .loaders import load_orthoxml_file, parse_orthoxml
from .exceptions import OrthoXMLParsingError
from lxml import etree
from .models import Gene, Species, OrthologGroup, ParalogGroup, Taxon

class OrthoXMLTree:
    def __init__(
        self,
        genes: list[Gene],
        species: list[Species],
        groups: list[OrthologGroup|ParalogGroup],
        taxonomy: Taxon,
        xml_tree: etree.ElementTree,
        orthoxml_version: str = None
    ):
        self.genes = genes
        self.species = species
        self.groups = groups
        self.taxonomy = taxonomy
        self.xml_tree = xml_tree
        self.orthoxml_version = orthoxml_version

    def __repr__(self):
        return f"OrthoXMLTree(genes={self.genes}, species={self.species}, groups={self.groups}), taxonomy={self.taxonomy}, orthoxml_version={self.orthoxml_version}"
        
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
            species_list, taxonomy, groups, orthoxml_version = parse_orthoxml(xml_tree)

            # TODO: Parse genes one time and avoid duplicate representations
            genes = []
            for species in species_list:
                for gene in species.genes:
                    genes.append(gene)

            return cls(
                genes=genes,
                species=species_list,
                groups=groups,
                taxonomy=taxonomy,
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