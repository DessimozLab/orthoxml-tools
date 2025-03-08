# tree.py

from collections import defaultdict
from .loaders import load_orthoxml_file, parse_orthoxml
from .exceptions import OrthoXMLParsingError
from lxml import etree
from .models import Gene, Species, OrthologGroup, ParalogGroup, Taxon
from .exporters import get_ortho_pairs_recursive

class OrthoXMLTree:
    def __init__(
        self,
        genes: dict[str, Gene],
        species: list[Species],
        groups: list[OrthologGroup|ParalogGroup|Gene],
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

    def debug_repr(self):
        return f"OrthoXMLTree(genes={self.genes}, species={self.species}, groups={self.groups}, taxonomy={self.taxonomy}, orthoxml_version={self.orthoxml_version})"
        
    def __repr__(self):
        return f"OrthoXMLTree(genes=[{len(self.genes)} genes], species=[{len(self.species)} species], groups=[{len(self.groups)} groups], taxonomy=[{len(self.taxonomy)} taxons], orthoxml_version={self.orthoxml_version})"
    
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
            genes = defaultdict(Gene)
            for species in species_list:
                for gene in species.genes:
                    genes[gene._id] = gene

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
        """
        Create an OrthoXMLTree instance from an OrthoXML string.

        Args:
            xml_str: OrthoXML string

        Returns:
            OrthoXMLTree: Initialized OrthoXMLTree instance

        Raises:
            OrthoXMLParsingError: If there's an error parsing the string
        """
        try:
            xml_tree = etree.fromstring(xml_str)
            species_list, taxonomy, groups, orthoxml_version = parse_orthoxml(xml_tree)

            genes = defaultdict(Gene)
            for species in species_list:
                for gene in species.genes:
                    genes[gene._id] = gene
            
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


    def to_orthoxml(self, filepath=None, pretty=True, use_source_tree=False):
        """
        Serialize the OrthoXMLTree to an OrthoXML file.

        Args:
            filepath: Path to write the OrthoXML file
            pretty: Pretty-print the XML output (default: True)
            use_source_tree: Use the source XML tree to generate the output (default: False)

        Returns:
            str: OrthoXML file content
        """
        if use_source_tree:
            xml_tree = self.xml_tree
        else:
            raise NotImplementedError("Generating OrthoXML from scratch is not yet supported")
        
        if filepath:
            with open(filepath, "wb") as f:
                f.write(etree.tostring(xml_tree, pretty_print=pretty))
        else:
            return etree.tostring(xml_tree, pretty_print=pretty).decode()

    def to_ortho_pairs(self, filepath=None, sep=",") -> list[(str, str)]:
        """
        Recursively traverse the tree and return all of the
        ortholog pairs in the tree.
        Specify a filepath if you want to write the pairs to file.

        Args:
            filepath: Path to write the pairs to
        Returns:
            list[(str, str)]: List of ortholog pairs
        """
        pairs = []
        for ortho in self.groups:
            if isinstance(ortho, OrthologGroup):
                _, valid_pairs = get_ortho_pairs_recursive(ortho)
                pairs.extend(valid_pairs)
        
        if filepath:
            with open(filepath, "w") as f:
                f.writelines(f"{a}{sep}{b}\n" for a, b in pairs)

        return pairs
