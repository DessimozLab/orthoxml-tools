# loaders.py
from lxml import etree
from .exceptions import OrthoXMLParsingError


def load_orthoxml_file(filepath) -> etree.ElementTree:
    """
    Load an OrthoXML file from disk.
    
    :param filepath: Path to the OrthoXML file.
    :return: An instance of the XML tree.
    """
    try:
        return etree.parse(filepath)
    except Exception as e:
        raise OrthoXMLParsingError(f"Failed to load OrthoXML file: {e}")

def parse_orthoxml(xml_tree) -> tuple:
    """
    Parse an OrthoXML document into genes, species, and groups.

    :param xml_tree: An instance of the XML tree.
    :return: A tuple of genes, species and groups.
    """
    genes = []
    species = []
    groups = []

    return genes, species, groups