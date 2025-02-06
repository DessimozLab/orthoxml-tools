# loaders.py
from importlib import resources
from lxml import etree
from .exceptions import OrthoXMLParsingError
from .logger import get_logger
import os

logger = get_logger(__name__)

def load_orthoxml_file(filepath: str, validate: bool = False) -> etree.ElementTree:
    """
    Load an OrthoXML file from disk.
    
    :param filepath: Path to the OrthoXML file.
    :return: An instance of the XML tree.
    """
    if not os.path.exists(filepath):
        raise OrthoXMLParsingError(f"OrthoXML file not found: {filepath}")
    
    try:
        tree = etree.parse(filepath)
    except Exception as e:
        raise OrthoXMLParsingError(f"Failed to load OrthoXML file: {e}")

    if validate:
        orthoxml_version = tree.getroot().attrib.get('version')
        
        if not validate_xml(tree, orthoxml_version):
            raise OrthoXMLParsingError(f"OrthoXML file is not valid for version {orthoxml_version}")
        else:
            logger.info(f"OrthoXML file is valid for version {orthoxml_version}")
            return tree
    else:
        return tree

def validate_xml(xml_tree, orthoxml_version):
    try:
        # Load XSD schema from package resources
        with resources.files('orthoxml.schemas').joinpath(f'orthoxml-{orthoxml_version}.xsd').open('rb') as schema_file:
            schema_root = etree.XML(schema_file.read())
            schema = etree.XMLSchema(schema_root)

        # Validate
        if schema.validate(xml_tree):
            return True
        else:
            logger.warning(schema.error_log)
            
    except Exception as e:
        logger.error(f"Error: {e}")

def parse_orthoxml(xml_tree) -> tuple:
    """
    Parse an OrthoXML document into genes, species, and groups.

    :param xml_tree: An instance of the XML tree.
    :return: A tuple of genes, species, groups, and the OrthoXML version.
    """
    orthoxml_version = xml_tree.getroot().attrib.get('version')
    
    genes = []
    species = []
    groups = []

    return genes, species, groups, orthoxml_version