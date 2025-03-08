# loaders.py
import os
from importlib import resources

from lxml import etree
from .exceptions import OrthoXMLParsingError
from .models import Gene, Species, OrthologGroup, ParalogGroup, Taxon, ORTHO_NS
from .logger import get_logger


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
            logger.info(f"OrthoXML file: {filepath} is valid for version {orthoxml_version}")
            return tree
    else:
        return tree

def validate_xml(xml_tree, orthoxml_version):
    """
    Validate an OrthoXML document against the scheme.
    
    :param xml_tree: An instance of the XML tree.
    :param orthoxml_version: The OrthoXML version.
    """
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

def parse_orthoxml(xml_tree) -> tuple[list[Species], Taxon, list[OrthologGroup|ParalogGroup|Gene], str]:
    """
    Parse an OrthoXML document into genes, species, and groups.

    :param xml_tree: An instance of the XML tree.
    :return: A tuple of species_list, taxonomy, groups, and the OrthoXML version.
    """    
    root = xml_tree.getroot()
    orthoxml_version = root.get("version", None)

    # Parse species.
    species_list = []
    for sp_el in root.xpath("ortho:species", namespaces={"ortho": ORTHO_NS}):
        species_list.append(Species.from_xml(sp_el))

    # Parse taxonomy if present.
    taxonomy = None
    taxonomy_el = root.find(f"{{{ORTHO_NS}}}taxonomy")
    if taxonomy_el is not None:
        taxon_el = taxonomy_el.find(f"{{{ORTHO_NS}}}taxon")
        if taxon_el is not None:
            taxonomy = Taxon.from_xml(taxon_el)

    # Parse groups.
    groups = []
    groups_el = root.find(f"{{{ORTHO_NS}}}groups")
    if groups_el is not None:
        # Iterate over all ortholog groups
        for ortholog_group_el in groups_el.findall(f"{{{ORTHO_NS}}}orthologGroup"):
            groups.append(OrthologGroup.from_xml(ortholog_group_el))

        # Iterate over all paralog groups
        for paralog_group_el in groups_el.findall(f"{{{ORTHO_NS}}}paralogGroup"):
            groups.append(ParalogGroup.from_xml(paralog_group_el))

        # Iterate over all gene references
        for gene_el in groups_el.findall(f"{{{ORTHO_NS}}}geneRef"):
            groups.append(Gene.from_xml(gene_el))

    return species_list, taxonomy, groups, orthoxml_version