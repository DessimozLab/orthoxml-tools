# loaders.py

import os
from importlib import resources
from typing import Union

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

def validate_xml(xml_tree, orthoxml_version) -> bool:
    """
    Validate an OrthoXML document against the scheme.
    
    :param xml_tree: An instance of the XML tree.
    :param orthoxml_version: The OrthoXML version.

    :return: True if the document is valid, False otherwise.
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
            return False
            
    except Exception as e:
        logger.error(f"Error: {e}")

def parse_orthoxml(xml_tree) -> tuple[list[Species], Taxon, list[Union[OrthologGroup, ParalogGroup, Gene]], str]:
    """
    Parse an OrthoXML document into species, taxonomy, groups, and version.

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

def filter_by_score(xml_tree, score_id, score_threshold, skip_no_scores=True, keep_low_score_parents=False) -> None:
    """
    Filter OrthoXML document by score. Works in-place.

    Code from https://github.com/DessimozLab/FastOMA/blob/main/FastOMA/zoo/hog/filter_orthoxml.py

    :param xml_tree: An instance of the XML tree.
    :param score_id: The score ID.
    :param score_threshold: The score threshold.
    :param skip_no_scores: If True, skip ortholog groups without scores. If False, remove them.
    :param keep_low_score_parents: If True, keep parents of low-scoring ortholog groups.
    """
    root = xml_tree.getroot()
    to_rem = []
    for hog in root.iterfind('.//{{{0}}}orthologGroup'.format(ORTHO_NS)):
        score = hog.find('./{{{0}}}score'.format(ORTHO_NS))
        if score is None:
            if skip_no_scores:
                continue
            else:
                to_rem.append(hog)
                continue
        if score.get('id') == score_id and float(score.get('value')) < score_threshold:
            to_rem.append(hog)
        if keep_low_score_parents and score.get('id') == score_id and float(score.get('value')) >= score_threshold:
            # remove all parents of this score from remove list
            for parent in hog.iterancestors():
                if parent.tag == "{{{0}}}orthologGroup".format(ORTHO_NS):
                    if parent in to_rem:
                        to_rem.remove(parent)

    for h in to_rem:
        parent = h.getparent()
        if parent is not None:
            parent.remove(h)
            if sum(c.tag == "{{{0}}}orthologGroup".format(ORTHO_NS) for c in parent) == 0:
                to_rem.append(parent)
