from pathlib import Path
from io import BytesIO
import xml.etree.ElementTree as ET
import pytest
import types
import dendropy

from orthoxml.converters.from_nhx import (
    orthoxml_from_newicktrees,
    nhx_species_encoded_leaf,
    label_with_species_end,
    nhx_to_event,
    OrthoXMLBuilder,
    GeneRefHelper
)



# -------------------------
# Label parsing tests
# -------------------------

def make_tree_from_str(newick_str):
    return dendropy.Tree.get(data=newick_str, schema="newick", preserve_underscores=True, extract_comment_metadata=True)


def test_label_with_species_end():
    tree = make_tree_from_str("(P53_HUMAN,MYC_MOUSE);")
    leaf = tree.find_node_with_taxon_label("P53_HUMAN")
    label, species = label_with_species_end(leaf)
    assert label == "P53_HUMAN"
    assert species == "HUMAN"


def test_nhx_species_encoded_leaf():
    tree = make_tree_from_str("(P53:0.32[&&NHX:S=Homo sapiens],MYC:0.2[&&NHX:S=Mus musculus]):0.1;")
    leaf = tree.find_node_with_taxon_label("P53")
    gene, species = nhx_species_encoded_leaf(leaf)
    assert gene == "P53"
    assert species == "Homo sapiens"


def test_nhx_species_encoded_leaf_real_node_T():
    newick = "(BRCA1[&&NHX:T=Pan troglodytes],MYC[&&NHX:S=Mus_musculus]);"
    tree = make_tree_from_str(newick)
    leaf = tree.find_node_with_taxon_label("BRCA1")
    label, species = nhx_species_encoded_leaf(leaf)
    assert label == "BRCA1"
    assert species == "Pan troglodytes"


def test_nhx_species_encoded_leaf_missing_annotations():
    tree = make_tree_from_str("(NO_ANNOT,FOO);")
    leaf = tree.find_node_with_taxon_label("NO_ANNOT")
    with pytest.raises(ValueError, match="no NHX annotations found"):
        nhx_species_encoded_leaf(leaf)


def test_species_end_without_underscore():
    tree = make_tree_from_str("(NO_ANNOT,FOO);")
    leaf = tree.find_node_with_taxon_label("FOO")
    with pytest.raises(ValueError, match="cannot extract species"):
        label, species = label_with_species_end(leaf)



# -------------------------
# NHX annotation parser tests
# -------------------------

@pytest.fixture
def dendro_node():
    from dendropy import Node
    return Node()

def test_nhx_to_event_D_duplication(dendro_node):
    dendro_node.annotations.add_new("D", "T")
    assert nhx_to_event(dendro_node) == "duplication"

def test_nhx_to_event_D_speciation(dendro_node):
    dendro_node.annotations.add_new("D", "0")
    assert nhx_to_event(dendro_node) == "speciation"

def test_nhx_to_event_Ev_duplication(dendro_node):
    dendro_node.annotations.add_new("Ev", "1>0>0>dup>t")
    assert nhx_to_event(dendro_node) == "duplication"

def test_nhx_to_event_Ev_speciation(dendro_node):
    dendro_node.annotations.add_new("Ev", "0>1>0>spec>t")
    assert nhx_to_event(dendro_node) == "speciation"


# -------------------------
# GeneRefHelper tests
# -------------------------

def test_gene_ref_and_species_node():
    root = ET.Element("orthoXML")
    helper = GeneRefHelper(root)

    ref_id = helper.gene("Gene1", "Human")
    assert ref_id == "1"

    # Should not create duplicate
    ref_id2 = helper.gene("Gene1", "Human")
    assert ref_id == ref_id2

    sp_nodes = root.findall("species")
    assert len(sp_nodes) == 1
    assert sp_nodes[0].attrib['name'] == "Human"

# -------------------------
# OrthoXMLBuilder integration
# -------------------------

def make_simple_tree():
    return dendropy.Tree.get(data="(geneA_HUMAN,geneB_MOUSE);", schema="newick", preserve_underscores=True)


def test_orthoxml_builder_writes_valid_xml():
    tree = make_simple_tree()
    builder = OrthoXMLBuilder(origin="unit_test")

    def dummy_label_to_event(node):
        return "speciation"

    builder.add_group(
        tree,
        label_to_event=dummy_label_to_event,
        label_to_id_and_species=label_with_species_end
    )

    output = BytesIO()
    builder.write(output)
    output.seek(0)

    doc = ET.parse(output)
    root = doc.getroot()
    assert root.tag.endswith("orthoXML")
    assert any(child.tag.endswith("groups") for child in root)
    print(output.getvalue())


# -------------------------
# End to end test
# -------------------------

@pytest.fixture
def single_nhx_example_file():
    return [Path(__file__).parent / "test-data" / "labeled_gene_trees.nwk"]


def test_convert_nhx_to_orthoxml(single_nhx_example_file):
    """ensure reading and writing of orthoxml works results in the same content"""
    out_stream = BytesIO()
    orthoxml_from_newicktrees(single_nhx_example_file, out_stream, label_to_id_and_species=nhx_species_encoded_leaf)
    out_stream.seek(0)

    print(out_stream.getvalue())
    oxml = ET.parse(out_stream)
    root = oxml.getroot()
    assert root.tag.endswith("orthoXML")
    NS = {'oxml': 'http://orthoXML.org/2011/'}
    assert len(oxml.findall("oxml:species", NS)) == 7, "Expected 7 species in the output"
    assert len(oxml.findall(".//oxml:gene", NS)) == 14, "Expected 14 genes in the output"
    assert len(oxml.findall(".//oxml:groups/oxml:orthologGroup", NS)) == 1, "Expected one ortholog group in the output"
    assert len(oxml.findall(".//oxml:geneRef", NS)) == 14, "Expected 14 gene references in the output"
    assert len(oxml.findall(".//oxml:paralogGroup", NS)) == 6, "Expected 6 paralog groups in the output"
