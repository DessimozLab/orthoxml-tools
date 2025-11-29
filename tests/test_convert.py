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
    nhx_taxonomy_node_name,
    TaxonomyBuilder,
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
# Taxonomy function tests
# -------------------------

def test_nhx_taxonomy_node_name():
    """Test extracting taxonomy node names from NHX annotations"""
    # Test T annotation
    tree = make_tree_from_str("(A,B)[&&NHX:T=NODE_1];")
    root = tree.seed_node
    assert nhx_taxonomy_node_name(root) == "NODE_1"

    # Test name annotation
    tree = make_tree_from_str("(A,B)[&&NHX:name=NODE_2];")
    root = tree.seed_node
    assert nhx_taxonomy_node_name(root) == "NODE_2"

    # Test S annotation (for leaf nodes)
    tree = make_tree_from_str("(A,B)[&&NHX:S=SPECIES_1];")
    root = tree.seed_node
    assert nhx_taxonomy_node_name(root) == "SPECIES_1"

    # Test T takes precedence over name
    tree = make_tree_from_str("(A,B)[&&NHX:T=NODE_T:name=NODE_NAME];")
    root = tree.seed_node
    assert nhx_taxonomy_node_name(root) == "NODE_T"

    # Test no annotations
    tree = make_tree_from_str("(A,B);")
    root = tree.seed_node
    assert nhx_taxonomy_node_name(root) is None


def test_taxonomy_builder():
    """Test TaxonomyBuilder functionality"""
    builder = TaxonomyBuilder()

    # Add some species
    builder.add_species("Human")
    builder.add_species("Mouse")

    # Check species IDs are assigned
    assert builder._species_to_taxon_id["Human"] == 1
    assert builder._species_to_taxon_id["Mouse"] == 2


def test_taxonomy_from_nhx_tree():
    """Test building taxonomy from NHX tree with T annotations"""
    newick_str = "((A[&&NHX:S=HUMAN]:0.1,B[&&NHX:S=MOUSE]:0.1)[&&NHX:T=MAMMALS]:0.1,C[&&NHX:S=FISH]:0.2)[&&NHX:T=VERTEBRATES];"
    tree = make_tree_from_str(newick_str)

    builder = TaxonomyBuilder()
    builder.build_from_tree(tree, nhx_species_encoded_leaf)

    # Check species were added
    assert "HUMAN" in builder._species_to_taxon_id
    assert "MOUSE" in builder._species_to_taxon_id
    assert "FISH" in builder._species_to_taxon_id

    # Check internal nodes were added
    assert "MAMMALS" in builder._internal_nodes
    assert "VERTEBRATES" in builder._internal_nodes

    # Check tree structure
    assert "VERTEBRATES" in builder._taxonomy_tree
    assert "MAMMALS" in builder._taxonomy_tree["VERTEBRATES"]
    assert "FISH" in builder._taxonomy_tree["VERTEBRATES"]
    assert "HUMAN" in builder._taxonomy_tree["MAMMALS"]
    assert "MOUSE" in builder._taxonomy_tree["MAMMALS"]


def test_taxonomy_from_nhx_tree_name():
    """Test building taxonomy from NHX tree with name annotations"""
    newick_str = "((A[&&NHX:S=HUMAN]:0.1,B[&&NHX:S=MOUSE]:0.1)[&&NHX:name=MAMMALS]:0.1,C[&&NHX:S=FISH]:0.2)[&&NHX:name=VERTEBRATES];"
    tree = make_tree_from_str(newick_str)

    builder = TaxonomyBuilder()
    builder.build_from_tree(tree, nhx_species_encoded_leaf)

    # Check internal nodes were added
    assert "MAMMALS" in builder._internal_nodes
    assert "VERTEBRATES" in builder._internal_nodes


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


@pytest.fixture
def multiple_nhx_example_file():
    return [Path(__file__).parent / "test-data" / "multiple_gene_trees.nwk"]


def test_convert_multi_nhx_to_orthoxml(multiple_nhx_example_file):
    """ensure reading and writing of orthoxml works results in the same content"""
    out_stream = BytesIO()
    orthoxml_from_newicktrees(multiple_nhx_example_file, out_stream, label_to_id_and_species=nhx_species_encoded_leaf)
    out_stream.seek(0)

    print(out_stream.getvalue())
    oxml = ET.parse(out_stream)
    root = oxml.getroot()
    assert root.tag.endswith("orthoXML")
    NS = {'oxml': 'http://orthoXML.org/2011/'}
    assert len(oxml.findall("oxml:species", NS)) == 6, "Expected 6 species in the output"
    assert len(oxml.findall(".//oxml:gene", NS)) == 14, "Expected 14 genes in the output"
    assert len(oxml.findall(".//oxml:groups/oxml:orthologGroup", NS)) == 2, "Expected one ortholog group in the output"
    assert len(oxml.findall(".//oxml:geneRef", NS)) == 14, "Expected 14 gene references in the output"
    assert len(oxml.findall(".//oxml:paralogGroup", NS)) == 4, "Expected 4 paralog groups in the output"


def test_taxonomy_xml_generation():
    """Test that OrthoXMLBuilder includes taxonomy section"""
    newick_str = "((A[&&NHX:S=HUMAN]:0.1,B[&&NHX:S=MOUSE]:0.1)[&&NHX:T=MAMMALS]:0.1,C[&&NHX:S=FISH]:0.2)[&&NHX:T=VERTEBRATES];"
    tree = make_tree_from_str(newick_str)

    builder = OrthoXMLBuilder(origin="test_taxonomy")
    builder.add_group(tree, label_to_event=nhx_to_event, label_to_id_and_species=nhx_species_encoded_leaf)

    output = BytesIO()
    builder.write(output)
    output.seek(0)

    doc = ET.parse(output)
    root = doc.getroot()

    # Check taxonomy section exists
    NS = {'oxml': 'http://orthoXML.org/2011/'}
    taxonomy_nodes = root.findall("oxml:taxonomy", NS)
    assert len(taxonomy_nodes) == 1, "Expected exactly one taxonomy section"

    # Check taxon elements exist
    taxon_nodes = root.findall(".//oxml:taxon", NS)
    assert len(taxon_nodes) >= 5, f"Expected at least 5 taxon nodes, got {len(taxon_nodes)}"

    # Check species are included
    taxon_names = [node.get('name') for node in taxon_nodes]
    assert 'HUMAN' in taxon_names
    assert 'MOUSE' in taxon_names
    assert 'FISH' in taxon_names
    assert 'MAMMALS' in taxon_names
    assert 'VERTEBRATES' in taxon_names
