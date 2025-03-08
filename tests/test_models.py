from orthoxml.models import Gene, ORTHO_NS
from lxml import etree


def test_gene_from_xml():
    xml = """
    <gene id="gene1" geneId="abc" protId="123"/>
    """
    gene = Gene.from_xml(etree.fromstring(xml))
    assert gene._id == "gene1"
    assert gene.geneId == "abc"
    assert gene.protId == "123"
    
    xml = """
    <gene id="gene2" geneId="def"/>
    """
    gene = Gene.from_xml(etree.fromstring(xml))
    assert gene._id == "gene2"
    assert gene.geneId == "def"
    assert gene.protId is None

    xml = """
    <gene id="gene3"/>
    """
    gene = Gene.from_xml(etree.fromstring(xml))
    assert gene._id == "gene3"
    assert gene.geneId is None
    assert gene.protId is None

    xml = """
    <gene id="gene4" protId="456"/>
    """
    gene = Gene.from_xml(etree.fromstring(xml))
    assert gene._id == "gene4"
    assert gene.geneId is None
    assert gene.protId == "456"

def test_gene_to_xml():
    pass
