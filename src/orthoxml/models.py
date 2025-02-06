# models.py
from lxml import etree

# Define the orthoXML namespace and namespace map.
ORTHO_NS = "http://orthoXML.org/2011/"
NSMAP = {None: ORTHO_NS}


class Species:
    __slots__ = ["name", "NCBITaxId", "genes"]
    def __init__(self, name, NCBITaxId, genes=None):
        self.name = name
        self.NCBITaxId = NCBITaxId
        self.genes = genes or []  # list of Gene objects
    
    def __repr__(self):
        return f"Species(name={self.name}, NCBITaxId={self.NCBITaxId}, genes={self.genes})"
    
    @classmethod
    def from_xml(cls, xml_element):
        # xml_element is a <species> element.
        name = xml_element.get("name")
        taxid = xml_element.get("NCBITaxId")
        genes = []
        # Find all gene elements (searching inside the species element).
        for gene_el in xml_element.xpath(".//ortho:gene", namespaces={"ortho": ORTHO_NS}):
            genes.append(Gene.from_xml(gene_el))
        return cls(name, taxid, genes)

    def to_xml(self):
        species_el = etree.Element(f"{{{ORTHO_NS}}}species")
        species_el.set("name", self.name)
        species_el.set("NCBITaxId", self.NCBITaxId)
        # Create a <database> element (adjust these attributes as needed).
        database_el = etree.SubElement(species_el, f"{{{ORTHO_NS}}}database")
        database_el.set("name", "someDB")
        database_el.set("version", "42")
        genes_el = etree.SubElement(database_el, f"{{{ORTHO_NS}}}genes")
        for gene in self.genes:
            genes_el.append(gene.to_xml())
        return species_el

class Gene:
    __slots__ = ["_id", "geneId"]
    def __init__(self, _id: str, geneId: str):
        self._id = _id
        self.geneId = geneId

    def __repr__(self):
        return f"Gene(id={self._id}, geneId={self.geneId})"
    
    @classmethod
    def from_xml(cls, xml_element):
        # xml_element is a <gene> element.
        return cls(
            _id=xml_element.get("id"),
            geneId=xml_element.get("geneId")
        )

    def to_xml(self):
        gene_el = etree.Element(f"{{{ORTHO_NS}}}gene")
        gene_el.set("id", self._id)
        gene_el.set("geneId", self.geneId)
        return gene_el

class Taxon:
    __slots__ = ["id", "name", "children"]
    def __init__(self, id, name, children=None):
        self.id = id
        self.name = name
        self.children = children or []  # list of Taxon objects

    def __repr__(self):
        return f"Taxon(id={self.id}, name={self.name}, children={self.children})"

    @classmethod
    def from_xml(cls, xml_element):
        # xml_element is a <taxon> element.
        taxon_id = xml_element.get("id")
        name = xml_element.get("name")
        children = []
        # Parse any nested <taxon> elements.
        for child in xml_element.xpath("./ortho:taxon", namespaces={"ortho": ORTHO_NS}):
            children.append(Taxon.from_xml(child))
        return cls(taxon_id, name, children)

    def to_xml(self):
        taxon_el = etree.Element(f"{{{ORTHO_NS}}}taxon")
        taxon_el.set("id", self.id)
        taxon_el.set("name", self.name)
        for child in self.children:
            taxon_el.append(child.to_xml())
        return taxon_el


class ParalogGroup:
    __slots__ = ["geneRefs"]
    def __init__(self, geneRefs=None):
        self.geneRefs = geneRefs or []  # list of gene id strings

    def __repr__(self):
        return f"ParalogGroup(geneRefs={self.geneRefs})"

    @classmethod
    def from_xml(cls, xml_element):
        # xml_element is a <paralogGroup> element.
        geneRefs = []
        for gene_ref in xml_element.xpath("./ortho:geneRef", namespaces={"ortho": ORTHO_NS}):
            geneRefs.append(gene_ref.get("id"))
        return cls(geneRefs)

    def to_xml(self):
        paralog_el = etree.Element(f"{{{ORTHO_NS}}}paralogGroup")
        for geneRef in self.geneRefs:
            gene_ref_el = etree.SubElement(paralog_el, f"{{{ORTHO_NS}}}geneRef")
            gene_ref_el.set("id", geneRef)
        return paralog_el


class OrthologGroup:
    __slots__ = ["taxonId", "geneRefs", "subgroups", "paralogGroups"]
    def __init__(self, taxonId=None, geneRefs=None, subgroups=None, paralogGroups=None):
        self.taxonId = taxonId  # optional attribute (as string)
        self.geneRefs = geneRefs or []        # list of gene id strings
        self.subgroups = subgroups or []      # list of OrthologGroup objects
        self.paralogGroups = paralogGroups or []  # list of ParalogGroup objects

    def __repr__(self):
        return f"OrthologGroup(taxonId={self.taxonId}, geneRefs={self.geneRefs}, subgroups={self.subgroups}, paralogGroups={self.paralogGroups})"

    @classmethod
    def from_xml(cls, xml_element):
        # xml_element is an <orthologGroup> element.
        taxonId = xml_element.get("taxonId")
        geneRefs = []
        subgroups = []
        paralogGroups = []
        # Process child elements.
        for child in xml_element:
            tag = etree.QName(child.tag).localname
            if tag == "geneRef":
                geneRefs.append(child.get("id"))
            elif tag == "orthologGroup":
                subgroups.append(OrthologGroup.from_xml(child))
            elif tag == "paralogGroup":
                paralogGroups.append(ParalogGroup.from_xml(child))
        return cls(taxonId, geneRefs, subgroups, paralogGroups)

    def to_xml(self):
        group_el = etree.Element(f"{{{ORTHO_NS}}}orthologGroup")
        if self.taxonId:
            group_el.set("taxonId", self.taxonId)
        # Note: If order matters you may want to store children in a single list.
        for subgroup in self.subgroups:
            group_el.append(subgroup.to_xml())
        for paralog in self.paralogGroups:
            group_el.append(paralog.to_xml())
        for geneRef in self.geneRefs:
            gene_ref_el = etree.SubElement(group_el, f"{{{ORTHO_NS}}}geneRef")
            gene_ref_el.set("id", geneRef)
        return group_el
    