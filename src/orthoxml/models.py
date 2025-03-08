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
    __slots__ = ["_id", "geneId", "protId"]
    def __init__(self, _id: str, geneId: str, protId: str):
        self._id = _id
        self.geneId = geneId
        self.protId = protId

    def __repr__(self):
        return f"Gene(id={self._id}, geneId={self.geneId}, protId={self.protId})"
    
    @classmethod
    def from_xml(cls, xml_element):
        # xml_element is a <gene> element.
        return cls(
            _id=xml_element.get("id"),
            geneId=xml_element.get("geneId"),
            protId=xml_element.get("protId")
        )

    def to_xml(self):
        gene_el = etree.Element(f"{{{ORTHO_NS}}}gene")
        gene_el.set("id", self._id)
        if self.geneId:
            gene_el.set("geneId", self.geneId)
        if self.protId:
            gene_el.set("protId", self.protId)
        return gene_el

class Taxon:
    __slots__ = ["id", "name", "children"]
    def __init__(self, id, name, children=None):
        self.id = id
        self.name = name
        self.children = children or []  # list of Taxon objects

    def __repr__(self):
        return f"Taxon(id={self.id}, name={self.name}, children={self.children})"

    def __len__(self):
        # TODO
        return 0
    
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

    def to_str(self):
        """
        Returns a string representation of the taxonomy tree in a hierarchical format.
        Example output:
        
        LUCA
        ├── Archaea
        │   ├── KORCO
        │   ├── Euryarchaeota
        │   │   ├── HALSA
        │   │   └── THEKO
        │   └── NITMS
        ├── Bacteria
        │   └── ... (and so on)
        """
        def _child_str(node, prefix, is_last):
            # Determine the branch marker.
            branch = "└── " if is_last else "├── "
            # Build the line for this node.
            line = prefix + branch + node.name
            lines = [line]
            # Update the prefix for the children.
            new_prefix = prefix + ("    " if is_last else "│   ")
            for i, child in enumerate(node.children):
                # Check if this child is the last one.
                is_child_last = (i == len(node.children) - 1)
                lines.append(_child_str(child, new_prefix, is_child_last))
            return "\n".join(lines)

        # Start with the root node (printed without any branch symbol).
        lines = [self.name]
        # Process each child of the root.
        for i, child in enumerate(self.children):
            lines.append(_child_str(child, "", i == len(self.children) - 1))
        return "\n".join(lines)

class ParalogGroup:
    __slots__ = ["taxonId", "geneRefs", "orthologGroups", "paralogGroups"]
    def __init__(self, taxonId=None, geneRefs=None, orthologGroups=None, paralogGroups=None):
        self.taxonId = taxonId  # optional attribute (as string)
        self.geneRefs = geneRefs or []        # list of gene id strings
        self.orthologGroups = orthologGroups or []      # list of OrthologGroup objects
        self.paralogGroups = paralogGroups or []  # list of ParalogGroup objects  

    def __repr__(self):
        return f"ParalogGroup(taxonId={self.taxonId}, geneRefs={self.geneRefs}, orthologGroups={self.orthologGroups}, paralogGroups={self.paralogGroups})"
    
    def __len__(self):
        # TODO
        return 0
    
    @classmethod
    def from_xml(cls, xml_element):
        # xml_element is a <paralogGroup> element.
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
        group_el = etree.Element(f"{{{ORTHO_NS}}}paralogGroup")
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


class OrthologGroup:
    __slots__ = ["taxonId", "geneRefs", "orthologGroups", "paralogGroups"]
    def __init__(self, taxonId=None, geneRefs=None, orthologGroups=None, paralogGroups=None):
        self.taxonId = taxonId  # optional attribute (as string)
        self.geneRefs = geneRefs or []        # list of gene id strings
        self.orthologGroups = orthologGroups or []      # list of OrthologGroup objects
        self.paralogGroups = paralogGroups or []  # list of ParalogGroup objects

    def __repr__(self):
        return f"OrthologGroup(taxonId={self.taxonId}, geneRefs={self.geneRefs}, orthologGroups={self.orthologGroups}, paralogGroups={self.paralogGroups})"

    def __len__(self):
        # TODO
        return 0
    
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
    