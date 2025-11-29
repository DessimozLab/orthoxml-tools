from xml.etree import ElementTree as ET
from typing import Optional
import dendropy
import re
import logging
import collections
from .. import __version__ as orthoxml_tools_version
logger = logging.getLogger(__name__)


def nhx_to_event(node):
    """extract event type from the node annotatins."""
    is_D = node.annotations.get_value('D')
    if is_D is not None:
        return "duplication" if is_D.lower() in "dty" else "speciation"
    is_Ev = node.annotations.get_value('Ev')
    if is_Ev is not None:
        # Ev=duplications>speciations>gene losses>eventtype>duplication type
        try:
            nr_dup, nr_spe, *_ = is_Ev.split(">")
            if int(nr_dup) > 0:
                return "duplication"
            if int(nr_spe) > 0:
                return "speciation"
        except ValueError:
            if is_Ev.lower().startswith("d"):
                return "duplication"
            if is_Ev.lower().startswith("s"):
                return "speciation"
            logger.error(f"Cannot parse annotation Ev={is_Ev}. Invalid format")
            raise ValueError(f"Cannot parse annotation Ev={is_Ev}")
    logger.debug("no annotations found for node %s. Assuming speciation", node)
    return "speciation"


def label_with_species_end(node):
    """extract gene name, species tuple from a label of a dendropy node that encodes
    the species at the end, separated with a '_' character, i.e.
    uniprot ids.

    Example:
        >>> nd = dendropy.Tree.get(data="(P53_HUMAN, TP53_MOUSE);", schema="newick",preserve_underscores=True).leaf_nodes()[0]
        >>> label_with_species_end(nd)
        ('P53_HUMAN', 'HUMAN')
    """
    label = node.taxon.label
    k = label.rfind('_')
    if k == -1:
        # no species information found
        raise ValueError("cannot extract species from label: {}".format(label))
    return label, label[k+1:]


def nhx_species_encoded_leaf(node):
    """extract the gene name, species tuple from a nhx encoded label.

    Example::
        >>> nd = dendropy.Tree.get(data="(P53[&&NHX:S=Homo sapiens],TP53[&&NHX:S=Mus musculus]);", schema="newick",preserve_underscores=True).leaf_nodes()[0]
        >>> nhx_species_encoded_leaf(nd)
        ('P53', 'Homo sapiens')
    """
    if not node.has_annotations:
        raise ValueError(f"no NHX annotations found in {node}")
    species = node.annotations.get_value('S') or node.annotations.get_value('T')
    return node.taxon.label, species


def nhx_taxonomy_node_name(node):
    """Extract taxonomy node name from NHX annotations.

    Looks for T, S, or name annotations to identify taxonomy nodes.
    Order of precedence: T -> S -> name

    Returns None if no taxonomy annotation is found.
    """
    if not node.has_annotations:
        return None

    # Order of precedence: T -> S -> name
    taxon_name = (node.annotations.get_value('T') or
                  node.annotations.get_value('S') or
                  node.annotations.get_value('name'))
    return taxon_name


class GeneRefHelper:
    def __init__(self, xml_root):
        self._xml_root = xml_root
        self._gene2ref = {}
        self._species2node = {}
        self._genref_cnt = 0
        self._species_cnt = 0

    def _new_species_node(self, species):
        self._species_cnt += 1
        sp_node = ET.Element("species", dict(name=species, NCBITaxId=str(self._species_cnt)))
        sp_db_node = ET.SubElement(sp_node, "database", dict(name=species, version="n/a"))
        sp_db_genes_node = ET.SubElement(sp_db_node, "genes")
        self._xml_root.insert(-1, sp_node)
        self._species2node[species] = sp_db_genes_node

    def get_species_node(self, species):
        if not species in self._species2node:
            self._new_species_node(species)
        return self._species2node[species]

    def next_generef(self, species):
        self._genref_cnt += 1
        return str(self._genref_cnt)

    def gene(self, id, species):
        if not id in self._gene2ref:
            geneRef = self.next_generef(species)
            sp_node = self.get_species_node(species)
            ET.SubElement(sp_node, "gene", dict(id=str(geneRef), protId=id))
            self._gene2ref[id] = geneRef
        return self._gene2ref[id]


class SpeciesConsecutiveGeneRefsHelper(GeneRefHelper):
    def __init__(self, xml_root):
        super().__init__(xml_root)
        self._generef_per_species = {}

    def _new_species_node(self, species):
        super()._new_species_node(species)
        self._generef_per_species[species] = self._species_cnt * 100000

    def next_generef(self, species):
        if species not in self._generef_per_species:
            self._new_species_node(species)
        self._generef_per_species[species] += 1
        return str(self._generef_per_species[species])


class TaxonomyBuilder:
    """Helper class to build taxonomy trees from NHX annotations."""

    def __init__(self):
        self._taxon_id_counter = 0
        self._taxon_name_to_id = {}
        self._species_to_taxon_id = {}
        self._internal_nodes = {}
        self._taxonomy_tree = {}

    def _next_taxon_id(self):
        self._taxon_id_counter += 1
        return self._taxon_id_counter

    def _get_or_create_taxon_id(self, taxon_name):
        """Get existing or create new taxon ID for given name."""
        if taxon_name not in self._taxon_name_to_id:
            self._taxon_name_to_id[taxon_name] = self._next_taxon_id()
        return self._taxon_name_to_id[taxon_name]

    def add_species(self, species_name):
        """Add a species to the taxonomy."""
        if species_name not in self._species_to_taxon_id:
            taxon_id = self._get_or_create_taxon_id(species_name)
            self._species_to_taxon_id[species_name] = taxon_id
        return self._species_to_taxon_id[species_name]

    def build_from_tree(self, tree, label_to_id_and_species):
        """Build taxonomy structure from a dendropy tree."""

        # First, collect all species from leaves
        for leaf in tree.leaf_node_iter():
            _, species = label_to_id_and_species(leaf)
            self.add_species(species)

        # Track species assignments to avoid circular references
        species_assigned = set()

        # Build internal node hierarchy
        def process_node(node, parent_taxon_name=None):
            if node.is_leaf():
                return

            # Get taxonomy name for this internal node
            taxon_name = nhx_taxonomy_node_name(node)

            if taxon_name:
                taxon_id = self._get_or_create_taxon_id(taxon_name)
                self._internal_nodes[taxon_name] = taxon_id

                # Build parent-child relationships
                if parent_taxon_name and parent_taxon_name != taxon_name:
                    if parent_taxon_name not in self._taxonomy_tree:
                        self._taxonomy_tree[parent_taxon_name] = []
                    if taxon_name not in self._taxonomy_tree[parent_taxon_name]:
                        self._taxonomy_tree[parent_taxon_name].append(taxon_name)

                # Process children with this node as parent
                for child in node.child_node_iter():
                    if child.is_leaf():
                        # Leaf nodes belong to this internal taxon
                        _, species = label_to_id_and_species(child)
                        species_key = f"{taxon_name}:{species}"
                        if species_key not in species_assigned:
                            if taxon_name not in self._taxonomy_tree:
                                self._taxonomy_tree[taxon_name] = []
                            if species not in self._taxonomy_tree[taxon_name]:
                                self._taxonomy_tree[taxon_name].append(species)
                            species_assigned.add(species_key)
                    else:
                        process_node(child, taxon_name)
            else:
                # No taxonomy name, process children with same parent
                for child in node.child_node_iter():
                    process_node(child, parent_taxon_name)

        # Start from root
        if tree.seed_node:
            process_node(tree.seed_node)

    def build_xml(self, parent_element):
        """Build XML taxonomy section."""
        if not self._taxonomy_tree and not self._internal_nodes:
            return None

        taxonomy_elem = ET.SubElement(parent_element, "taxonomy")

        # Find root nodes (nodes that are not children of any other node)
        all_children = set()
        for parent, children in self._taxonomy_tree.items():
            all_children.update(children)

        root_nodes = []
        for node in self._taxonomy_tree.keys():
            if node not in all_children:
                root_nodes.append(node)

        # If no clear root structure, create an artificial root
        if not root_nodes and self._taxonomy_tree:
            # Create a root node that contains all top-level taxa
            root_taxon_id = self._get_or_create_taxon_id("Root")
            root_elem = ET.SubElement(taxonomy_elem, "taxon",
                                    {"id": str(root_taxon_id), "name": "Root"})

            # Add all taxa as children of root
            for taxon_name in self._taxonomy_tree.keys():
                self._build_taxon_element(root_elem, taxon_name)
        else:
            # Build tree from identified roots
            for root_node in root_nodes:
                self._build_taxon_element(taxonomy_elem, root_node)

        return taxonomy_elem

    def _build_taxon_element(self, parent_elem, taxon_name, visited=None):
        """Recursively build taxon elements."""
        if visited is None:
            visited = set()

        # Prevent infinite recursion from circular references
        if taxon_name in visited:
            return
        visited.add(taxon_name)

        taxon_id = self._taxon_name_to_id.get(taxon_name)
        if not taxon_id:
            # This is a species
            taxon_id = self._species_to_taxon_id.get(taxon_name)

        if taxon_id:
            taxon_elem = ET.SubElement(parent_elem, "taxon",
                                     {"id": str(taxon_id), "name": taxon_name})

            # Add children if any
            children = self._taxonomy_tree.get(taxon_name, [])
            for child_name in children:
                self._build_taxon_element(taxon_elem, child_name, visited.copy())


class OrthoXMLBuilder:
    generef_helper = GeneRefHelper

    def __init__(self, origin="newick_tree_import", version: Optional[str] = None):
        version = version or f"orthoxml-tools:{orthoxml_tools_version}"
        self.NS = "http://orthoXML.org/2011/"
        ET.register_namespace("", self.NS)
        self.doc = ET.ElementTree()
        self.xml = ET.Element('orthoXML', dict(xmlns=self.NS, origin=origin, version="0.5", originVersion=version))
        self.grps = ET.SubElement(self.xml, "groups")
        self.doc._setroot(self.xml)
        self._genehelper = self.generef_helper(self.xml)
        self._grpCnt = 0
        self._taxonomy_builder = TaxonomyBuilder()

    def write(self, fh):
        # Add taxonomy section after species and before groups
        # Find the index of the groups element
        groups_index = list(self.xml).index(self.grps)

        # Build taxonomy and insert it before groups
        taxonomy_element = self._taxonomy_builder.build_xml(self.xml)
        if taxonomy_element is not None:
            # Move taxonomy to correct position (before groups)
            self.xml.remove(taxonomy_element)
            self.xml.insert(groups_index, taxonomy_element)

        self.add_loft_ids()
        ET.indent(self.xml)
        self.doc.write(fh, xml_declaration=True, encoding="UTF-8", default_namespace=None)

    def add_loft_ids(self, id_formatter="HOG:{:08d}"):

        def encodeParalogClusterId(prefix, nr):
            letters = []
            while nr // 26 > 0:
                letters.append(chr(97 + (nr % 26)))
                nr = nr // 26 - 1
            letters.append(chr(97 + (nr % 26)))
            return prefix + ''.join(letters[::-1])

        def nextSubHogId(idx):
            dups[idx] += 1
            return dups[idx]

        def rec_annotate(node, og, idx=0):
            if node.tag == "orthologGroup":
                #taxon_id = node.get('taxonId')
                node.set('id', og ) #+ "_" + taxon_id)
                for child in list(node):
                    rec_annotate(child, og, idx)
            elif node.tag == "paralogGroup":
                idx += 1
                next_og = "{}.{}".format(og, nextSubHogId(idx))
                for i, child in enumerate(list(node)):
                    rec_annotate(child, encodeParalogClusterId(next_og, i), idx)

        for i, el in enumerate(self.xml.findall(".//groups/orthologGroup")):
            og = id_formatter.format(i+1)
            dups = collections.defaultdict(int)
            rec_annotate(el, og)

    def add_group(self, tree, label_to_event, label_to_id_and_species):
        self._grpCnt += 1
        hogid = "HOG:{:07d}".format(self._grpCnt)

        # Build taxonomy information from this tree
        self._taxonomy_builder.build_from_tree(tree, label_to_id_and_species)

        for idx, node in enumerate(tree.preorder_node_iter(), start=1):

            if node.is_leaf():
                id, species = label_to_id_and_species(node)
                parent_hog = node.parent_node.hog_node
                geneRef = ET.SubElement(
                    parent_hog,
                    "geneRef",
                    dict(id=self._genehelper.gene(id, species)) #, loft=parent_hog['loft'])
                )
            else:
                parent_grp_node = None
                typ = label_to_event(node)
                if typ == "duplication" and node.parent_node is None:
                    # introduce a artificial topOG
                    hog = ET.Element('orthologGroup', dict(id="{}_0".format(hogid), loft=hogid))
                    parent_grp_node = hog
                    self.grps.append(hog)
                elif node.parent_node is not None:
                    parent_grp_node = node.parent_node.hog_node

                gtyp = "orthologGroup" if typ=="speciation" else "paralogGroup"
                if parent_grp_node is None:
                    hog_node = ET.SubElement(self.grps, gtyp, dict(id="{}".format(hogid)))
                else:
                    hog_node = ET.SubElement(parent_grp_node, gtyp)
                node.hog_node = hog_node


class FlatOrthoXMLBuilder(OrthoXMLBuilder):
    generef_helper = SpeciesConsecutiveGeneRefsHelper
    def add_group(self, tree, label_to_event, label_to_id_and_species):
        self._grpCnt += 1
        hogid = "HOG:{:07d}".format(self._grpCnt)
        hog_node = ET.SubElement(self.grps, "orthologGroup", dict(id="{}".format(hogid)))
        ET.SubElement(hog_node, "property", dict(name="TaxRange", value="Root"))
        for node in tree.leaf_node_iter():
            id, species = label_to_id_and_species(node.taxon.label)
            geneRef = ET.SubElement(
                hog_node,
                "geneRef",
                dict(id=self._genehelper.gene(id, species))  # , loft=parent_hog['loft'])
            )


class SpecialLevelOrthoXMLBuilder(OrthoXMLBuilder):
    def __init__(self, species_tree: dendropy.Tree, origin="newick_tree_import", version="n/a"):
        super().__init__(origin, version)


def orthoxml_from_newicktrees(nwkfiles, output, label_to_event=None, label_to_id_and_species=None):
    """Converts a list of labeled gene trees in newick format into a single orthoxml file.

    The gene trees are assumed to be correctly rooted. Speciation event nodes will be converted into
    orthologGroup elements, and duplication events into paralogGroup elements. In addition, the HOGs
    will be annotated with LOFT ids.

    The parameter `label_to_event` takes accepts a function which should return either the string
    "speciation" or "duplication" from the label of an internal node (only called with internal
    tree nodes). If not set, the default function will extract the event from a nhx annotation
    using the function :py:func:`nhx_to_event`.

    The parameter `label_to_id_and_species` accepts a function which should return a tuple with the
    gene id and it's species from the leaf label of a tree node (only called with leaf nodes).
    By default, the function :py:func:`label_with_species_end` is used. The function
    :py:func:`nhx_species_encoded_leaf` could be used if the species information is encoded in nhx
    format, or a user defined function is possible.

    :param nwkfiles: A list of paths containing gene trees in newick format
    :param output:   An output file path where the resulting orthoxml file is stored.
    :param label_to_event:  see above
    :param label_to_id_and_species:  see above
    """
    if label_to_event is None:
        label_to_event = nhx_to_event
    if label_to_id_and_species is None:
        label_to_id_and_species = label_with_species_end
    builder = OrthoXMLBuilder()
    for nwkfile in nwkfiles:
        for c, tree in enumerate(dendropy.TreeList.get_from_path(nwkfile, schema="newick", preserve_underscores=True), 1):
            builder.add_group(tree, label_to_event=label_to_event, label_to_id_and_species=label_to_id_and_species)
        logger.info("processed %d trees from %s", c, nwkfile)
    builder.write(output)