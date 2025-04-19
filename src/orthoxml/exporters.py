# exporters.py

from .models import OrthologGroup, ParalogGroup, UnionFind, Taxon, Species
from .logger import logger

def get_ortho_pairs_iterative(root):
    """
    Iteratively traverse the tree starting at the given root node and yield valid ortholog pairs.
    A valid pair is yielded when the lowest common ancestor of the genes is not a ParalogGroup.
    
    This function uses an explicit stack to perform a postorder traversal.
    Each stack frame is a dictionary with:
      - 'node': the current node,
      - 'children': a list of this node's children (both ortholog and paralog groups),
      - 'child_index': index for the next child to process,
      - 'child_gene_refs_list': a list where each element is a list of gene references from one child subtree.
    
    The function yields a tuple (geneRef1, geneRef2) for every valid pair discovered.
    """
    stack = []
    # Push the root node frame onto the stack.
    stack.append({
        'node': root,
        'children': root.orthologGroups + root.paralogGroups,
        'child_index': 0,
        'child_gene_refs_list': []
    })
    
    while stack:
        current = stack[-1]
        node = current['node']
        children = current['children']
        child_index = current['child_index']
        
        # If not all children have been processed, push the next child onto the stack.
        if child_index < len(children):
            child = children[child_index]
            current['child_index'] += 1  # increment the child index for the current node
            stack.append({
                'node': child,
                'children': child.orthologGroups + child.paralogGroups,
                'child_index': 0,
                'child_gene_refs_list': []
            })
        else:
            # All children for this node have been processed.
            # Combine the geneRefs of the current node with the geneRefs from all children.
            cur_gene_refs = list(node.geneRefs)  # starting with the current node's gene references
            for child_refs in current['child_gene_refs_list']:
                cur_gene_refs.extend(child_refs)
            
            # Only yield pairs if the current node is NOT a ParalogGroup.
            if not isinstance(node, ParalogGroup):
                # 1. Pair each geneRef in the current node with each geneRef from each child branch.
                for branch_refs in current['child_gene_refs_list']:
                    for r in node.geneRefs:
                        for s in branch_refs:
                            yield (r, s)
                            
                # 2. Pair geneRefs from different child branches.
                branch_lists = current['child_gene_refs_list']
                for i in range(len(branch_lists)):
                    for j in range(i + 1, len(branch_lists)):
                        for r in branch_lists[i]:
                            for s in branch_lists[j]:
                                yield (r, s)
                                
                # 3. Pair geneRefs within the current node itself.
                own_refs = list(node.geneRefs)
                for i in range(len(own_refs)):
                    for j in range(i + 1, len(own_refs)):
                        yield (own_refs[i], own_refs[j])
            
            # Finished processing the current node; pop the frame off the stack.
            stack.pop()
            
            # If there is a parent frame, add the aggregated gene references for the current node
            # to the parent's 'child_gene_refs_list'. This passes the information upward.
            if stack:
                stack[-1]['child_gene_refs_list'].append(cur_gene_refs)


def get_ortho_pairs_recursive(node: OrthologGroup) -> tuple[list[str], list[(str, str)]]:
    """
    Recursively traverse the tree and return a tuple:
      (all_gene_refs_in_subtree, valid_pairs)
    where valid_pairs is a list of tuple pairs (r, s) of geneRefs for which
    the lowest common ancestor is not a ParalogGroup.

    :param node: The OrthologGroup node to start the traversal from.

    :return: A tuple of all orthologous pairs.
    """
    # Start with geneRefs at the current node.
    gene_refs = list(node.geneRefs)
    pairs = []
    
    # Process both ortholog and paralog children.
    # We also want to keep track of the geneRefs from each child separately
    # so that we only pair refs coming from different branches at a non-paralog node.
    child_gene_refs_list = []
    for child in node.orthologGroups + node.paralogGroups:
        child_refs, child_pairs = get_ortho_pairs_recursive(child)
        pairs.extend(child_pairs)
        child_gene_refs_list.append(child_refs)
        gene_refs.extend(child_refs)

    # If the current node is not a ParalogGroup, then the geneRefs coming from
    # different child branches (or from the current node vs. a child branch)
    # have their closest common ancestor at this node.
    # We only form these pairs if this node is non-paralog.
    if not isinstance(node, ParalogGroup):
        # Pair the current node's own geneRefs with each child's refs.
        for child_refs in child_gene_refs_list:
            for r in node.geneRefs:
                for s in child_refs:
                    pairs.append((r, s))

        # Also pair geneRefs coming from different children branches.
        for i in range(len(child_gene_refs_list)):
            for j in range(i+1, len(child_gene_refs_list)):
                for r in child_gene_refs_list[i]:
                    for s in child_gene_refs_list[j]:
                        pairs.append((r, s))

        # Also pair geneRefs coming from the current node with each other.
        for i in range(len(node.geneRefs)):
            for j in range(i+1, len(node.geneRefs)):
                pairs.append((list(node.geneRefs)[i], list(node.geneRefs)[j]))

    # If the current node is a ParalogGroup, then its children are “merged”
    # but we do not count pairs at this level because then the LCA would be a ParalogGroup.
    
    return gene_refs, pairs


def get_paralog_pairs_recursive(node: OrthologGroup) -> tuple[list[str], list[(str, str)]]:
    """
    Recursively traverse the tree and return a tuple:
      (all_gene_refs_in_subtree, valid_pairs)
    where valid_pairs is a list of tuple pairs (r, s) of geneRefs for which
    the lowest common ancestor is a ParalogGroup.

    :param node: The OrthologGroup node to start the traversal from.
    :return: A tuple of all paralogous pairs.
    """
    # Start with geneRefs at the current node.
    gene_refs = list(node.geneRefs)
    pairs = []

    # Process both ortholog and paralog children.
    # We also keep track of the geneRefs from each child separately so that we only
    # pair refs coming from different branches when the current node is a ParalogGroup.
    child_gene_refs_list = []
    for child in node.orthologGroups + node.paralogGroups:
        child_refs, child_pairs = get_paralog_pairs_recursive(child)
        pairs.extend(child_pairs)
        child_gene_refs_list.append(child_refs)
        gene_refs.extend(child_refs)

    # If the current node is a ParalogGroup, then geneRefs coming from different
    # branches (or from the current node vs. a child branch) have their lowest common
    # ancestor at this node, and are considered paralogous.
    if isinstance(node, ParalogGroup):
        # Pair the current node's own geneRefs with each child's refs.
        for child_refs in child_gene_refs_list:
            for r in node.geneRefs:
                for s in child_refs:
                    pairs.append((r, s))

        # Pair geneRefs coming from different children branches.
        for i in range(len(child_gene_refs_list)):
            for j in range(i + 1, len(child_gene_refs_list)):
                for r in child_gene_refs_list[i]:
                    for s in child_gene_refs_list[j]:
                        pairs.append((r, s))

        # Also pair geneRefs coming from the current node with each other.
        for i in range(len(node.geneRefs)):
            for j in range(i + 1, len(node.geneRefs)):
                pairs.append((list(node.geneRefs)[i], list(node.geneRefs)[j]))

    return gene_refs, pairs


def get_maximal_og(group: OrthologGroup, species_dic: dict[str, str]) -> list[str]:
    """
    Given a tree of OrthologGroup (and ParalogGroup) nodes, process any duplication events by
    keeping only the branch (child group) that has the largest number of distinct species, and then
    return a list of gene references (leaves) that remain.
    
    The function uses a global species_dic mapping gene names (after adjusting if they start with "no_")
    to species names.
    
    Assumptions:
      - A node is considered a leaf if it has no children (i.e. both orthologGroups and paralogGroups are empty).
      - A duplication event is indicated by a node having an attribute Ev set to "duplication".
      - Gene references in a node are stored in the list geneRefs.
    """
    # TODO: remove species_dic

    excluded_genes = set()
    
    def is_leaf(node) -> bool:
        """Return True if the node has no children groups."""
        return not (node.orthologGroups or node.paralogGroups)
    
    def traverse_nodes(node):
        """Generator to yield nodes in pre-order (current node then all children)."""
        yield node
        for child in node.orthologGroups:
            yield from traverse_nodes(child)
        for child in node.paralogGroups:
            yield from traverse_nodes(child)
    
    # First pass: mark (exclude) gene leaves for non-maximal children of duplication events.
    for node in traverse_nodes(group):
        # Check if the node is internal and a duplication event.
        if not is_leaf(node) and isinstance(node, ParalogGroup):
            # For children we combine both the ortholog and paralog children.
            children = node.orthologGroups + node.paralogGroups
            list_num_species = []
            # Evaluate distinct species counts for each child branch.
            for child in children:
                species_set = set()
                # Use the convenience function get_all_leaves() to fetch gene ids from the child
                for gene in child.get_all_leaves():
                    if gene in species_dic:
                        species_set.add(species_dic[gene])
                    else:
                        logger.debug(f"species not in the dic {gene}")
                list_num_species.append(len(species_set))
            # Optionally, warn about potential polytomy (more than one child has the same maximum species count).
            if list_num_species.count(max(list_num_species)) > 1:
                logger.debug(f"Polytomy detected: more than one child has the maximal number of species. node: {node}, list: {list_num_species}")
            
            # Choose the child with the maximum number of species.
            if children:
                max_species = max(list_num_species)
                index_max = list_num_species.index(max_species)
                # For all children other than the one selected, exclude all their gene leaves.
                for i, child in enumerate(children):
                    if i != index_max:
                        for gene in child.get_all_leaves():
                            excluded_genes.add(gene)
            if node.geneRefs:
                for gene in node.geneRefs:
                    excluded_genes.add(gene)

        if is_leaf(node) and isinstance(node, ParalogGroup):
            # handle paralog groups without ortholog or paralog groups and only generefs
            # Keep the first gene and remove the rest
            for gene in node.geneRefs[1:]:
                excluded_genes.add(gene)

    
    def collect_leaves(node) -> list[str]:
        """
        Recursively collect gene references from the tree that have not been marked as excluded.
        """
        leaves = []
        # Add gene references from the current node if they have not been excluded.
        for gene in node.geneRefs:
            if gene not in excluded_genes:
                leaves.append(gene)
        # Recurse into the children groups.
        for child in node.orthologGroups:
            leaves.extend(collect_leaves(child))
        for child in node.paralogGroups:
            leaves.extend(collect_leaves(child))
        return leaves

    # The list of maximal ortholog gene references is just the collection
    # of all leaves (geneRefs) from the tree that were not excluded.
    return collect_leaves(group)


def get_ogs(pairs: list[(str, str)]) -> dict[str, list[str]]:
    """
    Given a list of valid gene pairs, return a dictionary mapping of representative gene to the orthologous group genes.
    Uses Union-Find to group genes.

    :param pairs: A list of orthologous gene pairs.

    :return: A dictionary of representative gene to orthologous group genes.
    """
    raise Warning("This algorithm behaves not as expected!")
    # Create Union-Find structure
    uf = UnionFind()

    # Process all pairs
    for a, b in pairs:
        uf.union(a, b)

    # Collect groups based on root parent
    groups = {}
    for x in uf.parent:
        root = uf.find(x)
        groups.setdefault(root, []).append(x)
    
    return groups

def compute_gene_counts_per_level(taxonomy: Taxon, species: list[Species]) -> dict[str, int]:
    """
    Compute the number of genes per taxon level.

    :param taxonomy: The taxonomy tree
    :param species: The list of species

    :return: A dictionary with the taxonId as key and the number of genes as value
    """
    # Create a dictionary to hold gene counts
    gene_counts = {}

    # Initialize counts for species using your species_list
    for sp in species:
        gene_counts[sp.taxonId] = len(sp.genes)

    # Traverse the tree using DFS and collect nodes in postorder
    stack = [taxonomy]
    postorder = []  # this will store nodes in the order they are finished

    while stack:
        node = stack.pop()
        postorder.append(node)
        for child in node.children:
            stack.append(child)

    # Process nodes in reverse postorder (bottom-up)
    while postorder:
        node = postorder.pop()  # processing from leaves upward
        own_count = gene_counts.get(node.id, 0)
        total_count = own_count
        for child in node.children:
            total_count += gene_counts.get(child.id, 0)
        gene_counts[node.id] = total_count  # store the total gene count externally

    return gene_counts


# Code from zoo/hog/convert.py for exporting OrthoXML as Gene Tree with Newick format

class TaxonNHXMixin:
    def get_tax_nhx(self):
        tags = []
        if self.level:
            tags.append(":S={}".format(self.level))
        if self.taxid:
            tags.append(":T={}".format(self.taxid))
        return tags

class Speciation:
    type = None

    def __init__(self, parent=None):
        self.level = ""
        self.taxid = None
        self.children = []
        self.parent = parent
        if parent is not None:
            parent.add_child(self)

    def add_child(self, e):
        self.children.append(e)

    def set_level(self, level):
        self.level = level

    def set_taxid(self, taxid):
        self.taxid = taxid

    def get_newick_node_name(self):
        if not hasattr(self, 'get_tax_nhx'):
            return self.level.replace(' ', '_')
        return ""

    def as_nhx(self):
        nhx = "[&&NHX"
        t = ",".join([c.as_nhx() for c in self.children])
        if t != "":
            t = "({})".format(t)
        tags = self.get_newick_node_name()

        if self.type:
            nhx += ":Ev={}".format(self.type)
        if hasattr(self, "get_tax_nhx"):
            nhx += "".join(self.get_tax_nhx())
        nhx += "]"
        if len(nhx) > 7:
            tags += nhx
        return "{}{}".format(t, tags)

class Duplication(Speciation):
    type = "duplication"

class Leaf(Speciation):
    def __init__(self, xref, species, parent=None):
        super().__init__(parent=parent)
        self.name = xref
        self.level = species

    def get_newick_node_name(self):
        return self.name

class NHXSpeciation(Speciation, TaxonNHXMixin):
    pass

class NHXDuplication(Duplication, TaxonNHXMixin):
    pass

class NHXLeaf(Leaf, TaxonNHXMixin):
    pass

class OrthoxmlToNewick:

    def __init__(self, xref_tag="protId", encode_levels_as_nhx=True, return_gene_to_species=False):
        self.xref_tag = xref_tag
        self.gene2xref = {}
        self.trees = {}
        self.depth = 0
        self.famid = None
        self.cur_event = None
        self.cur_species = None
        self._use_nhx = encode_levels_as_nhx
        self._return_gene_to_species = return_gene_to_species

    def start(self, tag, attrib):
        if tag == "{http://orthoXML.org/2011/}species":
            self.cur_species = attrib['name']
        if tag == "{http://orthoXML.org/2011/}gene":
            self.gene2xref[attrib['id']] = (attrib[self.xref_tag], self.cur_species)
        elif tag == "{http://orthoXML.org/2011/}geneRef":
            leaf_cls = NHXLeaf if self._use_nhx else Leaf
            self.cur_event.add_child(leaf_cls(*self.gene2xref[attrib['id']]))
        elif tag == "{http://orthoXML.org/2011/}orthologGroup":
            if self.depth == 0:
                self.famid = attrib.get('id', f"family_{len(self.trees.keys()) + 1}")
            speciation_cls = NHXSpeciation if self._use_nhx else Speciation
            self.cur_event = speciation_cls(self.cur_event)
            self.depth += 1
        elif tag == "{http://orthoXML.org/2011/}paralogGroup":
            dupl_cls = NHXDuplication if self._use_nhx else Duplication
            self.cur_event = dupl_cls(self.cur_event)
        elif tag == "{http://orthoXML.org/2011/}property":
            if attrib['name'] == "TaxRange":
                self.cur_event.set_level(attrib['value'])
            elif attrib['name'].lower() in ("taxid", "taxonid", "taxon_id", "ncbi_taxon_id"):
                self.cur_event.set_taxid(attrib['value'])

    def end(self, tag):
        if tag == "{http://orthoXML.org/2011/}paralogGroup":
            self.cur_event = self.cur_event.parent
        elif tag == "{http://orthoXML.org/2011/}orthologGroup":
            self.depth -= 1
            if self.depth == 0:
                assert(self.cur_event.parent is None)
                self.trees[self.famid] = self.cur_event.as_nhx() + ";"
            self.cur_event = self.cur_event.parent

    def close(self):
        if self._return_gene_to_species:
            gene2species = {k[0]: k[1] for k in self.gene2xref.values()}
            return self.trees, gene2species
        return self.trees
