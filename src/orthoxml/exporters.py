# exporters.py

from .models import OrthologGroup, ParalogGroup, UnionFind, Taxon, Species


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


def get_ogs(pairs: list[(str, str)]) -> dict[str, list[str]]:
    """
    Given a list of valid gene pairs, return a dictionary mapping of representative gene to the orthologous group genes.
    Uses Union-Find to group genes.

    :param pairs: A list of orthologous gene pairs.

    :return: A dictionary of representative gene to orthologous group genes.
    """
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
