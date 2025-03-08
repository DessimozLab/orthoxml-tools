# exporters.py

from .models import OrthologGroup, ParalogGroup


def get_ortho_pairs_recursive(node) -> list[(str, str)]:
    """
    Recursively traverse the tree and return a tuple:
      (all_gene_refs_in_subtree, valid_pairs)
    where valid_pairs is a list of tuple pairs (r, s) of geneRefs for which
    the lowest common ancestor is not a ParalogGroup.
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