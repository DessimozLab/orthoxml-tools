# streamfilters.py

import abc
import enum
import logging
import inspect
import xml.etree.ElementTree as ET

from .parsers import StreamOrthoXMLParser, process_stream_orthoxml

logger = logging.getLogger(__name__)


class NodePredicate(abc.ABC):
    """Testable predicate for a node, used in filters."""
    @abc.abstractmethod
    def __call__(self, node) -> bool:
        pass

class CompositePredicate(NodePredicate):
    """Composite predicate that tests all conditions for the node"""
    def __init__(self, predicates: list[NodePredicate]):
        self.predicates = predicates

    def __call__(self, node) -> bool:
        return all(pred(node) for pred in self.predicates)


def node_to_str(node):
    id_field = 0
    ids = node.xpath("@id")
    if len(ids) > 0:
        id_field = ids[0]
    return f"{node.tag} {id_field}"

class ScoreCheck(NodePredicate):
    """
    Checks if the node has a direct child of type:
        <score id="SCORE_ID" value="VALUE">
    and its value >= threshold
    """
    def __init__(self, score_id: str, threshold: float):
        self.score_id = score_id
        self.threshold = threshold

    def __call__(self, node):
        found = False
        for score in node.iterfind(f'./{{http://orthoXML.org/2011/}}score[@id="{self.score_id}"]'):
            found = True
            if float(score.get('value')) < self.threshold:
                return False

        if not found:
            logger.warning(f'Field {self.score_id} not found for node "{node_to_str(node)}". Node remains')
        return True

class GeneNumberCheck(NodePredicate):
    """
    Checks if the node has at least the given number of direct children of type:
        <geneRef ...>
    """
    def __init__(self, min_size: int):
        self.min_size = min_size

    def __call__(self, node):
        size = sum(1 for c in node.iterfind("./{{http://orthoXML.org/2011/}}geneRef"))
        return size >= self.min_size


def enum_to_str(e):
    lower = e.name.lower()
    return lower.replace('_', '-')


class FilterStrategy(enum.Enum):
    CASCADE_REMOVE = 0
    EXTRACT = 1


class CascadeRemoveFilter(StreamOrthoXMLParser):
    """
    Cascade remove (top-down) filtering streamer class

    This stream filter class trims sub-orthologGroups in a top-down fashion.
    Whenever a (sub) orthologGroup element should be removed, the entire
    subtree will be removed.
    """
    def __init__(self, source, predicate: NodePredicate):
        super().__init__(source)
        self.predicate = predicate

    def process_toplevel_group(self, elem):
        # check if rootnode needs to be removed
        if not self.predicate(elem):
            # root element need to be removed, don't return anything
            return None

        to_rem = []
        for hog in elem.iterfind(f'.//{{{self._ns}}}orthologGroup'):
            if not self.predicate(hog):
                to_rem.append(hog)

        logger.info(f"will remove {len(to_rem)} hogs")
        pos_childs = tuple(f"{{{self._ns}}}{z}" for z in ('orthologGroup', 'paralogGroup', 'geneRef'))
        for h in to_rem:
            parent = h.getparent()
            parent.remove(h)
            if sum(c.tag in pos_childs for c in parent) == 0:
                logger.info("removing also empty hog {}".format(parent))
                if parent == elem:
                    return None
                to_rem.append(parent)
        return elem


class ExtractFilter(StreamOrthoXMLParser):
    """
    Extract ("bottom-up") filtering streamer class

    This stream filter class trims sub-orthologGroups in a bottom-up fashion.
    Whenever a (sub) orthologGroup element should be removed, that subtree is
    added as a new toplevel orthologGroup if it contains at least min_hog_size
    genes.
    """
    def __init__(self, source, predicate: NodePredicate, min_hog_size):
        super().__init__(source)
        self.predicate = predicate
        self.min_hog_size = min_hog_size

    def process_toplevel_group(self, elem):
        def get_children_groups(node):
            return node.xpath(
                "./ox:orthologGroup | ./ox:paralogGroup | ./ox:geneRef", namespaces={"ox": self._ns}
            )

        def get_hog_size(node):
            return len(node.xpath(".//ox:geneRef", namespaces={"ox": self._ns}))

        def dfs(node):
            to_extract = set()
            for child in get_children_groups(node):
                child_subtrees = dfs(child)
                to_extract = to_extract.union(child_subtrees)

            # Now decide whether to keep the current element
            if self.strip_ns(node.tag) == "orthologGroup" and not self.predicate(node):
                # If not, remove it from the parent and return valid subtrees
                parent = node.getparent()
                parent.remove(node)
                return to_extract

            if self.strip_ns(node.tag) == "paralogGroup":
                children = get_children_groups(node)
                # Paralog group might not exist anymore if its
                # children were filtered out / extracted
                if len(children) < 2:
                    parent = node.getparent()
                    parent.remove(node)

                    # if only one orthologous group left, reparent it
                    if len(children):
                        parent.insert(-1, children[0])

            # if orthologGroup
            else:
                # Promote parent, drop direct children. If more descended nodes were
                # extracted, we keep them
                for child in get_children_groups(node):
                    if child in to_extract:
                        to_extract.remove(child)

                if get_hog_size(node) >= self.min_hog_size:
                    to_extract.add(node)
            return to_extract

        new_roothogs = list(dfs(elem))
        return new_roothogs

def _strategy_to_filterclass(strategy: str) -> StreamOrthoXMLParser:
    to_enum_val = { enum_to_str(e): e for e in FilterStrategy }

    filter_select = {
        FilterStrategy.CASCADE_REMOVE: CascadeRemoveFilter,
        FilterStrategy.EXTRACT: ExtractFilter
    }

    if strategy in to_enum_val:
        strategy_value = to_enum_val[strategy]
        if strategy_value in filter_select:
            return filter_select[strategy_value]

    valid_choices = ", ".join(v for v in to_enum_val.keys())
    raise ValueError(f"Unsupported strategy {strategy}. Choices are: {valid_choices}")


def _prune_hog_by_genes(elem, kept_gene_ids, ns):
    """
    Remove geneRefs whose id is not in kept_gene_ids, then prune groups
    that become empty or invalid (paralogGroup with <2 children) bottom-up.
    Returns the (modified) elem, or None if no content remains.
    """
    group_tags = {f"{{{ns}}}orthologGroup", f"{{{ns}}}paralogGroup"}
    content_tags = group_tags | {f"{{{ns}}}geneRef"}

    for generef in list(elem.iterfind(f".//{{{ns}}}geneRef")):
        if generef.get("id") not in kept_gene_ids:
            generef.getparent().remove(generef)

    def has_content(node):
        return any(c.tag in content_tags for c in node)

    def prune(node):
        for child in list(node):
            if child.tag in group_tags:
                prune(child)
        for child in list(node):
            if child.tag not in group_tags:
                continue
            children = [c for c in child if c.tag in content_tags]
            local_tag = child.tag.split('}', 1)[-1]
            if not children:
                node.remove(child)
            elif local_tag == 'paralogGroup' and len(children) < 2:
                node.remove(child)
                if children:
                    node.append(children[0])

    prune(elem)
    return elem if has_content(elem) else None


class SpeciesFilter(StreamOrthoXMLParser):
    """
    Single-pass streaming filter: keeps only the listed species.
    geneRefs pointing to removed species are pruned from HOGs,
    and groups that become empty (or paralogGroups with <2 children) are removed.
    """
    def __init__(self, source, species_names):
        super().__init__(source)
        self.species_names = set(species_names)
        self._kept_gene_ids = set()

    def process_species(self, elem):
        name = elem.get("name")
        if name not in self.species_names:
            return None
        for gene in elem.iterfind(f".//{{{self._ns}}}gene"):
            self._kept_gene_ids.add(gene.get("id"))
        return elem

    def process_toplevel_group(self, elem):
        return _prune_hog_by_genes(elem, self._kept_gene_ids, self._ns)


class IndexHOGsByIds(StreamOrthoXMLParser):
    """
    First pass for HOG-ID subsetting.
    Finds every gene ID referenced by the requested HOG IDs.
    If a parent and a child HOG ID are both requested, only the parent is kept.
    """
    def __init__(self, source, hog_ids):
        super().__init__(source)
        self.hog_ids = set(hog_ids)
        self.present_genes: set = set()
        self.matched_hog_ids: set = set()

    def process_toplevel_group(self, elem):
        ns = {"ox": self._ns}
        matches = []

        if elem.get("id") in self.hog_ids:
            matches.append(elem)
        for hog_id in self.hog_ids:
            found = elem.find(f'.//ox:orthologGroup[@id="{hog_id}"]', ns)
            if found is not None and found not in matches:
                matches.append(found)

        if not matches:
            return None

        # keep only non-descendants: if A is ancestor of B, drop B
        kept = []
        for m in matches:
            ancestors = set(m.iterancestors())
            if not any(other in ancestors for other in matches if other is not m):
                kept.append(m)

        for m in kept:
            hog_id = m.get("id")
            if hog_id:
                self.matched_hog_ids.add(hog_id)
            for generef in m.iterfind(".//ox:geneRef", ns):
                self.present_genes.add(generef.get("id"))

        return None


class OutputHOGsById(StreamOrthoXMLParser):
    """
    Second pass for HOG-ID subsetting.
    Emits matched sub-HOGs as new root HOGs, retaining only their genes.
    When species_names is given, also filters to those species.
    """
    def __init__(self, source, matched_hog_ids, present_genes, species_names=None):
        super().__init__(source)
        self.matched_hog_ids = set(matched_hog_ids)
        self.present_gene_ids = set(present_genes)
        self.species_names = set(species_names) if species_names is not None else None
        self._emitted_gene_ids: set = set()

    def process_species(self, elem):
        name = elem.get("name")
        if self.species_names is not None and name not in self.species_names:
            return None
        for gene in list(elem.iterfind(f".//{{{self._ns}}}gene")):
            gid = gene.get("id")
            if gid not in self.present_gene_ids:
                gene.getparent().remove(gene)
            else:
                self._emitted_gene_ids.add(gid)
        genes_left = elem.findall(f".//{{{self._ns}}}gene")
        return elem if genes_left else None

    def process_toplevel_group(self, elem):
        ns = {"ox": self._ns}
        results = []

        if elem.get("id") in self.matched_hog_ids:
            results.append(elem)
        else:
            for hog_id in self.matched_hog_ids:
                found = elem.find(f'.//ox:orthologGroup[@id="{hog_id}"]', ns)
                if found is not None:
                    results.append(found)

        if not results:
            return None

        if self.species_names is not None:
            results = [
                r for r in (
                    _prune_hog_by_genes(sub, self._emitted_gene_ids, self._ns)
                    for sub in results
                )
                if r is not None
            ]

        return results if results else None


def subset_orthoxml(source_orthoxml, out, species_names=None, hog_ids=None):
    """
    Produce a subset OrthoXML containing only the specified species and/or HOGs.

    :param source_orthoxml: path to input OrthoXML file
    :param out: path to write the subset OrthoXML
    :param species_names: iterable of species names to keep (None = keep all)
    :param hog_ids: iterable of HOG IDs to extract as new root HOGs (None = keep all)

    HOG IDs may refer to any nesting level. When a parent and a child HOG ID
    are both listed, only the parent is extracted. When both filters are given
    they are applied together (species filter applied inside each extracted HOG).
    """
    if species_names is None and hog_ids is None:
        raise ValueError("Provide at least one of species_names or hog_ids")

    if hog_ids is not None:
        with IndexHOGsByIds(source_orthoxml, hog_ids) as indexer:
            indexer.parse_through()

        process_stream_orthoxml(
            source_orthoxml,
            out,
            parser_cls=OutputHOGsById,
            parser_kwargs={
                "matched_hog_ids": indexer.matched_hog_ids,
                "present_genes": indexer.present_genes,
                "species_names": list(species_names) if species_names is not None else None,
            },
        )
    else:
        process_stream_orthoxml(
            source_orthoxml,
            out,
            parser_cls=SpeciesFilter,
            parser_kwargs={"species_names": list(species_names)},
        )


def filter_kwargs(cls, kwargs):
    """Selects kwargs that the class expects from the given kwargs array"""
    sig = inspect.signature(cls.__init__)
    valid_params = set(sig.parameters.keys())
    valid_params.discard("self")
    return {k: v for k, v in kwargs.items() if k in valid_params}


def filter_hogs(source_orthoxml, out,
                score_threshold: float,
                min_hog_size: int = 2,
                strategy: str = FilterStrategy.CASCADE_REMOVE):
    """Filter hogs according to the given strategy with a given HOGFilter"""

    parser_cls = _strategy_to_filterclass(strategy)
    score_filter = ScoreCheck("CompletenessScore", score_threshold)

    all_kwargs = {"predicate": score_filter, "min_hog_size": min_hog_size}
    parser_kwargs = filter_kwargs(parser_cls, all_kwargs)

    process_stream_orthoxml(source_orthoxml,
                            out,
                            parser_cls=parser_cls,
                            parser_kwargs=parser_kwargs)

