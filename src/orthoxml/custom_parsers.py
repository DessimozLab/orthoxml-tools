# custom_parsers.py

from collections import defaultdict
from .parsers import StreamOrthoXMLParser
from .logger import get_logger
from .legacy.models import Taxon, ORTHO_NS
from lxml import etree

logger = get_logger(__name__)

class BasicStats(StreamOrthoXMLParser):
    def __init__(self, source):
        super().__init__(source)
        self.gene_count = 0
        self.rhog_count = 0
        self.species_count = 0
        self.leave_taxon_count = 0
        self.all_taxa_count = 0

    def process_species(self, elem):
        """Count how many species and genes we have in the orthoxml file"""

        self.species_count += 1

        gene_tag = f"{{{self._ns}}}gene"
        genes_in_this_species = elem.findall(f".//{gene_tag}")
        num_genes = len(genes_in_this_species)
        self.gene_count += num_genes

        return None
    
    def process_taxonomy(self, elem):
        """Count how many leave taxon we have in the taxonomy"""

        taxon_tag = f"{{{self._ns}}}taxon"
        all_taxa = elem.findall(f".//{taxon_tag}")
        self.all_taxa_count = len(all_taxa)
        
        count = 0
        for taxon in all_taxa:
            has_child_taxon = any(child.tag == taxon_tag for child in taxon)
            if not has_child_taxon:
                count += 1
        self.leave_taxon_count = count

        return None

    def process_scores(self, elem):
        return None
    def process_toplevel_group(self, elem):
        self.rhog_count += 1
        return None

class GenePerTaxonStats(StreamOrthoXMLParser):
    def __init__(self, source):
        super().__init__(source)
        self.gene_count_per_taxon = defaultdict(int)
        self.header_gene_count_per_species = {}
        self.gene_to_species_name = {}
        self.taxonomy_counts = {}
        self.taxonomy_tree = None

    def process_species(self, elem):
        """Count how many genes we have per species in the orthoxml file"""

        species_name = elem.get("name")

        gene_tag = f"{{{self._ns}}}gene"
        genes_in_this_species = elem.findall(f".//{gene_tag}")
        num_genes = len(genes_in_this_species)

        self.header_gene_count_per_species[species_name] = num_genes

        for gene in genes_in_this_species:
            gene_id = gene.get("id")
            self.gene_to_species_name[gene_id] = species_name

        return None
    
    def process_toplevel_group(self, elem):
        """
        Called once for each top-level <orthologGroup> or <paralogGroup>.
        Count all geneRef's per species under this group.
        """
        gene_ref_tag = f"{{{self._ns}}}geneRef"

        # find every geneRef anywhere inside this group
        for gr in elem.findall(f".//{gene_ref_tag}"):
            gid = gr.get("id")
            species = self.gene_to_species_name.get(gid)
            if not species:
                logger.warning(
                    f"GeneRef with id '{gid}' not found in species mapping. "
                    "This may indicate a mismatch in gene IDs between header and groups."
                )            
                continue

            # accumulate into the global tally
            self.gene_count_per_taxon[species] = (
                self.gene_count_per_taxon.get(species, 0) + 1
            )

        return None

    def process_taxonomy(self, elem):
            """Build an in‐memory tree of nested <taxon> elements."""
            taxon_tag = f"{{{self._ns}}}taxon"
            def build_node(tx_elem):
                return {
                    "id":      tx_elem.get("id"),
                    "name":    tx_elem.get("name"),
                    "children":[ build_node(c) 
                                for c in tx_elem 
                                if isinstance(c, etree._Element) and c.tag==taxon_tag ]
                }

            roots = [ build_node(c) for c in elem 
                    if isinstance(c, etree._Element) and c.tag==taxon_tag ]

            if len(roots)==1:
                self.taxonomy_tree = roots[0]
            else:
                self.taxonomy_tree = {"id":None, "name":"<root>", "children":roots}
            return None

    def compute_taxon_counts(self):
        """Walk the taxonomy_tree and sum up gene_count_per_taxon into every node."""
        def recurse(node):
            if not node["children"]:
                cnt = self.gene_count_per_taxon.get(node["name"], 0)
            else:
                cnt = sum(recurse(ch) for ch in node["children"])
            self.taxonomy_counts[node["name"]] = cnt
            return cnt

        if self.taxonomy_tree is None:
            logger.warning("No taxonomy tree found. Cannot compute taxon counts.")
            return 0
        recurse(self.taxonomy_tree)

class PrintTaxonomy(StreamOrthoXMLParser):
    def __init__(self, source):
        super().__init__(source)
        self.taxonomy = None

    def process_taxonomy(self, elem):
        """Build an in‐memory tree of nested <taxon> elements."""

        if elem is not None:
            taxon_el = elem.find(f"{{{ORTHO_NS}}}taxon")
            if taxon_el is not None:
                self.taxonomy = Taxon.from_xml(taxon_el)

        return None


class RootHOGCounter(StreamOrthoXMLParser):
    def __init__(self, source, **kwargs):
        super().__init__(source, **kwargs)
        self.rhogs_count = 0

    def process_toplevel_group(self, elem):
        self.rhogs_count += 1

        return None

class SplitterByRootHOGS(StreamOrthoXMLParser):
    def __init__(self, source, rhogs_number):
        super().__init__(source)
        self.rhogs_number = rhogs_number
        self.current_rhog = 0

    def process_toplevel_group(self, elem):
        self.current_rhog += 1

        if self.current_rhog == self.rhogs_number:
            return elem

class GetGene2IdMapping(StreamOrthoXMLParser):
    """Get the mapping between id and geneId or protId, ..."""
    def __init__(self, source, id):
        super().__init__(source)
        self.gene_id2id_mapping = {}
        self.id = id

    def process_species(self, elem):
        gene_tag = f"{{{self._ns}}}gene"
        genes_in_this_species = elem.findall(f".//{gene_tag}")

        for gene in genes_in_this_species:
            self.gene_id2id_mapping[gene.attrib.get("id")] = gene.attrib.get(self.id)

        return None


class StreamPairsParser(StreamOrthoXMLParser):
    """
    Extends StreamOrthoXMLParser with a streaming ortholog or para-pair extractor.
    """
    def __init__(self, source, ortho_para):
        super().__init__(source)
        self.ortho_para = ortho_para # orthologGroup or paralogGroup

    def iter_pairs(self):
        """
        Yield (r_id, s_id) for every ortholog pair in the file,
        using only O(tree-depth) memory.
        """
        # each frame: { type, own_refs, child_refs, child_pairs }
        group_stack = []

        for event, elem in self._context:
            tag = self.strip_ns(elem.tag)

            # 1) on group start, push a fresh frame
            if event == 'start' and tag in ('orthologGroup','paralogGroup'):
                group_stack.append({
                    "type":        tag,
                    "own_refs":    [],
                    "child_refs":  [],
                    "child_pairs": []
                })

            # 2) on geneRef end, stash its id
            elif event == 'end' and tag == 'geneRef':
                if group_stack:
                    group_stack[-1]["own_refs"].append(elem.get("id"))

            # 3) on group end, pop & process
            elif event == 'end' and tag in ('orthologGroup','paralogGroup'):
                frame   = group_stack.pop()
                own     = frame["own_refs"]
                cref_l  = frame["child_refs"]
                cpair_l = frame["child_pairs"]

                # build the flat list of all gene-refs under this node
                gene_refs = own.copy()
                for cr in cref_l:
                    gene_refs.extend(cr)

                # flatten child-computed pairs so we can pass them up
                flat_child_pairs = [p for cp in cpair_l for p in cp]

                # if this is an orthologGroup, compute *new* pairs here
                new_pairs = []
                if frame["type"] == self.ortho_para:
                    # (a) all‐own‐refs
                    for i in range(len(own)):
                        for j in range(i+1, len(own)):
                            new_pairs.append((own[i], own[j]))
                    # (b) own‐vs‐each‐child
                    for cr in cref_l:
                        for r in own:
                            for s in cr:
                                new_pairs.append((r, s))
                    # (c) between‐different‐children
                    for i in range(len(cref_l)):
                        for j in range(i+1, len(cref_l)):
                            for r in cref_l[i]:
                                for s in cref_l[j]:
                                    new_pairs.append((r, s))

                    # **only** yield the new pairs for this node
                    for pair in new_pairs:
                        yield pair

                # aggregate pairs to hand up to parent
                pairs_to_pass_up = flat_child_pairs + new_pairs

                # 4) hand results up to parent frame (if any)
                if group_stack:
                    parent = group_stack[-1]
                    parent["child_refs"].append(gene_refs)
                    parent["child_pairs"].append(pairs_to_pass_up)

                # 5) free memory under this element
                elem.clear()
                while elem.getprevious() is not None:
                    del elem.getparent()[0]


class StreamMaxOGParser(StreamOrthoXMLParser):
    def __init__(self, source):
        super().__init__(source)
        # map from geneId (string) → species name
        self.species_map: dict[str,str] = {}

    def process_species(self, elem):
        """Called on </species>: collect gene→species mapping."""
        sp_name = elem.get("name")
        # walk down to <gene> elements
        for gene in elem.findall(".//{%s}gene" % self._ns, namespaces=self.nsmap):
            gid = gene.get("geneId") or gene.get("id")
            if gid:
                self.species_map[gid] = sp_name
        # we return None so nothing is yielded for species
        return None

    def process_toplevel_group(self, elem):
        """
        Called on each top-level <orthologGroup> or <paralogGroup> under <groups>.
        Here `elem` is the root of one OG/PG subtree.
        We compute and return the list of gene IDs to keep.
        """
        def local_strip(tag):
            return tag.split("}",1)[-1]

        def recurse(node) -> list[str]:
            # gather direct geneRef IDs
            direct_refs = [gr.get("id")
                           for gr in node
                           if local_strip(gr.tag) == "geneRef"
                          ]
            # gather all group‐children
            child_groups = [c for c in node
                            if local_strip(c.tag) in ("orthologGroup","paralogGroup")]

            # if there are child groups, process them first
            if child_groups:
                child_kept = [recurse(c) for c in child_groups]

                # Duplication event = a <paralogGroup> that has child groups
                if local_strip(node.tag) == "paralogGroup":
                    # compute species‐counts for each branch
                    counts = []
                    for genes in child_kept:
                        # only count those we know species for
                        sps = { self.species_map[g] 
                                for g in genes 
                                if g in self.species_map }
                        counts.append(len(sps))
                    # pick the branch with the max distinct species
                    idx = counts.index(max(counts))
                    return child_kept[idx]

                else:
                    # an <orthologGroup> with children: union everything + any direct refs
                    out = []
                    for genes in child_kept:
                        out.extend(genes)
                    out.extend(direct_refs)
                    return out

            else:
                # leaf group (no sub-groups)
                if local_strip(node.tag) == "paralogGroup":
                    # keep only the first geneRef in a leaf paralogGroup
                    return direct_refs[:1]
                else:
                    # orthologGroup leaf: keep them all
                    return direct_refs

        # run our bottom-up pass and return its result
        return [recurse(elem)]
