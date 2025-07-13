# custom_parsers.py
import enum
from collections import defaultdict
from .parsers import StreamOrthoXMLParser
from .logger import get_logger
from .models import Taxon, ORTHO_NS
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
