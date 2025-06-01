# custom_parsers.py

from .parsers import StreamOrthoXMLParser

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
