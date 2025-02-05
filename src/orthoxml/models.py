# models.py

class Gene:
    def __init__(self, gene_id: str, gene_name: str, species: str):
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.species = species

    def __repr__(self):
        return f"<Gene: {self.gene_id} ({self.gene_name})>"
    
    def __str__(self):
        return f"Gene: {self.gene_id} ({self.gene_name})"
    
class Species:
    def __init__(self, species_id: str, species_name: str):
        self.species_id = species_id
        self.species_name = species_name

    def __repr__(self):
        return f"<Species: {self.species_id} ({self.species_name})>"
    
    def __str__(self):
        return f"Species: {self.species_id} ({self.species_name})"
    
class OrthologyGroup:
    def __init__(self, group_id: str, genes: list[Gene]):
        self.group_id = group_id
        self.genes = genes

    def __repr__(self):
        return f"<OrthologyGroup: {self.group_id} ({len(self.genes)} genes)>"
    
    def __str__(self):
        return f"OrthologyGroup: {self.group_id} ({len(self.genes)} genes)"
