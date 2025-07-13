import csv
from lxml import etree
from ..parsers import OrthoXMLStreamWriter

def convert_csv_to_orthoxml(
    csv_path: str,
    xml_path: str,
    species_metadata: list[dict] = None,
    xmlns="http://orthoXML.org/2011/",
    root_attrib=None
):
    # 1) Read CSV, collect per‐species gene sets and raw OG rows
    with open(csv_path, newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader)
        species_keys = header[1:]

        # auto-generate minimal metadata if none provided
        if species_metadata is None:
            species_metadata = []
            for sp in species_keys:
                species_metadata.append({
                    "name":        sp,
                    "NCBITaxId":   "0",                # placeholder tax-ID
                    "db_name":     "orthogroups_csv",  # generic DB name
                    "db_version":  "1"                 # generic version
                })

        if len(species_keys) != len(species_metadata):
            raise ValueError("CSV header/species_metadata length mismatch")

        species_genes = {k: set() for k in species_keys}
        og_rows: list[tuple[str, list[list[str]]]] = []

        for row in reader:
            og_id = row[0]
            gene_lists = []
            for idx, cell in enumerate(row[1:]):
                genes = [g.strip() for g in cell.split(',') if g.strip()]
                species_genes[species_keys[idx]].update(genes)
                gene_lists.append(genes)
            og_rows.append((og_id, gene_lists))

    # 2) Assign unique numeric IDs to every gene string
    gene2id: dict[str,int] = {}
    next_id = 1
    for sp in species_keys:
        for gene in sorted(species_genes[sp]):
            gene2id[gene] = next_id
            next_id += 1

    # 3) Write the OrthoXML
    root_attrib = root_attrib or {}
    with OrthoXMLStreamWriter(
         xml_path,
         root_tag="orthoXML",
         xmlns=xmlns,
         attrib=root_attrib
    ) as writer:

        # 3a) species section
        for sp_key, meta in zip(species_keys, species_metadata):
            sp_elem = etree.Element(
                "species",
                name=meta["name"],
                NCBITaxId=str(meta["NCBITaxId"])
            )
            db_elem = etree.SubElement(
                sp_elem, "database",
                name=meta["db_name"],
                version=str(meta["db_version"])
            )
            genes_elem = etree.SubElement(db_elem, "genes")
            for gene in sorted(species_genes[sp_key], key=lambda g: gene2id[g]):
                etree.SubElement(
                    genes_elem, "gene",
                    id=str(gene2id[gene]),
                    geneId=gene
                )
            writer.write_element("species", sp_elem)

        # 3b) groups
        writer.write_element("start_groups", None)
        for og_id, gene_lists in og_rows:
            og_elem = etree.Element("orthologGroup", id=og_id)
            flat_genes = []
            for genes in gene_lists:
                flat_genes.extend(genes)
            
            if not flat_genes:
                continue
            
            for gene in flat_genes:
                etree.SubElement(
                    og_elem, "geneRef",
                    id=str(gene2id[gene])
                )

            writer.write_element("orthologGroup", og_elem)
        writer.write_element("end_groups", None)

    print(f"Wrote OrthoXML to {xml_path}")
