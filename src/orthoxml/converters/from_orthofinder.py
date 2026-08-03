import csv
from lxml import etree
from ..parsers import OrthoXMLStreamWriter


def convert_orthoxml_to_csv(
    xml_path: str,
    csv_path: str,
    id_attr: str = "geneId",
    delimiter: str = "\t",
):
    """Export a simplified OrthoXML file to the OrthoFinder-style TSV/CSV format.

    The output is a tab-delimited table with one row per top-level ortholog group.
    The first column is the group id, followed by one column per species. Each cell
    contains a comma-separated list of gene identifiers for that species in that group.
    """

    tree = etree.parse(xml_path)
    root = tree.getroot()
    ns = root.nsmap.get(None, "")
    ns_tag = f"{{{ns}}}" if ns else ""

    species_nodes = root.findall(f"{ns_tag}species")
    species_names = []
    gene_id_to_label = {}
    gene_id_to_species = {}

    for species in species_nodes:
        species_name = species.get("name")
        species_names.append(species_name)

        for gene in species.findall(f".//{ns_tag}gene"):
            gene_internal_id = gene.get("id")
            gene_label = gene.get(id_attr) or gene.get("id")
            if not gene_internal_id:
                continue
            gene_id_to_label[gene_internal_id] = gene_label
            gene_id_to_species[gene_internal_id] = species_name

    group_rows = []
    groups_node = root.find(f"{ns_tag}groups")
    if groups_node is not None:
        for group in groups_node.findall(f"{ns_tag}orthologGroup"):
            genes_by_species = {name: [] for name in species_names}
            for gene_ref in group.findall(f".//{ns_tag}geneRef"):
                gene_internal_id = gene_ref.get("id")
                if not gene_internal_id:
                    continue
                species_name = gene_id_to_species.get(gene_internal_id)
                if species_name is None:
                    continue
                genes_by_species[species_name].append(gene_id_to_label.get(gene_internal_id, gene_internal_id))

            group_rows.append([
                group.get("id") or "",
                *[
                    ", ".join(genes_by_species[species_name])
                    for species_name in species_names
                ],
            ])

    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter=delimiter, lineterminator="\n")
        writer.writerow(["Orthogroup", *species_names])
        writer.writerows(group_rows)

    print(f"Wrote CSV/TSV to {csv_path}")


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
