from io import BytesIO
from pathlib import Path

import lxml.etree as etree
import pytest

from orthoxml.streamfilters import subset_orthoxml

SAMPLE = Path(__file__).parent / "test-data" / "sample-for-subset.orthoxml"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_output(buf):
    buf.seek(0)
    return etree.parse(buf, etree.XMLParser(remove_blank_text=True))


def species_names(tree):
    ns = {"ox": "http://orthoXML.org/2011/"}
    return {e.get("name") for e in tree.findall(".//ox:species", ns)}


def gene_ids_in_species(tree, sp_name):
    ns = {"ox": "http://orthoXML.org/2011/"}
    for sp in tree.findall(".//ox:species", ns):
        if sp.get("name") == sp_name:
            return {g.get("id") for g in sp.findall(".//ox:gene", ns)}
    return set()


def root_hog_ids(tree):
    ns = {"ox": "http://orthoXML.org/2011/"}
    groups = tree.find(".//ox:groups", ns)
    if groups is None:
        return set()
    return {g.get("id") for g in groups if etree.QName(g.tag).localname == "orthologGroup"}


def all_hog_ids(tree):
    ns = {"ox": "http://orthoXML.org/2011/"}
    return {e.get("id") for e in tree.findall(".//ox:orthologGroup", ns) if e.get("id")}


def gene_refs_in_groups(tree):
    ns = {"ox": "http://orthoXML.org/2011/"}
    return {e.get("id") for e in tree.findall(".//ox:geneRef", ns)}


def declared_gene_ids(tree):
    ns = {"ox": "http://orthoXML.org/2011/"}
    return {e.get("id") for e in tree.findall(".//ox:gene", ns)}


# ---------------------------------------------------------------------------
# Species filter tests
# ---------------------------------------------------------------------------

class TestSpeciesFilter:
    def test_keeps_only_requested_species(self):
        buf = BytesIO()
        subset_orthoxml(SAMPLE, buf, species_names=["Homo sapiens", "Mus musculus"])
        tree = parse_output(buf)
        assert species_names(tree) == {"Homo sapiens", "Mus musculus"}

    def test_removes_genes_of_excluded_species_from_groups(self):
        buf = BytesIO()
        subset_orthoxml(SAMPLE, buf, species_names=["Homo sapiens"])
        tree = parse_output(buf)
        kept = declared_gene_ids(tree)
        refs = gene_refs_in_groups(tree)
        # every geneRef must point to a declared gene
        assert refs <= kept
        # no ref should point to genes of other species
        assert not refs & {"mm1", "mm2", "rn1", "rn2", "gg1", "gg2", "cm1", "cm2", "at1", "at2", "sc1", "sc2"}

    def test_single_species_removes_empty_hogs(self):
        buf = BytesIO()
        # Arabidopsis is only in HOG_Viridiplantae / HOG_Pentapetalae
        subset_orthoxml(SAMPLE, buf, species_names=["Arabidopsis thaliana"])
        tree = parse_output(buf)
        hogs = all_hog_ids(tree)
        # Opistokonta branch should be gone entirely
        assert "HOG_Opistokonta" not in hogs
        assert "HOG_Fungi" not in hogs

    def test_produces_valid_orthoxml_no_dangling_generefs(self):
        buf = BytesIO()
        subset_orthoxml(SAMPLE, buf, species_names=["Gallus gallus", "Chelonia mydas"])
        tree = parse_output(buf)
        declared = declared_gene_ids(tree)
        refs = gene_refs_in_groups(tree)
        assert refs <= declared

    def test_all_species_is_identity(self):
        all_sp = [
            "Homo sapiens", "Mus musculus", "Rattus norvegicus",
            "Gallus gallus", "Chelonia mydas", "Arabidopsis thaliana",
            "Saccharomyces cerevisiae",
        ]
        buf = BytesIO()
        subset_orthoxml(SAMPLE, buf, species_names=all_sp)
        tree = parse_output(buf)
        assert species_names(tree) == set(all_sp)
        assert gene_refs_in_groups(tree) == declared_gene_ids(tree)


# ---------------------------------------------------------------------------
# HOG ID filter tests
# ---------------------------------------------------------------------------

class TestHOGIdFilter:
    def test_extract_root_hog(self):
        buf = BytesIO()
        subset_orthoxml(SAMPLE, buf, hog_ids=["HOG_Eukaryota"])
        tree = parse_output(buf)
        assert root_hog_ids(tree) == {"HOG_Eukaryota"}
        # all species should be present
        assert species_names(tree) == {
            "Homo sapiens", "Mus musculus", "Rattus norvegicus",
            "Gallus gallus", "Chelonia mydas", "Arabidopsis thaliana",
            "Saccharomyces cerevisiae",
        }

    def test_extract_nested_hog_becomes_root(self):
        buf = BytesIO()
        subset_orthoxml(SAMPLE, buf, hog_ids=["HOG_Opistokonta"])
        tree = parse_output(buf)
        assert root_hog_ids(tree) == {"HOG_Opistokonta"}
        # Arabidopsis is NOT in Opistokonta
        assert "Arabidopsis thaliana" not in species_names(tree)

    def test_extract_deeply_nested_hog(self):
        buf = BytesIO()
        subset_orthoxml(SAMPLE, buf, hog_ids=["HOG_Sauria"])
        tree = parse_output(buf)
        assert root_hog_ids(tree) == {"HOG_Sauria"}
        # only Gallus + Chelonia should appear
        assert species_names(tree) == {"Gallus gallus", "Chelonia mydas"}

    def test_extract_two_independent_hogs(self):
        buf = BytesIO()
        subset_orthoxml(SAMPLE, buf, hog_ids=["HOG_Sauria", "HOG_Viridiplantae"])
        tree = parse_output(buf)
        assert root_hog_ids(tree) == {"HOG_Sauria", "HOG_Viridiplantae"}
        assert species_names(tree) == {"Gallus gallus", "Chelonia mydas", "Arabidopsis thaliana"}

    def test_parent_child_dedup_keeps_parent(self):
        """If both a parent and a child HOG ID are requested, only the parent is kept."""
        buf = BytesIO()
        subset_orthoxml(SAMPLE, buf, hog_ids=["HOG_Opistokonta", "HOG_Sauria"])
        tree = parse_output(buf)
        # HOG_Sauria is inside HOG_Opistokonta → only one root HOG
        assert root_hog_ids(tree) == {"HOG_Opistokonta"}

    def test_unknown_hog_id_produces_empty_groups(self):
        buf = BytesIO()
        subset_orthoxml(SAMPLE, buf, hog_ids=["HOG_DOES_NOT_EXIST"])
        tree = parse_output(buf)
        assert gene_refs_in_groups(tree) == set()

    def test_no_dangling_generefs(self):
        buf = BytesIO()
        subset_orthoxml(SAMPLE, buf, hog_ids=["HOG_Tetrapoda"])
        tree = parse_output(buf)
        declared = declared_gene_ids(tree)
        refs = gene_refs_in_groups(tree)
        assert refs <= declared


# ---------------------------------------------------------------------------
# Combined HOG ID + species filter tests
# ---------------------------------------------------------------------------

class TestCombinedFilter:
    def test_hog_and_species_combined(self):
        buf = BytesIO()
        subset_orthoxml(
            SAMPLE, buf,
            hog_ids=["HOG_Opistokonta"],
            species_names=["Homo sapiens", "Mus musculus"],
        )
        tree = parse_output(buf)
        assert root_hog_ids(tree) == {"HOG_Opistokonta"}
        assert species_names(tree) == {"Homo sapiens", "Mus musculus"}
        # reptile/bird/yeast refs should be gone
        refs = gene_refs_in_groups(tree)
        assert not refs & {"gg1", "gg2", "cm1", "cm2", "sc1", "sc2"}

    def test_combined_removes_now_empty_sub_hogs(self):
        buf = BytesIO()
        # Sauria has only Gallus + Chelonia; filtering to just Homo sapiens
        # should eliminate Sauria branch entirely
        subset_orthoxml(
            SAMPLE, buf,
            hog_ids=["HOG_Tetrapoda"],
            species_names=["Homo sapiens"],
        )
        tree = parse_output(buf)
        assert "HOG_Sauria" not in all_hog_ids(tree)
        assert "Homo sapiens" in species_names(tree)

    def test_combined_no_dangling_generefs(self):
        buf = BytesIO()
        subset_orthoxml(
            SAMPLE, buf,
            hog_ids=["HOG_Eukaryota"],
            species_names=["Homo sapiens", "Arabidopsis thaliana"],
        )
        tree = parse_output(buf)
        declared = declared_gene_ids(tree)
        refs = gene_refs_in_groups(tree)
        assert refs <= declared
