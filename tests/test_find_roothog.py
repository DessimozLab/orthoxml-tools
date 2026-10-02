import gzip
from pathlib import Path

import pytest

from orthoxml import roothog_lookup
from orthoxml.roothog_lookup import RootHOGLookup

EXAMPLES = Path(__file__).parent.parent / "examples" / "data"
FASTOMA = Path(__file__).parent / "test-data" / "fastoma-0.3-roothogs.orthoxml"
SUBSET = Path(__file__).parent / "test-data" / "sample-for-subset.orthoxml"
MULTI_RHOGS = EXAMPLES / "ex4-int-taxon-multiple-rhogs.orthoxml"


def find(path, queries, id="id", **kwargs):
    lookup = RootHOGLookup(path, queries, id=id, **kwargs).run()
    return lookup, {r["query"]: r for r in lookup.results}


def test_default_query_is_internal_id():
    _, res = find(FASTOMA, ["1018005004", "1000000001"])
    assert res["1018005004"]["roothog_id"] == "HOG_D0637107_sub13680"
    assert res["1000000001"]["roothog_id"] == "HOG_D0655833_sub10556"


def test_protid_only_matched_when_requested():
    _, res = find(FASTOMA, ["A0A8M1N6K4"])
    assert res == {}
    _, res = find(FASTOMA, ["A0A8M1N6K4"], id="protId")
    assert res["A0A8M1N6K4"]["gene_id"] == "1000000001"
    assert res["A0A8M1N6K4"]["roothog_id"] == "HOG_D0655833_sub10556"


def test_taxrange_of_roothog_not_nested_group():
    _, res = find(FASTOMA, ["1030000384"])
    assert res["1030000384"]["taxon_level"] == "Euteleostomi"


def test_num_genes_counts_nested_and_paralog_genes():
    _, res = find(FASTOMA, ["1018005004", "1000000001"])
    assert res["1018005004"]["num_genes"] == 3
    assert res["1000000001"]["num_genes"] == 3


def test_missing_attribute_detected():
    lookup, res = find(FASTOMA, ["A0A8M1N6K4"], id="geneId")
    assert res == {}
    assert not lookup.id_attr_seen


def test_missing_query_is_not_reported():
    _, res = find(SUBSET, ["P00001", "DOES_NOT_EXIST"], id="protId")
    assert set(res) == {"P00001"}
    assert res["P00001"]["roothog_id"] == "HOG_Eukaryota"
    assert res["P00001"]["num_genes"] == 14


def test_taxonid_resolved_through_taxonomy():
    _, res = find(MULTI_RHOGS, ["1", "8"])
    assert res["1"]["roothog_id"] is None
    assert res["1"]["taxon_level"] == "Root"
    assert res["1"]["num_genes"] == 6
    assert res["8"]["num_genes"] == 2


@pytest.mark.parametrize("chunk_size", [7, 64, 512])
def test_tiny_chunks_give_same_result(chunk_size):
    queries = ["1000000002", "1030000384", "1018005005", "1000000001"]
    _, expected = find(FASTOMA, queries)
    _, res = find(FASTOMA, queries, chunk_size=chunk_size)
    assert res == expected
    assert len(res) == 4


def test_many_queries_regex_path(monkeypatch):
    monkeypatch.setattr(roothog_lookup, "_FIND_THRESHOLD", 0)
    _, res = find(FASTOMA, ["P12345", "Q22222", "NOPE"], id="protId")
    assert res["P12345"]["roothog_id"] == "HOG_D0637107_sub13680"
    assert res["Q22222"]["roothog_id"] == "HOG_D0655833_sub10556"
    assert "NOPE" not in res


def test_gzipped_input(tmp_path):
    gz = tmp_path / "fastoma.orthoxml.gz"
    gz.write_bytes(gzip.compress(FASTOMA.read_bytes()))
    _, res = find(gz, ["P12345"], id="protId")
    assert res["P12345"]["taxon_level"] == "Euteleostomi"
