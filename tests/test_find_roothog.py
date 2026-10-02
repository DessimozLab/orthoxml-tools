from pathlib import Path

from orthoxml.custom_parsers import FindRootHOG

EXAMPLES = Path(__file__).parent.parent / "examples" / "data"
FASTOMA = Path(__file__).parent / "test-data" / "fastoma-0.3-roothogs.orthoxml"
SUBSET = Path(__file__).parent / "test-data" / "sample-for-subset.orthoxml"
MULTI_RHOGS = EXAMPLES / "ex4-int-taxon-multiple-rhogs.orthoxml"


def find(path, queries, id="id"):
    with FindRootHOG(path, queries, id=id) as parser:
        parser.parse_through()
    return parser, {r["query"]: r for r in parser.results}


def test_default_query_is_internal_id():
    _, res = find(FASTOMA, ["1018005004", "1000000001"])
    assert res["1018005004"]["roothog_id"] == "HOG_D0637107_sub13680"
    assert res["1018005004"]["roothog_index"] == 1
    assert res["1000000001"]["roothog_id"] == "HOG_D0655833_sub10556"
    assert res["1000000001"]["roothog_index"] == 2


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
    parser, res = find(FASTOMA, ["A0A8M1N6K4"], id="geneId")
    assert res == {}
    assert not parser.id_attr_seen


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
    assert res["8"]["roothog_index"] == 2
    assert res["8"]["num_genes"] == 2


def test_done_after_all_found():
    with FindRootHOG(MULTI_RHOGS, ["1"]) as parser:
        for _ in parser.parse():
            if parser.done:
                break
    assert parser.done
    assert parser.rhog_index == 1
