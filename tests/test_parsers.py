from pathlib import Path
from io import BytesIO
import lxml.etree as etree
import pytest

from orthoxml.parsers import process_stream_orthoxml


###########################################
# Helper functions for comparing xml files
###########################################
def remove_all_comments(tree: etree._ElementTree) -> None:
    comments = tree.xpath('//comment()')
    for comment in comments:
        parent = comment.getparent()
        if parent is not None:
            parent.remove(comment)


def canonicalize_xml(xml_path):
    """this function returns the canonical form of an input xml file as a string,
    such that two xml files can be compared on the content. Whitespace differences
    are ignored"""
    close = False
    if isinstance(xml_path, str):
        f = open(xml_path, 'rb')
        close = True
    else:
        f = xml_path

    parser = etree.XMLParser(remove_blank_text=True)
    tree = etree.parse(f, parser)
    remove_all_comments(tree)
    if close:
        f.close()
    return etree.tostring(tree, method="c14n")  # Canonical form


################################################
# pytest fixtures
################################################

@pytest.fixture
def single_group_example_orthoxml():
    return Path(__file__).parent / "test-data" / "case_filtering.orthoxml"

@pytest.fixture
def multiple_group_example_orthoxml():
    return Path(__file__).parent / "test-data" / "multiple_groups.orthoxml"


def test_process_stream_round_robin(single_group_example_orthoxml):
    """ensure reading and writing of orthoxml works results in the same content"""
    out_stream = BytesIO()
    process_stream_orthoxml(single_group_example_orthoxml, out_stream)
    out_stream.seek(0)
    assert canonicalize_xml(out_stream) == canonicalize_xml(single_group_example_orthoxml)

