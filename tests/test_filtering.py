from orthoxml.legacy.tree import OrthoXMLTree

otree = OrthoXMLTree.from_file(
    "tests/test-data/case_filtering.orthoxml",
    score_threshold=0.75,
    score_id="CompletenessScore",
    high_child_as_rhogs=True,
)

otree