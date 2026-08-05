from pathlib import Path
from types import SimpleNamespace

from orthoxml.cli import handle_taxonomy
from orthoxml.legacy.models import Taxon


def test_taxon_to_nhx():
    root = Taxon(
        id="5",
        name="Root",
        children=[
            Taxon(id="3", name="Mus musculus"),
            Taxon(
                id="4",
                name="Primates",
                children=[
                    Taxon(id="1", name="Homo sapiens"),
                    Taxon(id="2", name="Pan troglodytes"),
                ],
            ),
        ],
    )

    expected = (
        "(Mus_musculus[&&NHX:S=Mus musculus],"
        "(Homo_sapiens[&&NHX:S=Homo sapiens],Pan_troglodytes[&&NHX:S=Pan troglodytes])"
        "Primates[&&NHX:T=Primates])Root[&&NHX:T=Root];"
    )
    assert root.to_nhx() == expected


def test_handle_taxonomy_writes_nhx(tmp_path, capsys):
    infile = Path(__file__).resolve().parent.parent / "examples" / "data" / "ex1-int-taxon.orthoxml"
    outfile = tmp_path / "taxonomy.nhx"
    args = SimpleNamespace(infile=str(infile), outfile=str(outfile))

    handle_taxonomy(args)

    captured = capsys.readouterr()
    assert "Root" in captured.out
    assert outfile.exists()
    nhx_text = outfile.read_text(encoding="utf-8").strip()
    assert nhx_text.endswith(";")
    assert "[&&NHX:T=Root]" in nhx_text
    assert "[&&NHX:S=Mus musculus]" in nhx_text
