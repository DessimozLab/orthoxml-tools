from pathlib import Path

import lxml.etree as etree

from orthoxml.custom_parsers import OrthoXMLSplitter


def test_splitter_only_yields_ortholog_or_paralog_groups(tmp_path):
    xml = '''<?xml version="1.0" encoding="UTF-8"?>
<orthoXML xmlns="http://orthoXML.org/2011/">
  <species name="s1">
    <database name="db1">
      <genes>
        <gene id="g1"/>
      </genes>
    </database>
  </species>
  <groups>
    <notAGroup id="bad"/>
    <orthologGroup id="HOG_1">
      <geneRef id="g1"/>
    </orthologGroup>
    <paralogGroup id="HOG_2">
      <geneRef id="g1"/>
    </paralogGroup>
  </groups>
</orthoXML>
'''
    infile = tmp_path / "mini.orthoxml"
    infile.write_text(xml)
    outdir = tmp_path / "out"
    outdir.mkdir()

    splitter = OrthoXMLSplitter(str(infile), cache_dir=str(outdir))
    yielded = [node.get("id") for node in splitter._iter_toplevel_groups()]

    assert yielded == ["HOG_1", "HOG_2"]


def test_splitter_creates_one_file_per_selected_root_hog(tmp_path):
    xml = '''<?xml version="1.0" encoding="UTF-8"?>
<orthoXML xmlns="http://orthoXML.org/2011/">
  <species name="s1">
    <database name="db1">
      <genes>
        <gene id="g1"/>
        <gene id="g2"/>
      </genes>
    </database>
  </species>
  <groups>
    <orthologGroup id="HOG_1">
      <geneRef id="g1"/>
    </orthologGroup>
    <orthologGroup id="HOG_2">
      <geneRef id="g2"/>
    </orthologGroup>
  </groups>
</orthoXML>
'''
    infile = tmp_path / "mini.orthoxml"
    infile.write_text(xml)
    outdir = tmp_path / "split-output"
    outdir.mkdir()

    splitter = OrthoXMLSplitter(str(infile), cache_dir=str(outdir))
    splitter(hogs_to_extract={"HOG_1"}, single_hog_files=False, basename="subset.orthoxml")

    out_file = outdir / "HOG_1.orthoxml"
    assert out_file.exists()

    tree = etree.parse(str(out_file), etree.XMLParser(remove_blank_text=True))
    groups = tree.findall("{http://orthoXML.org/2011/}groups/{http://orthoXML.org/2011/}orthologGroup")
    assert [g.get("id") for g in groups] == ["HOG_1"]


def test_splitter_creates_single_bundle_file_for_single_hog_mode(tmp_path):
    xml = '''<?xml version="1.0" encoding="UTF-8"?>
<orthoXML xmlns="http://orthoXML.org/2011/">
  <species name="s1">
    <database name="db1">
      <genes>
        <gene id="g1"/>
      </genes>
    </database>
  </species>
  <groups>
    <paralogGroup id="HOG_9">
      <geneRef id="g1"/>
    </paralogGroup>
  </groups>
</orthoXML>
'''
    infile = tmp_path / "mini.orthoxml"
    infile.write_text(xml)
    outdir = tmp_path / "split-output"
    outdir.mkdir()

    splitter = OrthoXMLSplitter(str(infile), cache_dir=str(outdir))
    splitter(hogs_to_extract={"HOG_9"}, single_hog_files=True, basename="subset.orthoxml")

    out_file = outdir / "subset.orthoxml"
    assert out_file.exists()

    tree = etree.parse(str(out_file), etree.XMLParser(remove_blank_text=True))
    groups = tree.findall("{http://orthoXML.org/2011/}groups/{http://orthoXML.org/2011/}paralogGroup")
    assert [g.get("id") for g in groups] == ["HOG_9"]
