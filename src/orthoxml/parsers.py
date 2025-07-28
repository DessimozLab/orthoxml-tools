# parsers.py
import os.path
from logging import getLogger
from os import PathLike
from typing import Union, IO, Any, Optional

from lxml import etree

from .utils import auto_open

logger = getLogger(__name__)


def strip_namespace(element):
    """Remove namespace from the given element and its descendants."""
    new_element = etree.Element(etree.QName(element).localname, nsmap=None)
    # Copy attributes
    new_element.attrib.update(element.attrib)
    # Recursively copy children
    for child in element:
        new_child = child if isinstance(child, etree._Comment) else strip_namespace(child)
        new_element.append(new_child)

    # If element has text, add it
    if element.text:
        new_element.text = element.text
    return new_element


class OrthoXMLStreamWriter:
    def __init__(self, target: Union[str, PathLike, IO[bytes]], root_tag="orthoXML", xmlns=None, attrib=None, **kwargs):
        self._target = target
        self.root_tag = root_tag
        self.xmlns = xmlns or ''
        self.attrib = attrib or {}
        self._depth = 0

    def __enter__(self):
        self._should_close = False
        if isinstance(self._target, (str, PathLike)):
            self.stream = auto_open(self._target, "wb")
            self._should_close = True
        else:
            self.stream = self._target
        self._write_header()
        return self

    def _write_header(self):
        self.stream.write(b'<?xml version="1.0" encoding="UTF-8"?>\n')
        self.stream.write(b"<%s" % self.root_tag.encode("utf-8"))
        if self.xmlns:
            self.stream.write(b' xmlns="%s"' % self.xmlns.encode("utf-8"))
        for key, value in self.attrib.items():
            self.stream.write(b' %s="%s"' % (key.encode("utf-8"), value.encode("utf-8")))
        self.stream.write(b">\n")
        self._depth += 1

    def write_element(self, tag: str, elem: Optional[etree._Element]):
        if tag == "start_groups":
            self.stream.write(b"  " * self._depth)
            self.stream.write(b"<groups>\n")
            self._depth += 1
        elif tag == "end_groups":
            self._depth -= 1
            self.stream.write(b"  " * self._depth)
            self.stream.write(b"</groups>\n")
        elif elem is not None and isinstance(elem, etree._Element):
            elem = strip_namespace(elem)
            etree.indent(elem, level=self._depth)
            self.stream.write(b"  " * self._depth)
            self.stream.write(etree.tostring(elem, encoding="utf-8"))
            self.stream.write(b"\n")
        else:
            raise ValueError("Unsupported data type")
        self.stream.flush()

    def close(self):
        self.stream.write(f"</{self.root_tag}>\n".encode("utf-8"))
        if self._should_close:
            self.stream.close()

    def __exit__(self, *exc):
        self.close()


class StreamOrthoXMLParser:
    def __init__(self, source: Union[str, PathLike, IO[Any]], **kwargs):
        self.source = source

    def strip_ns(self, tag):
        return tag.split('}', 1)[-1]

    def __enter__(self):
        self._should_close = False
        if isinstance(self.source, (str, PathLike)):
            self.stream = auto_open(self.source, "rb")
            self._should_close = True
        else:
            self.stream = self.source

        self._context = etree.iterparse(self.stream, events=('start', 'end'))
        _, root = next(self._context)
        self._ns, tag = root.tag[1:].split('}')
        self.nsmap = {'': self._ns}
        self.root_tag = tag
        self.root_attribs = dict(root.attrib)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._should_close:
            self.stream.close()

    def parse_through(self):
        for _ in self.parse():
            pass

    def parse(self):
        og_depth, pg_depth = 0, 0
        for event, elem in self._context:
            tag = self.strip_ns(elem.tag)

            if event == 'end':
                if tag == "species":
                    if (result := self.process_species(elem)) is not None:
                        yield ("species", result)
                    self._cleanup(elem)
                elif tag == "taxonomy":
                    if (result := self.process_taxonomy(elem)) is not None:
                        yield ("taxonomy", result)
                    self._cleanup(elem)
                elif tag == "scores":
                    if (result := self.process_scores(elem)) is not None:
                        yield ("scores", result)
                    self._cleanup(elem)
                elif tag == "orthologGroup":
                    og_depth -= 1
                    if og_depth == 0 and pg_depth == 0:
                        if (result := self.process_toplevel_group(elem)) is not None:
                            if isinstance(result, list):
                                yield from (("orthologGroup", r) for r in result)
                            else:
                                yield ("orthologGroup", result)
                        self._cleanup(elem)
                elif tag == "paralogGroup":
                    pg_depth -= 1
                    if og_depth == 0 and pg_depth == 0:
                        if (result := self.process_toplevel_group(elem)) is not None:
                            yield ("paralogGroup", result)
                        self._cleanup(elem)
                elif tag == "groups":
                    yield ("end_groups", None)
                elif tag == "notes":
                    if (result := self.process_notes(elem)) is not None:
                        yield ("notes", result)
                    self._cleanup(elem)

            elif event == 'start':
                if tag == "groups":
                    yield ("start_groups", None)
                elif tag == "orthologGroup":
                    og_depth += 1
                elif tag == "paralogGroup":
                    pg_depth += 1

    def _cleanup(self, elem):
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]

    # ---- Default hooks to override ----
    def process_species(self, elem): return elem
    def process_taxonomy(self, elem): return elem
    def process_scores(self, elem): return elem
    def process_toplevel_group(self, elem): return elem
    def process_notes(self, elem): return elem


def process_stream_orthoxml(
        infile: Union[str, PathLike, IO],
        outfile: Union[str, PathLike, IO],
        parser_cls=StreamOrthoXMLParser,
        writer_cls=OrthoXMLStreamWriter,
        parser_kwargs=dict(),
        writer_kwargs=dict(),
):
    # only mkdir if outfile is actually a path
    if isinstance(outfile, (str, PathLike)):
        outpath = str(outfile)
        parent_dir = os.path.dirname(outpath)
        if parent_dir:
            os.makedirs(parent_dir, exist_ok=True)

    with parser_cls(infile, **parser_kwargs) as parser:
        root_tag, nsmap, attrib = parser.root_tag, parser.nsmap, parser.root_attribs
        with writer_cls(outfile, root_tag=root_tag, xmlns=nsmap[''], attrib=attrib, **writer_kwargs) as writer:
            for tag, elem in parser.parse():
                writer.write_element(tag, elem)
