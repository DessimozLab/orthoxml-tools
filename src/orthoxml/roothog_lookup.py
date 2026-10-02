# roothog_lookup.py
"""Fast lookup of the rootHOG containing given genes.

Instead of building an XML tree for the whole file, the file is scanned as
raw bytes in large chunks:

1. In the header, the queried attribute (``id``, ``protId``, ``geneId``) is
   searched with ``bytes.find`` and the internal gene ``id`` is read from the
   matching ``<gene>`` tag.
2. In ``<groups>``, ``<geneRef id="...">`` occurrences are searched the same
   way. The nesting depth at a hit is obtained by counting group open/close
   tags (C-speed ``bytes.count``), and the enclosing top-level group is found
   by walking the group tags backwards until depth 0.
3. Only that rootHOG is parsed with lxml to get its id, taxonomic level and
   number of genes.

Assumptions: groups use the default OrthoXML namespace (no tag prefix) and
are not self-closing (``<orthologGroup/>``).
"""
from __future__ import annotations

import re
from os import PathLike
from typing import IO, Iterable, Union

from lxml import etree

from .utils import auto_open

DEFAULT_CHUNK_SIZE = 1 << 25  # 32 MiB

# below this number of queries each one is searched with bytes.find,
# above it a single regex pass with a set lookup is cheaper
_FIND_THRESHOLD = 32

_OPEN_TAGS = (b"<orthologGroup", b"<paralogGroup")
_CLOSE_TAGS = (b"</orthologGroup>", b"</paralogGroup>")
_WS = b" \t\r\n"
_ID_IN_TAG_RE = re.compile(rb'\sid="([^"]*)"')
_GENEREF_ID_RE = re.compile(rb'<geneRef\s+id="([^"]*)"')
_GROUPS_START_RE = re.compile(rb"<groups[\s>]")


def _depth_delta(buf: bytes, start: int, end: int) -> int:
    opens = sum(buf.count(t, start, end) for t in _OPEN_TAGS)
    closes = sum(buf.count(t, start, end) for t in _CLOSE_TAGS)
    return opens - closes


def _group_tag_kind(buf: bytes, i: int):
    """Given the index of b"Group" in buf, return (tag_start, +1 open / -1 close) or None."""
    for prefix, kind in ((b"<ortholog", 1), (b"<paralog", 1), (b"</ortholog", -1), (b"</paralog", -1)):
        start = i - len(prefix)
        if start >= 0 and buf.startswith(prefix, start):
            return start, kind
    return None


def _toplevel_start(buf: bytes, pos: int, depth: int, lo: int = 0) -> int:
    """Start of the top-level group enclosing `pos`, `depth` being the nesting depth at `pos`."""
    while depth > 0:
        i = buf.rfind(b"Group", lo, pos)
        if i < 0:
            raise ValueError("unbalanced group tags while looking for the rootHOG start")
        tag = _group_tag_kind(buf, i)
        pos = i
        if tag is None:
            continue
        pos, kind = tag
        depth -= kind
    return pos


def _toplevel_end(buf: bytes, pos: int, depth: int, hi: int) -> int:
    """End (exclusive) of the top-level group enclosing `pos`."""
    while depth > 0:
        i = buf.find(b"Group", pos, hi)
        if i < 0:
            raise ValueError("unbalanced group tags while looking for the rootHOG end")
        tag = _group_tag_kind(buf, i)
        pos = i + len(b"Group")
        if tag is None:
            continue
        depth += tag[1]
        if tag[1] < 0:
            pos = buf.index(b">", pos) + 1
    return pos


def _safe_cut(buf: bytes) -> int:
    """Largest prefix of buf that ends right after a tag at group depth 0.

    buf is assumed to start at group depth 0.
    """
    base = buf.rfind(b">") + 1
    depth = _depth_delta(buf, 0, base)
    if depth == 0:
        return base
    return _toplevel_start(buf, base, depth)


def _find_attr_values(buf: bytes, end: int, attr: bytes, needles: dict, tag: bytes):
    """Yield (pos, value) of `attr="value"` inside `<tag ...>` with value in needles."""
    if len(needles) <= _FIND_THRESHOLD:
        for value in needles:
            pattern = attr + b'="' + value + b'"'
            pos = buf.find(pattern, 0, end)
            while pos >= 0:
                if pos > 0 and buf[pos - 1] in _WS:
                    tag_start = buf.rfind(b"<", 0, pos)
                    if buf.startswith(tag, tag_start) and buf[tag_start + len(tag)] in _WS:
                        yield pos, value
                pos = buf.find(pattern, pos + 1, end)
    else:
        regex = re.compile(rb"<" + tag[1:] + rb'\s[^>]*?(?<=\s)' + re.escape(attr) + rb'="([^"]*)"')
        for m in regex.finditer(buf, 0, end):
            if m.group(1) in needles:
                yield m.start(1), m.group(1)


class RootHOGLookup:
    """Find the rootHOG that contains each of the queried genes.

    Genes are matched on the given <gene> attribute: the internal OrthoXML
    ``id`` by default, or ``protId`` / ``geneId`` when explicitly requested.
    After :meth:`run`, ``results`` holds one dict per found gene with the
    rootHOG id, its taxonomic level and its number of genes.

    The taxonomic level is read from the rootHOG's ``TaxRange`` property
    (OrthoXML <= 0.4, e.g. FastOMA/OMA output) or from its ``taxonId``
    attribute (OrthoXML 0.5), resolved to a name through <taxonomy> if
    possible.
    """

    def __init__(self, source: Union[str, PathLike, IO[bytes]], queries: Iterable[str],
                 id: str = "id", chunk_size: int = DEFAULT_CHUNK_SIZE):
        self.source = source
        self.queries = {q.encode("utf-8"): q for q in queries}
        self.id = id
        self.attr = id.encode("utf-8")
        self.chunk_size = chunk_size
        self.id_attr_seen = False  # whether any <gene> carries the `id` attribute
        self.target2query = {}     # internal gene id (bytes) -> query value
        self.taxon_names = {}      # taxon id -> taxon name
        self.results = []
        self.found_queries = set()
        self._found_targets = set()
        self._attr_re = re.compile(rb"<gene\s[^>]*?(?<=\s)" + re.escape(self.attr) + rb'="')

    @property
    def done(self):
        """True once all queries that exist in the header were located."""
        return len(self._found_targets) == len(self.target2query)

    # ---- header ----
    def _scan_header(self, buf: bytes, end: int):
        if not self.id_attr_seen and self._attr_re.search(buf, 0, end):
            self.id_attr_seen = True
        for pos, value in _find_attr_values(buf, end, self.attr, self.queries, b"<gene"):
            if self.id == "id":
                gene_id = value
            else:
                tag_start = buf.rfind(b"<", 0, pos)
                m = _ID_IN_TAG_RE.search(buf, tag_start, buf.index(b">", pos))
                if m is None:
                    continue
                gene_id = m.group(1)
            self.target2query[gene_id] = self.queries[value]

    def _load_taxonomy(self, xml: bytes):
        root = etree.fromstring(xml)
        for taxon in root.iter("{*}taxon", "taxon"):
            self.taxon_names[taxon.get("id")] = taxon.get("name")

    # ---- groups ----
    def _scan_groups(self, buf: bytes, end: int):
        remaining = {t: None for t in self.target2query if t not in self._found_targets}
        hits = sorted(_find_attr_values(buf, end, b"id", remaining, b"<geneRef"))
        depth, counted_to, done_until = 0, 0, -1
        for pos, _ in hits:
            if pos < done_until:
                continue  # rootHOG already handled
            depth += _depth_delta(buf, counted_to, pos)
            counted_to = pos
            start = _toplevel_start(buf, pos, depth)
            stop = _toplevel_end(buf, pos, depth, end)
            self._report(etree.fromstring(buf[start:stop]))
            done_until = stop
            if self.done:
                return

    def _report(self, rhog):
        gene_ids = {gr.get("id").encode("utf-8") for gr in rhog.iter("{*}geneRef", "geneRef")}
        hits = gene_ids.intersection(self.target2query)
        level = self._level_of(rhog)
        for gene_id in hits:
            self._found_targets.add(gene_id)
            query = self.target2query[gene_id]
            self.found_queries.add(query)
            self.results.append({
                "query": query,
                "gene_id": gene_id.decode("utf-8"),
                "roothog_id": rhog.get("id"),
                "taxon_level": level,
                "num_genes": len(gene_ids),
            })

    def _level_of(self, rhog):
        for child in rhog:
            if isinstance(child.tag, str) and etree.QName(child).localname == "property" \
                    and child.get("name") == "TaxRange":
                return child.get("value")
        taxon_id = rhog.get("taxonId")
        if taxon_id is not None:
            return self.taxon_names.get(taxon_id, taxon_id)
        return None

    # ---- driver ----
    def run(self):
        with auto_open(self.source, "rb") as fh:
            buf = b""
            in_groups = False
            taxonomy = None  # bytes of <taxonomy> while being captured
            while True:
                chunk = fh.read(self.chunk_size)
                buf += chunk
                eof = not chunk

                if not in_groups:
                    m = _GROUPS_START_RE.search(buf)
                    end = m.start() if m else (len(buf) if eof else buf.rfind(b">") + 1)

                    # capture <taxonomy> (used to name taxonId levels)
                    if taxonomy is None:
                        t0 = buf.find(b"<taxonomy", 0, end)
                        if t0 >= 0:
                            taxonomy = buf[t0:end]
                    elif taxonomy is not False:
                        taxonomy += buf[:end]
                    if taxonomy:
                        t1 = taxonomy.find(b"</taxonomy>")
                        if t1 >= 0:
                            self._load_taxonomy(taxonomy[:t1 + len(b"</taxonomy>")])
                            taxonomy = False  # done

                    self._scan_header(buf, end)
                    buf = buf[end:]
                    if m:
                        in_groups = True
                        if not self.id_attr_seen or not self.target2query:
                            return self  # nothing can match
                    elif eof:
                        return self

                if in_groups:
                    end = len(buf) if eof else _safe_cut(buf)
                    self._scan_groups(buf, end)
                    buf = buf[end:]
                    if self.done or eof:
                        return self
