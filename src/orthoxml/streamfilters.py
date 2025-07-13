# streamfilters.py

import abc
import logging
from typing import Iterable

from .parsers import StreamOrthoXMLParser, process_stream_orthoxml

logger = logging.getLogger(__name__)


class HOGFilter(metaclass=abc.ABCMeta):
    def remove(self, node):
        pass


class ExistingScoreBasedHOGFilter(HOGFilter):
    def __init__(self, score:str, value:float):
        self.score = score
        self.value = value

    def remove(self, node):
        for score in node.iterfind(f'./{{http://orthoXML.org/2011/}}score[@id="{self.score}"]'):
            if float(score.get('value')) < self.value:
                return True
        return False


class FilterAllSubHOGTopDown(StreamOrthoXMLParser):
    """Top-Down filtering streamer class

    This Stream Filter class trims sub-orthologGroups in a top-down fashion.
    Whenever a (sub) orthologGroup element should be removed, the entire
    sub-tree will be removed."""
    def __init__(self, source, filters:Iterable[HOGFilter]=None):
        super().__init__(source)
        self.filters = list(filters)

    def process_toplevel_group(self, elem):
        # check if rootnode needs to be removed
        if any(filt.remove(elem) for filt in self.filters):
            # root element need to be removed, don't return anything
            return None
        to_rem = []
        for hog in elem.iterfind(f'.//{{{self._ns}}}orthologGroup'):
            if any(filt.remove(hog) for filt in self.filters):
                to_rem.append(hog)

        logger.info(f"will remove {len(to_rem)} hogs")
        pos_childs = tuple(f"{{{self._ns}}}{z}" for z in ('orthologGroup', 'paralogGroup', 'geneRef'))
        for h in to_rem:
            parent = h.getparent()
            parent.remove(h)
            if sum(c.tag in pos_childs for c in parent) == 0:
                logger.info("removing also empty hog {}".format(parent))
                if parent == elem:
                    return None
                to_rem.append(parent)
        return elem


class FilterHOGBottomUp(StreamOrthoXMLParser):
    """Bottom-up filtering streamer class

    This Stream Filter class trims sub-orthologGroups in a bottom-up fashion.
    Whenever a (sub) orthologGroup element should be removed, that sub-tree is
    added as a new toplevel orthologGroup if it contains at least min_hog_size
    genes.
    """
    def __init__(self, source, filters: Iterable[HOGFilter] = None, min_hog_size=2):
        super().__init__(source)
        self.filters = list(filters)
        self.min_hog_size = min_hog_size

    def process_toplevel_group(self, elem):

        def should_remove(node):
            if any(filt.remove(node) for filt in self.filters):
                return True
            return False

        def add_subhog_to_be_yield(node):
            while True:
                children = [child for child in node.iterchildren()
                            if self.strip_ns(child.tag) in ("orthologGroup", "paralogGroup", "geneRef")]
                if len(children) != 1:
                    break
                node = children[0]
            if len(children) > 1 and sum(1 for c in node.iterfind(f'{{{self._ns}}}geneRef')) >= self.min_hog_size:
                new_roothogs.append(node)

        def _filter_ortholog_group(cur):
            # TODO: This function does not yet properly work! paralogGroups will not be processed!
            # Recursively process child orthologGroups
            children = list(cur.findall(f'./{{{self._ns}}}orthologGroup'))
            for child in children:
                result = _filter_ortholog_group(child)
                if result is None:
                    cur.remove(child)

            # Now decide whether to keep the current element
            if should_remove(cur):
                return None
            return elem

        logger.critical("THIS FUNCTION DOES NOT YET PROPERLY WORK!")
        new_roothogs = []
        res = _filter_ortholog_group(elem)
        if res is not None:
            new_roothogs.insert(0, res)
        return new_roothogs


def filter_hogs(source_orthoxml, out, filter: HOGFilter, strategy: str = "top-down"):
    """Filter hogs according to the given strategy with a given HOGFilter"""
    if strategy == "bottom-up":
        parser_cls = FilterHOGBottomUp
    elif strategy == "top-down":
        parser_cls = FilterAllSubHOGTopDown
    else:
        raise ValueError(f"unknown strategy {strategy}")
    process_stream_orthoxml(source_orthoxml,
                            out,
                            parser_cls=parser_cls,
                            parser_kwargs={"filters": [filter]})

