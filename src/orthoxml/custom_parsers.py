# custom_parsers.py
from __future__ import annotations, unicode_literals, division

import os
import re
import errno
from collections import defaultdict
from .parsers import StreamOrthoXMLParser
from .logger import get_logger
from .legacy.models import Taxon, ORTHO_NS
from lxml import etree

logger = get_logger(__name__)

class BasicStats(StreamOrthoXMLParser):
    def __init__(self, source):
        super().__init__(source)
        self.gene_count = 0
        self.rhog_count = 0
        self.species_count = 0
        self.leave_taxon_count = 0
        self.all_taxa_count = 0

    def process_species(self, elem):
        """Count how many species and genes we have in the orthoxml file"""

        self.species_count += 1

        gene_tag = f"{{{self._ns}}}gene"
        genes_in_this_species = elem.findall(f".//{gene_tag}")
        num_genes = len(genes_in_this_species)
        self.gene_count += num_genes

        return None
    
    def process_taxonomy(self, elem):
        """Count how many leave taxon we have in the taxonomy"""

        taxon_tag = f"{{{self._ns}}}taxon"
        all_taxa = elem.findall(f".//{taxon_tag}")
        self.all_taxa_count = len(all_taxa)
        
        count = 0
        for taxon in all_taxa:
            has_child_taxon = any(child.tag == taxon_tag for child in taxon)
            if not has_child_taxon:
                count += 1
        self.leave_taxon_count = count

        return None

    def process_scores(self, elem):
        return None

    def process_toplevel_group(self, elem):
        self.rhog_count += 1
        return None


class GenePerTaxonStats(StreamOrthoXMLParser):
    def __init__(self, source):
        super().__init__(source)
        self.gene_count_per_taxon = defaultdict(int)
        self.header_gene_count_per_species = {}
        self.gene_to_species_name = {}
        self.taxonomy_counts = {}
        self.taxonomy_tree = None

    def process_species(self, elem):
        """Count how many genes we have per species in the orthoxml file"""

        species_name = elem.get("name")

        gene_tag = f"{{{self._ns}}}gene"
        genes_in_this_species = elem.findall(f".//{gene_tag}")
        num_genes = len(genes_in_this_species)

        self.header_gene_count_per_species[species_name] = num_genes

        for gene in genes_in_this_species:
            gene_id = gene.get("id")
            self.gene_to_species_name[gene_id] = species_name

        return None
    
    def process_toplevel_group(self, elem):
        """
        Called once for each top-level <orthologGroup> or <paralogGroup>.
        Count all geneRef's per species under this group.
        """
        if self.taxonomy_tree is None:
            return None

        gene_ref_tag = f"{{{self._ns}}}geneRef"

        # find every geneRef anywhere inside this group
        for gr in elem.findall(f".//{gene_ref_tag}"):
            gid = gr.get("id")
            species = self.gene_to_species_name.get(gid)
            if not species:
                logger.warning(
                    f"GeneRef with id '{gid}' not found in species mapping. "
                    "This may indicate a mismatch in gene IDs between header and groups."
                )            
                continue

            # accumulate into the global tally
            self.gene_count_per_taxon[species] = (
                self.gene_count_per_taxon.get(species, 0) + 1
            )

        return None

    def process_taxonomy(self, elem):
            """Build an in‐memory tree of nested <taxon> elements."""
            taxon_tag = f"{{{self._ns}}}taxon"
            def build_node(tx_elem):
                return {
                    "id":      tx_elem.get("id"),
                    "name":    tx_elem.get("name"),
                    "children":[ build_node(c) 
                                for c in tx_elem 
                                if isinstance(c, etree._Element) and c.tag==taxon_tag ]
                }

            roots = [ build_node(c) for c in elem 
                    if isinstance(c, etree._Element) and c.tag==taxon_tag ]

            if len(roots)==1:
                self.taxonomy_tree = roots[0]
            else:
                self.taxonomy_tree = {"id":None, "name":"<root>", "children":roots}
            return None

    def compute_taxon_counts(self):
        """Walk the taxonomy_tree and sum up gene_count_per_taxon into every node."""
        def recurse(node):
            if not node["children"]:
                cnt = self.gene_count_per_taxon.get(node["name"], 0)
            else:
                cnt = sum(recurse(ch) for ch in node["children"])
            self.taxonomy_counts[node["name"]] = cnt
            return cnt

        if self.taxonomy_tree is None:
            logger.warning("No taxonomy tree found. Falling back to per-species gene counts.")
            self.taxonomy_counts = dict(self.header_gene_count_per_species)
            return
        recurse(self.taxonomy_tree)

class PrintTaxonomy(StreamOrthoXMLParser):
    def __init__(self, source):
        super().__init__(source)
        self.taxonomy = None

    def process_taxonomy(self, elem):
        """Build an in‐memory tree of nested <taxon> elements."""

        if elem is not None:
            taxon_el = elem.find(f"{{{ORTHO_NS}}}taxon")
            if taxon_el is not None:
                self.taxonomy = Taxon.from_xml(taxon_el)

        return None


class GetGene2IdMapping(StreamOrthoXMLParser):
    """Get the mapping between id and geneId or protId, ..."""
    def __init__(self, source, id):
        super().__init__(source)
        self.gene_id2id_mapping = {}
        self.id = id

    def process_species(self, elem):
        gene_tag = f"{{{self._ns}}}gene"
        genes_in_this_species = elem.findall(f".//{gene_tag}")

        for gene in genes_in_this_species:
            self.gene_id2id_mapping[gene.attrib.get("id")] = gene.attrib.get(self.id)

        return None


class StreamPairsParser(StreamOrthoXMLParser):
    """
    Extends StreamOrthoXMLParser with a streaming ortholog or para-pair extractor.
    """
    def __init__(self, source, ortho_para: str):
        """
        :param source: path or file-like object for the orthoXML.
        :param ortho_para: either 'orthologGroup' or 'paralogGroup'.
        """
        super().__init__(source)
        if ortho_para not in ('orthologGroup', 'paralogGroup'):
            raise ValueError("ortho_para must be 'orthologGroup' or 'paralogGroup'")
        self.ortho_para = ortho_para

    def iter_pairs(self):
        """
        Yield (r_id, s_id) for every pair of the specified type in the file,
        in a single pass using only O(tree-depth × average‐refs‐per‐group) memory.
        """
        # Each frame holds:
        #   type       = tag name ('orthologGroup' or 'paralogGroup')
        #   own_refs   = list of geneRef IDs directly under this group
        #   child_refs = list of lists of gene IDs from each finished child group
        group_stack = []

        for event, elem in self._context:
            tag = self.strip_ns(elem.tag)

            # 1) On group start: push a new frame
            if event == 'start' and tag in ('orthologGroup', 'paralogGroup'):
                group_stack.append({
                    "type":       tag,
                    "own_refs":   [],
                    "child_refs": []
                })

            # 2) On geneRef end: record its ID, then immediately clear it
            elif event == 'end' and tag == 'geneRef':
                if group_stack:
                    group_stack[-1]["own_refs"].append(elem.get("id"))
                # free the <geneRef> element from memory
                elem.clear()
                while elem.getprevious() is not None:
                    del elem.getparent()[0]

            # 3) On group end: pop the frame, compute & yield pairs, pass up refs
            elif event == 'end' and tag in ('orthologGroup', 'paralogGroup'):
                frame      = group_stack.pop()
                own_refs   = frame["own_refs"]
                child_refs = frame["child_refs"]

                # Build the full list of IDs under this group
                all_refs = own_refs.copy()
                for cr in child_refs:
                    all_refs.extend(cr)

                # If this is the group type we're extracting, yield its pairs:
                if frame["type"] == self.ortho_para:
                    # (a) own-vs-own
                    for i in range(len(own_refs)):
                        for j in range(i + 1, len(own_refs)):
                            yield (own_refs[i], own_refs[j])

                    # (b) own-vs-each-child
                    for cr in child_refs:
                        for r in own_refs:
                            for s in cr:
                                yield (r, s)

                    # (c) between-different-children
                    for i in range(len(child_refs)):
                        for j in range(i + 1, len(child_refs)):
                            for r in child_refs[i]:
                                for s in child_refs[j]:
                                    yield (r, s)

                # Pass the aggregated ID list up to the parent frame (if any)
                if group_stack:
                    group_stack[-1]["child_refs"].append(all_refs)

                # 4) Free memory for this group element
                elem.clear()
                while elem.getprevious() is not None:
                    del elem.getparent()[0]


class StreamMaxOGParser(StreamOrthoXMLParser):
    def __init__(self, source):
        super().__init__(source)
        # map from geneId (string) → species name
        self.species_map: dict[str,str] = {}

    def process_species(self, elem):
        """Called on </species>: collect gene→species mapping."""
        sp_name = elem.get("name")
        # walk down to <gene> elements
        for gene in elem.findall(".//{%s}gene" % self._ns, namespaces=self.nsmap):
            gid = gene.get("geneId") or gene.get("id")
            if gid:
                self.species_map[gid] = sp_name
        # we return None so nothing is yielded for species
        return None

    def process_toplevel_group(self, elem):
        """
        Called on each top-level <orthologGroup> or <paralogGroup> under <groups>.
        Here `elem` is the root of one OG/PG subtree.
        We compute and return the list of gene IDs to keep.
        """
        def local_strip(tag):
            return tag.split("}",1)[-1]

        def recurse(node) -> list[str]:
            # gather direct geneRef IDs
            direct_refs = [gr.get("id")
                           for gr in node
                           if local_strip(gr.tag) == "geneRef"
                          ]
            # gather all group‐children
            child_groups = [c for c in node
                            if local_strip(c.tag) in ("orthologGroup","paralogGroup")]

            # if there are child groups, process them first
            if child_groups:
                child_kept = [recurse(c) for c in child_groups]

                # Duplication event = a <paralogGroup> that has child groups
                if local_strip(node.tag) == "paralogGroup":
                    # compute species‐counts for each branch
                    counts = []
                    for genes in child_kept:
                        # only count those we know species for
                        sps = { self.species_map[g] 
                                for g in genes 
                                if g in self.species_map }
                        counts.append(len(sps))
                    # pick the branch with the max distinct species
                    idx = counts.index(max(counts))
                    return child_kept[idx]

                else:
                    # an <orthologGroup> with children: union everything + any direct refs
                    out = []
                    for genes in child_kept:
                        out.extend(genes)
                    out.extend(direct_refs)
                    return out

            else:
                # leaf group (no sub-groups)
                if local_strip(node.tag) == "paralogGroup":
                    # keep only the first geneRef in a leaf paralogGroup
                    return direct_refs[:1]
                else:
                    # orthologGroup leaf: keep them all
                    return direct_refs

        # run our bottom-up pass and return its result
        return [recurse(elem)]



class OrthoXMLSplitter(object):
    """Convert orthoxml files with several families.

    This class provides the means to extract a subset of root HOGs (i.e.
    families) into a new output orthoxml file, or to split it and create
    for each family an individual file.

    The object should be instantiated with the input orthoxml file and
    optionally a cache_dir argument where the output orthoxml files will
    be stored. This later parameter can be overwritten in the __call__
    method call that does the work.

    .. note::

       Calls to the splitter will remove the created families from the
       loaded input file, so subsequent calls that contain a family in
       common will miss them from the second call onwards.


    :Example:

      splitter = OrthoXMLSplitter("data.orthoxml", cache_dir="./splits")
      splitter()

    will create files HOGxxxxxx.orthoxml in the ./splits directory.
    
    Code from FastOMA/utils/OrthoXMLSplitter.py"""

    def __init__(self, xml_file, cache_dir=None, release_char=None):
        self.xml_file = xml_file
        if cache_dir is not None:
            self._assert_cache_dir(cache_dir)
        if release_char is None:
            self.release_char = ""
        elif re.match(r"^[A-Z]?$", release_char):
            self.release_char = release_char
        else:
            raise ValueError(
                "unexpected value for release_char: '{}'. Needs to be a single capital ascii letter".format(
                    release_char
                )
            )
        logger.info("loading xml file {}...".format(xml_file))
        parser = etree.XMLParser(remove_blank_text=True)
        self.Etree_XML = etree.parse(self.xml_file, parser=parser)
        self.Etree_root = self.Etree_XML.getroot()
        logger.info("building lookup table for genes")
        self.gene_lookup = {gene.get("id"): gene for gene in self._iter_gene_elements()}
        logger.info("init of OrthoXMLSplitter finished")

    def _assert_cache_dir(self, cache_dir):
        # Ensure existance of cache directory (py2 compat)
        try:
            os.makedirs(cache_dir, exist_ok=True)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(cache_dir):
                pass
            else:
                raise
        self.cache_dir = cache_dir

    def _iter_gene_elements(self):
        """This method is a faster version of xpath '//ns:gene'.

        It iterates the element in sequential order"""
        for node in self.Etree_root:
            if node.tag == "{http://orthoXML.org/2011/}species":
                for gene in node.iter("{http://orthoXML.org/2011/}gene"):
                    yield gene

    def _iter_toplevel_groups(self):
        """This method yields all the root hogs sequentially."""
        for node in self.Etree_root:
            if node.tag == "{http://orthoXML.org/2011/}groups":
                for root_hog in node:
                    if root_hog.tag in {"{http://orthoXML.org/2011/}orthologGroup",
                                         "{http://orthoXML.org/2011/}paralogGroup"}:
                        yield root_hog

    def __call__(
        self,
        hogs_to_extract=None,
        single_hog_files=False,
        basename=None,
        cache_dir=None,
    ):
        """Split/extract hogs from orthoxml file based on root hogs ids.

        Split the input orthoxml or extract a subset of root hogs. If no
        argument is passed, one orthoxml file per root hog is created,
        named as 'HOGxxxxxx.orthoxml', where xxxxxx is the numeric id of
        each hog.

        The set of root hogs to be extracted can be limited by specifying
        a subset of hog ids in the hogs_to_extract parameter. If
        single_hog_files is set to true, each of these hogs will be converted
        into a single orthoxml file named as explained above. If single_hog_files
        is set to false, the whole subset of hogs will be stored in one
        orthoxml file named as specified in `basename`.

        The file(s) will be stored in the cache_dir folder which can be
        specified in the constructor or overwritten as an argument in
        this method.

        :param hogs_to_extract: list or set that contains the set of root
            hogs to be extracted. If set to None, all hogs are extracted.
        :param bool single_hog_files: whether or not to build one orthoxml
            file for all the selected hogs or individual ones.
        :param str basename: name of the output file if a subset of hogs
            is extracted into a single file.
        :param str cache_dir: folder where to store the output files.
        """
        if cache_dir is not None:
            self._assert_cache_dir(cache_dir)
        elif self.cache_dir is None:
            raise RuntimeError("cache dir to output files to is not set")

        if single_hog_files:
            if hogs_to_extract is None:
                raise RuntimeError(
                    "useless to extract all hogs into single output file"
                )
            if basename is None or not isinstance(basename, (str, bytes)):
                raise ValueError("basename needs to be specified: {}".format(basename))
            ogs = [
                og
                for og in self._iter_toplevel_groups()
                if og.get("id") in hogs_to_extract
            ]
            fn = os.path.join(self.cache_dir, basename)
            logger.debug("extracting {:d} hogs into {:s}".format(len(ogs), fn))
            self.create_new_orthoxml(fn, ogs)
        else:
            logger.info("extracting roothogs into individual files in: {:s}".format(self.cache_dir))
            for counter, og in enumerate(self._iter_toplevel_groups()):
                if hogs_to_extract is None or og.get("id") in hogs_to_extract:
                    hog_nr = og.get("id")
                    if hog_nr:
                        hog_id = hog_nr+".orthoxml" #"HOG{:07d}.orthoxml".format(hog_nr)
                    else:
                        hog_id = "HOG{:07d}.orthoxml".format(counter)
                    fname = os.path.join(self.cache_dir, hog_id)
                    logger.debug("extracting {} into {}".format(hog_id, fname))
                    self.create_new_orthoxml(fname, [og])

    def iter_generefs_in_og(self, og_node):
        for node in og_node.iterdescendants("{http://orthoXML.org/2011/}geneRef"):
            yield node

    def get_gene_via_generef(self, genesref_ids):
        # Deduplicate while preserving first-seen order from geneRef traversal.
        seen = set()
        ordered_gene_ids = []
        for gene_id in genesref_ids:
            if gene_id not in seen:
                seen.add(gene_id)
                ordered_gene_ids.append(gene_id)
        return [self.gene_lookup[gene_id] for gene_id in ordered_gene_ids]

    def create_new_orthoxml(self, fn, OGs):
        """create a new orthoxml file for the passed orthologGroup elements.

        :param fn: the filename of the output file. The path needs to exists
            prior to calling this method.
        :param OGs: the orthologGroup elements that should be included in the
            new output file."""
        # Get element to store
        for og_node in OGs:
            gene_ids = [
                gene_ref_elem.get("id")
                for gene_ref_elem in self.iter_generefs_in_og(og_node)
            ]
        gene_els = self.get_gene_via_generef(gene_ids)

        # Get all information to store
        zoo = {}  # <- {key:sp_etree || value: {key:db_el || values:[list_genes]}}
        for gene_el in gene_els:  # <- for all gene el
            db_el = gene_el.getparent().getparent()
            sp_el = db_el.getparent()
            if sp_el in zoo.keys():  # <- if species already visited
                if db_el in zoo[sp_el].keys():  # <- if db already visited so add gene
                    zoo[sp_el][db_el].append(gene_el)
                else:  # <- if db not visited so add db,genes
                    zoo[sp_el][db_el] = []
                    zoo[sp_el][db_el].append(gene_el)
            else:  # <- if species not visited so add sp,db,gene
                zoo[sp_el] = {}
                zoo[sp_el][db_el] = []
                zoo[sp_el][db_el].append(gene_el)

        etree_2_dump = etree.Element("orthoXML", nsmap=self.Etree_root.nsmap)
        for attr, value in self.Etree_root.items():
            etree_2_dump.set(attr, value)

        for species_el in zoo.keys():
            species_xml = etree.Element("species")
            for attr, value in species_el.items():
                species_xml.set(attr, value)
            etree_2_dump.append(species_xml)

            for db_el in zoo[species_el].keys():
                # Add <database> into <species>
                database_xml = etree.SubElement(species_xml, "database")
                for attr, value in db_el.items():
                    database_xml.set(attr, value)

                # Add <genes> TAG into <database>
                genes_xml = etree.SubElement(database_xml, "genes")

                # Fill <genes> with <gene>
                for gene_el in zoo[species_el][db_el]:
                    gene_xml = etree.SubElement(genes_xml, "gene")
                    for attr, value in gene_el.attrib.items():
                        gene_xml.set(attr, value)

        groupsxml = etree.SubElement(etree_2_dump, "groups")
        for og_et in OGs:
            # if not og_et.get("id").startswith("HOG:{:s}".format(self.release_char)):
            #     og_et.set(
            #         "id",
            #         og_et.get("id"),
            #         #"HOG:{:s}{:07d}".format(self.release_char, int(og_et.get("id"))),
            #     )
            groupsxml.append(og_et)

        tree = etree.ElementTree(etree_2_dump)
        tree.write(
            fn, xml_declaration=True, encoding="utf-8", method="xml", pretty_print=True
        )
