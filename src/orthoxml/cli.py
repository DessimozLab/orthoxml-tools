# orthoxml/cli.py

import os
import argparse
import sys
import json
from orthoxml import OrthoXMLTree
from orthoxml import __version__
from orthoxml.parsers import process_stream_orthoxml
from orthoxml.converters.to_nhx import orthoxml_to_newick
from orthoxml.converters.from_nhx import orthoxml_from_newicktrees
from orthoxml.custom_parsers import BasicStats, GenePerTaxonStats, PrintTaxonomy, RootHOGCounter, SplitterByRootHOGS
from orthoxml.logger import get_logger

logger = get_logger(__name__)

def load_tree(filepath, validate=False, score_id=None, score_threshold=None, filter_strategy=None):
    """Load OrthoXML tree from file without applying any completeness filter."""
    try:
        if score_id and not all([score_id, score_threshold, filter_strategy]):
            raise ValueError("If score_id is provided, score_threshold and filter_strategy must also be provided.")
        if score_id and score_threshold and filter_strategy:
            if filter_strategy == "bottomup":
                tree = OrthoXMLTree.from_file(filepath,
                                              validate=validate,
                                              score_id=score_id,
                                              score_threshold=score_threshold,
                                              high_child_as_rhogs=True,
                                              keep_low_score_parents=False)
            elif filter_strategy == "topdown":
                tree = OrthoXMLTree.from_file(filepath,
                                          validate=validate,
                                          score_id=score_id,
                                          score_threshold=score_threshold,
                                          high_child_as_rhogs=False,
                                          keep_low_score_parents=False)
            else:
                raise ValueError("Invalid filter strategy. Use 'bottomup' or 'topdown'.")
        else:
            tree = OrthoXMLTree.from_file(filepath,
                                          validate=validate)
        return tree
    except Exception as e:
        print(f"Error loading file: {e}")
        sys.exit(1)

def handle_stats(args):
    with BasicStats(args.infile) as parser:
        for _ in parser.parse():
            pass
        logger.info(f"Number of species: {parser.species_count}")
        logger.info(f"Number of genes: {parser.gene_count}")
        logger.info(f"Number of rootHOGs: {parser.rhog_count}")
        logger.info(f"Number of leave taxa: {parser.leave_taxon_count}")
        logger.info(f"Total number of taxa: {parser.all_taxa_count}")

def handle_gene_stats(args):
    with GenePerTaxonStats(args.infile) as parser:
        for _ in parser.parse():
            pass
        parser.compute_taxon_counts()
        
        if args.outfile:
            with open(args.outfile, 'w') as outfile:
                json.dump(parser.taxonomy_counts, outfile, indent=4)
            logger.info(f"Gene count per taxon written to {args.outfile}")
        else:
            logger.info(parser.taxonomy_counts)

def handle_taxonomy(args):
    with PrintTaxonomy(args.infile) as parser:
        for _ in parser.parse():
            pass
        print(parser.taxonomy.to_str())

def handle_export(args):
    tree = load_tree(args.infile)
    if args.type == "pairs":
        pairs = tree.to_ortho_pairs(filepath=args.outfile if args.outfile else None)
        for pair in pairs:
            print(pair)
    elif args.type == "groups":
        groups = tree.to_ogs(filepath=args.outfile if args.outfile else None)
        for group in groups:
            print(group)
    else:
        print("Unknown export type specified.")

def handle_split_streaming(args):
    infile_name = args.infile.split("/")[-1]

    with RootHOGCounter(args.infile) as parser:
        for _ in parser.parse():
            pass
        print(f"Count {parser.rhogs_count}")

    for rhog in range(1, parser.rhogs_count + 1):
        process_stream_orthoxml(args.infile,
                        os.path.join(args.outdir, f"{rhog}_{infile_name}"),
                        parser_cls=SplitterByRootHOGS,
                        parser_kwargs={"rhogs_number": rhog})

def handle_conversion_to_nhx(args):
    infile = args.infile
    outdir = args.outdir
    xref_tag = args.xref_tag

    trees = orthoxml_to_newick(infile, xref_tag=xref_tag)

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # write trees to files
    for treeid_hog, tree in trees.items():
        tree_file_i = os.path.join(outdir, f"tree_{treeid_hog}.nwk")
        logger.debug(f"Writing tree {treeid_hog} to {tree_file_i}")
        with open(tree_file_i,'w') as handle:
            handle.write(tree)
        handle.close()

    logger.info(f"We wrote {len(trees)} trees  in nhx format from the input HOG orthoxml {infile} in {outdir}.")
    logger.info("You can visualise each tree using https://beta.phylo.io/viewer/ as extended newick format.")

def handle_conversion_from_nhx(args):
    orthoxml_from_newicktrees(
        args.infile,
        args.outfile,
        label_to_event=None, 
        label_to_id_and_species=None
    )

def handle_filter(args):

    try:
        tree = load_tree(args.infile,
                         score_id=args.score_name,
                         score_threshold=args.threshold,
                         filter_strategy=args.strategy)

    except Exception as e:
        print(f"Error filtering tree: {e}")
        sys.exit(1)

    if args.outfile:
        try:
            tree.to_orthoxml(args.outfile)
            print(f"Filtered file written to {args.outfile}")
        except Exception as e:
            print(f"Error writing filtered tree: {e}")
            sys.exit(1)
    else:
        try:
            print(tree.to_orthoxml(args.outfile))
        except Exception as e:
            print(f"Error serializing filtered tree to string: {e}")
            sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="Command Line Interface for orthoxml-tools")

    parser.add_argument("-v", "--version", action="version",
                        version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(
        title="subcommands", dest="command", required=True)

    # Validate subcommand
    validate_parser = subparsers.add_parser("validate", help="Validate an OrthoXML file")
    validate_parser.add_argument("--infile", required=True, help="Path to the OrthoXML file")
    validate_parser.set_defaults(func=lambda args: load_tree(args.infile, validate=True))

    # Stats subcommand
    stats_parser = subparsers.add_parser("stats", help="Show statistics of the OrthoXML tree")
    stats_parser.add_argument("--infile", required=True, help="Path to the OrthoXML file")
    stats_parser.set_defaults(func=handle_stats)

    # Gene Stats
    gene_stats_parser = subparsers.add_parser("gene-stats", help="Show gene statistics of the OrthoXML tree")
    gene_stats_parser.add_argument("--infile", required=True, help="Path to the OrthoXML file")
    gene_stats_parser.add_argument(
        "--outfile",
        help="If provided, write the gene statistics to this file; otherwise, print to stdout"
    )
    gene_stats_parser.set_defaults(func=handle_gene_stats)

    # Taxonomy subcommand
    tax_parser = subparsers.add_parser("taxonomy", help="Print the taxonomy tree")
    tax_parser.add_argument("--infile", required=True, help="Path to the OrthoXML file")
    tax_parser.set_defaults(func=handle_taxonomy)

    # Conversions
    ## OrthoXML to Newick (NHX)
    converter_to_nhx_parser = subparsers.add_parser("to-nhx", help="Convert OrthoXML to Newick (NHX) format")
    converter_to_nhx_parser.add_argument("--infile", required=True, help="Path to the OrthoXML file")
    converter_to_nhx_parser.add_argument("--outdir", required=True, help="Path to the folder where the trees will be saved")
    converter_to_nhx_parser.add_argument(
        "--xref-tag",
        default="protId",
        help="the attribute of the <gene> element that should be used to get as label for the leaves labels."
    )
    converter_to_nhx_parser.set_defaults(func=handle_conversion_to_nhx)

    ## Newick (NHX) to OrthoXML
    converter_from_nhx_parser = subparsers.add_parser("from-nhx", help="Convert Newick (NHX) to OrthoXML format")
    converter_from_nhx_parser.add_argument(
        "--infile",
        nargs="+",  # Accept one or more input files
        required=True,
        help="Paths to one or more Newick (NHX) files"
    )
    converter_from_nhx_parser.add_argument("--outfile", required=True, help="Path to the output OrthoXML file")
    converter_from_nhx_parser.set_defaults(func=handle_conversion_from_nhx)

    # Export subcommand
    export_parser = subparsers.add_parser("export", help="Export orthologous pairs or groups")
    export_parser.add_argument("--infile", required=True, help="Path to the OrthoXML file")
    export_parser.add_argument("type", choices=["pairs", "groups"], help="Type of export")
    export_parser.add_argument("--outfile", help="Output file to write the export")
    export_parser.set_defaults(func=handle_export)

    # Split subcommand
    split_parser = subparsers.add_parser("split", help="Split the tree by rootHOGs")
    split_parser.add_argument("--infile", required=True, help="Path to the OrthoXML file")
    split_parser.add_argument("--outdir", required=True, help="Path to the folder where the splitted rootHOGs will be saved")
    split_parser.set_defaults(func=handle_split_streaming)

    # Filter subcommand
    filter_parser = subparsers.add_parser("filter", help="Filter the OrthoXML tree by a completeness score")
    filter_parser.add_argument("--infile", required=True, help="Path to the OrthoXML file")
    filter_parser.add_argument(
        "--score-name",
        required=True,
        help="Name of the completeness score annotation (e.g. 'CompletenessScore')"
    )
    filter_parser.add_argument(
        "--threshold",
        type=float,
        required=True,
        help="Threshold value for the completeness score"
    )
    filter_parser.add_argument(
        "--strategy",
        choices=["bottomup", "topdown"],
        default="topdown",
        help="Filtering strategy (bottomup or topdown)"
    )
    filter_parser.add_argument(
        "--outfile",
        help="If provided, write the filtered OrthoXML to this file; otherwise, print to stdout"
)
    filter_parser.set_defaults(func=handle_filter)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
