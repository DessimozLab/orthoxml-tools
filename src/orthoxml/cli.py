# orthoxml/cli.py

import os
import io
import argparse
import json
from lxml import etree
from orthoxml import __version__
from orthoxml.parsers import process_stream_orthoxml
from orthoxml.converters.to_nhx import orthoxml_to_newick
from orthoxml.converters.from_nhx import orthoxml_from_newicktrees
from orthoxml.converters.from_orthofinder import convert_csv_to_orthoxml
from orthoxml.custom_parsers import (
    BasicStats,
    GenePerTaxonStats,
    PrintTaxonomy,
    RootHOGCounter,
    IndexNthRootHOG,
    OutputNthRootHOG,
    StreamPairsParser,
    GetGene2IdMapping,
    StreamMaxOGParser,
)
from orthoxml.streamfilters import filter_hogs, FilterStrategy, enum_to_str
from orthoxml.logger import get_logger, set_logger_level
from orthoxml.utils import validate_xml

logger = get_logger(__name__)

def handle_validation(args):
    try:
        # Get version from root attribute using fast/partial parse
        for event, elem in etree.iterparse(args.infile, events=('start',), remove_comments=True):
            orthoxml_version = elem.attrib.get('version')
            break  # Only need root element
    except (etree.XMLSyntaxError, OSError) as e:
        raise Exception(f"Failed to parse OrthoXML file '{args.infile}': {e}")

    if not orthoxml_version:
        raise Exception("Missing OrthoXML version attribute in the root element.")

    if not validate_xml(args.infile, orthoxml_version):
        raise Exception(
            f"OrthoXML file '{args.infile}' is not valid for version {orthoxml_version}"
        )

    logger.info(
        f"OrthoXML file '{args.infile}' is valid for version {orthoxml_version}"
    )

def handle_stats(args):
    with BasicStats(args.infile) as parser:
        for _ in parser.parse():
            pass
        print(f"Number of species: {parser.species_count}")
        print(f"Number of genes: {parser.gene_count}")
        print(f"Number of rootHOGs: {parser.rhog_count}")
        print(f"Number of leave taxa: {parser.leave_taxon_count}")
        print(f"Total number of taxa: {parser.all_taxa_count}")

def handle_gene_stats(args):
    with GenePerTaxonStats(args.infile) as parser:
        for _ in parser.parse():
            pass
        parser.compute_taxon_counts()
        
        if args.outfile:
            with open(args.outfile, 'w') as outfile:
                json.dump(parser.taxonomy_counts, outfile, indent=4)
            print(f"Gene count per taxon written to {args.outfile}")
        else:
            print(parser.taxonomy_counts)

def handle_taxonomy(args):
    with PrintTaxonomy(args.infile) as parser:
        for _ in parser.parse():
            pass
        print(parser.taxonomy.to_str())

def handle_export_pairs(args):
    if args.type == "ortho":
        ortho_para = "orthologGroup"
    elif args.type == "para":
        ortho_para = "paralogGroup"
    else:
        print("Unknown export type specified.")

    chunk_size  = args.chunk_size
    buffer_size = args.buffer_size

    if args.id != "id":
        with GetGene2IdMapping(args.infile, args.id) as parser:
            for _ in parser.parse():
                pass
            mapping = parser.gene_id2id_mapping

    with StreamPairsParser(args.infile, ortho_para) as parser, \
        open(args.outfile, 'wb') as raw, \
        io.BufferedWriter(raw, buffer_size=buffer_size) as buf:  # 4 MiB

        write = buf.write
        lines = []
        for count, (r_id, s_id) in enumerate(parser.iter_pairs(), 1):
            # prepare bytes once
            if args.id != "id":
                print(r_id, s_id)
                r_id, s_id = mapping.get(r_id, r_id), mapping.get(s_id, s_id)
                print(r_id, s_id)
                print()
            lines.append(f"{r_id}\t{s_id}\n".encode('utf8'))
            if count % chunk_size == 0:
                write(b''.join(lines))
                lines.clear()
        if lines:
            write(b''.join(lines))

def handle_export_ogs(args):
    if args.id != "id":
        with GetGene2IdMapping(args.infile, args.id) as parser:
            for _ in parser.parse():
                pass
            mapping = parser.gene_id2id_mapping
    else:
        mapping = None

    with StreamMaxOGParser(args.infile) as parser, \
         open(args.outfile, 'w', encoding='utf-8') as out:
        out.write("Group\tProtein\n")
        c = 1
        for tag, kept_gene_list in parser.parse():
            if tag == "orthologGroup":
                group_id = f"OG_{c:07d}"
                for gene in kept_gene_list:
                    protein_id = mapping[gene] if mapping else gene
                    out.write(f"{group_id}\t{protein_id}\n")
                c += 1

def handle_split_streaming(args):
    infile_name = args.infile.split("/")[-1]

    with RootHOGCounter(args.infile) as counter:
        counter.parse_through()
        logger.info(f"Processing {counter.rhogs_count} root-level groups...")

    for rhog in range(1, counter.rhogs_count + 1):

        with IndexNthRootHOG(args.infile, rhog) as index:
            index.parse_through()
            logger.debug(f"Group {rhog} has {len(index.present_genes)} gene refs")

            process_stream_orthoxml(args.infile,
                                    os.path.join(args.outdir, f"{rhog}_{infile_name}"),
                                    parser_cls=OutputNthRootHOG,
                                    parser_kwargs={
                                        "rhogs_number": rhog,
                                        "present_genes": index.present_genes})

def handle_conversion_to_nhx(args):
    infile = args.infile
    outdir = args.outdir
    xref_tag = args.xref_tag
    encode_levels_as_nhx = args.encode_levels

    trees = orthoxml_to_newick(infile, xref_tag=xref_tag, encode_levels_as_nhx=encode_levels_as_nhx)

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

def handle_conversion_from_orthofinder(args):
    convert_csv_to_orthoxml(
        csv_path=args.infile,
        xml_path=args.outfile,
        xmlns="http://orthoXML.org/2011/",
        root_attrib={"version":"0.4","origin":"orthoXML.org","originVersion":"1"}
    )

def handle_filter(args):
    filter_hogs(args.infile, args.outfile, args.threshold,
                strategy=args.strategy)


def main():
    # Parser for shared options between commands
    shared_args_parser = argparse.ArgumentParser(add_help=False)
    shared_args_parser.add_argument(
        "--log",
        default="WARNING",
        help="Set the logging level [DEBUG | INFO | WARNING | ERROR | CRITICAL]",
    )

    parser = argparse.ArgumentParser(
        description="Command Line Interface for orthoxml-tools",
    )
    parser.add_argument("-v", "--version", action="version",
                        version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(
        title="subcommands", dest="command", required=True)

    # Validate subcommand
    validate_parser = subparsers.add_parser("validate",
                                            parents=[shared_args_parser],
                                            help="Validate an OrthoXML file")
    validate_parser.add_argument("--infile", required=True, help="Path to the OrthoXML file")
    validate_parser.set_defaults(func=handle_validation)

    # Stats subcommand
    stats_parser = subparsers.add_parser("stats",
                                         parents=[shared_args_parser],
                                         help="Show statistics of the OrthoXML tree")
    stats_parser.add_argument("--infile", required=True, help="Path to the OrthoXML file")
    stats_parser.set_defaults(func=handle_stats)

    # Gene Stats
    gene_stats_parser = subparsers.add_parser("gene-stats",
                                              parents=[shared_args_parser],
                                              help="Show gene statistics of the OrthoXML tree")
    gene_stats_parser.add_argument("--infile", required=True, help="Path to the OrthoXML file")
    gene_stats_parser.add_argument(
        "--outfile",
        help="If provided, write the gene statistics to this file; otherwise, print to stdout"
    )
    gene_stats_parser.set_defaults(func=handle_gene_stats)

    # Taxonomy subcommand
    tax_parser = subparsers.add_parser("taxonomy",
                                       parents=[shared_args_parser],
                                       help="Print the taxonomy tree")
    tax_parser.add_argument("--infile", required=True, help="Path to the OrthoXML file")
    tax_parser.set_defaults(func=handle_taxonomy)

    # Conversions
    ## OrthoXML to Newick (NHX)
    converter_to_nhx_parser = subparsers.add_parser("to-nhx",
                                                    parents=[shared_args_parser],
                                                    help="Convert OrthoXML to Newick (NHX) format")
    converter_to_nhx_parser.add_argument("--infile", required=True, help="Path to the OrthoXML file")
    converter_to_nhx_parser.add_argument("--outdir", required=True, help="Path to the folder where the trees will be saved")
    converter_to_nhx_parser.add_argument(
        "--xref-tag",
        default="protId",
        help="the attribute of the <gene> element that should be used to get as label for the leaves labels."
    )
    converter_to_nhx_parser.add_argument(
        "--encode-levels",
        action="store_true",
        help="If set, encode group levels as NHX comments in the output tree."
    )
    converter_to_nhx_parser.set_defaults(func=handle_conversion_to_nhx)

    ## Newick (NHX) to OrthoXML
    converter_from_nhx_parser = subparsers.add_parser("from-nhx",
                                                      parents=[shared_args_parser],
                                                      help="Convert Newick (NHX) to OrthoXML format")
    converter_from_nhx_parser.add_argument(
        "--infile",
        nargs="+",  # Accept one or more input files
        required=True,
        help="Paths to one or more Newick (NHX) files"
    )
    converter_from_nhx_parser.add_argument("--outfile", required=True, help="Path to the output OrthoXML file")
    converter_from_nhx_parser.set_defaults(func=handle_conversion_from_nhx)

    ## OrthoGroup CSV to OrthoXML
    converter_from_ortho_parser = subparsers.add_parser("from-csv",
                                                        parents=[shared_args_parser],
                                                        help="Convert OrthoGroup CSV to OrthoXML format")
    converter_from_ortho_parser.add_argument("--infile", required=True, help="Paths to OrthoGroup CSV file")
    converter_from_ortho_parser.add_argument("--outfile", required=True, help="Path to the output OrthoXML file")
    converter_from_ortho_parser.set_defaults(func=handle_conversion_from_orthofinder)

    # Export pairs subcommand
    export_parser = subparsers.add_parser("export-pairs",
                                          parents=[shared_args_parser],
                                          help="Export orthologous pairs")
    export_parser.add_argument("--infile", required=True, help="Path to the OrthoXML file")
    export_parser.add_argument("type", choices=["ortho", "para"], help="Type of export")
    export_parser.add_argument("--outfile", required=True, help="Output file to write the export")
    export_parser.add_argument(
        "--id", default="id",
        help="the identifier used in output, default to id. other values: geneId, protId"
    )
    export_parser.add_argument(
        "--chunk-size", type=int, default=20_000,
        help="Number of pairs to buffer before each write (default: 20,000)"
    )
    export_parser.add_argument(
        "--buffer-size", type=int, default=(4 << 20),
        help="Internal buffer size (in bytes) for writing (default: 4 MiB)"
    )
    export_parser.set_defaults(func=handle_export_pairs)

    # Export OGs subcommand
    export_og_parser = subparsers.add_parser("export-ogs",
                                             parents=[shared_args_parser],
                                             help="Export orthologous groups")
    export_og_parser.add_argument("--infile", required=True, help="Path to the OrthoXML file")
    export_og_parser.add_argument("--outfile", required=True, help="Output file to write the export")
    export_og_parser.add_argument(
        "--id", default="id",
        help="the identifier used in output, default to id. other values: geneId, protId"
    )
    export_og_parser.set_defaults(func=handle_export_ogs)

    # Split subcommand
    split_parser = subparsers.add_parser("split",
                                         parents=[shared_args_parser],
                                         help="Split the tree by rootHOGs")
    split_parser.add_argument("--infile", required=True, help="Path to the OrthoXML file")
    split_parser.add_argument("--outdir", required=True, help="Path to the folder where the splitted rootHOGs will be saved")
    split_parser.set_defaults(func=handle_split_streaming)

    # Filter subcommand
    filter_parser = subparsers.add_parser("filter",
                                          parents=[shared_args_parser],
                                          help="Filter the OrthoXML tree by CompletenessScore.")
    filter_parser.add_argument("--infile", required=True, help="Path to the OrthoXML file")
    filter_parser.add_argument(
        "--threshold",
        type=float,
        required=True,
        help="Threshold value for the completeness score"
    )
    filter_parser.add_argument(
        "--strategy",
        choices=[enum_to_str(e) for e in FilterStrategy],
        default=enum_to_str(FilterStrategy.CASCADE_REMOVE),
        help="Filtering strategy (cascade-remove, extract)"
    )
    filter_parser.add_argument(
        "--outfile",
        required=True,
        help="Write the filtered OrthoXML to this file"
    )
    filter_parser.set_defaults(func=handle_filter)

    args = parser.parse_args()
    set_logger_level(args.log)
    args.func(args)

if __name__ == "__main__":
    main()
