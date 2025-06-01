# orthoxml/cli.py

import argparse
import sys
from orthoxml import OrthoXMLTree
from orthoxml import __version__

def load_tree(filepath, validate, score_id=None, score_threshold=None, filter_strategy=None):
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
    tree = load_tree(args.file, args.validate)
    base_stats = tree.base_stats()
    gene_stats = tree.gene_stats()
    print("Base Stats:")
    for key, value in base_stats.items():
        print(f"  {key}: {value}")
    print("\nGene Stats:")
    for taxon_id, count in gene_stats.items():
        print(f"  Taxon {taxon_id}: {count} genes")
    if args.outfile:
        with open(args.outfile, 'w') as f:
            f.write("Metric,Value\n")
            for key, value in base_stats.items():
                f.write(f"{key},{value}\n")
        print(f"\nStats written to {args.outfile}")

def handle_taxonomy(args):
    tree = load_tree(args.file, args.validate)
    print("Taxonomy Tree:")
    print(tree.taxonomy.to_str())

def handle_export(args):
    tree = load_tree(args.file, args.validate)
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

def handle_split(args):
    tree = load_tree(args.file, args.validate)
    trees = tree.split_by_rootHOGs()
    print(f"Split into {len(trees)} trees based on rootHOGs.")
    for idx, t in enumerate(trees):
        print(f"\nTree {idx + 1}:")
        print(t.groups)

def handle_filter(args):

    try:
        tree = load_tree(args.file, args.validate,
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
    parser.add_argument("--validate", action="store_true",
                        help="Validate the OrthoXML file")

    # Global argument for the file path
    parser.add_argument("file", help="Path to the OrthoXML file")

    subparsers = parser.add_subparsers(
        title="subcommands", dest="command", required=True)

    # Stats subcommand
    stats_parser = subparsers.add_parser("stats", help="Show statistics of the OrthoXML tree")
    stats_parser.add_argument("--outfile", help="Output file to write stats")
    stats_parser.set_defaults(func=handle_stats)

    # Taxonomy subcommand
    tax_parser = subparsers.add_parser("taxonomy", help="Print the taxonomy tree")
    tax_parser.set_defaults(func=handle_taxonomy)

    # Export subcommand
    export_parser = subparsers.add_parser("export", help="Export orthologous pairs or groups")
    export_parser.add_argument("type", choices=["pairs", "groups"], help="Type of export")
    export_parser.add_argument("--outfile", help="Output file to write the export")
    export_parser.set_defaults(func=handle_export)

    # Split subcommand
    split_parser = subparsers.add_parser("split", help="Split the tree by rootHOGs")
    split_parser.set_defaults(func=handle_split)

    # Filter subcommand
    filter_parser = subparsers.add_parser("filter", help="Filter the OrthoXML tree by a completeness score")
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
