# orthoxml/cli.py

import argparse
import sys
from orthoxml import OrthoXMLTree

def load_tree(filepath, validate, completeness):
    """Load OrthoXML tree from file."""
    kwargs = {'validate': validate}
    if completeness is not None:
        kwargs['CompletenessScore_threshold'] = completeness
    try:
        tree = OrthoXMLTree.from_file(filepath, **kwargs)
        return tree
    except Exception as e:
        print(f"Error loading file: {e}")
        sys.exit(1)

def handle_stats(args):
    tree = load_tree(args.file, args.validate, args.completeness)
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
    tree = load_tree(args.file, args.validate, args.completeness)
    print("Taxonomy Tree:")
    print(tree.taxonomy.to_str())

def handle_export(args):
    tree = load_tree(args.file, args.validate, args.completeness)
    if args.type == "pairs":
        pairs = tree.to_ortho_pairs(filepath=args.outfile if args.outfile else None)
        print("Orthologous Pairs:")
        for pair in pairs:
            print(pair)
    elif args.type == "groups":
        groups = tree.to_ogs(filepath=args.outfile if args.outfile else None)
        print("Orthologous Groups:")
        for key, group in groups.items():
            print(f"{key}: {group}")
    else:
        print("Unknown export type specified.")

def handle_split(args):
    tree = load_tree(args.file, args.validate, args.completeness)
    trees = tree.split_by_rootHOGs()
    print(f"Split into {len(trees)} trees based on rootHOGs.")
    for idx, t in enumerate(trees):
        print(f"\nTree {idx + 1}:")
        print(t.groups)

def main():
    parser = argparse.ArgumentParser(
        description="Command Line Interface for orthoxml-tools")
    parser.add_argument("--validate", action="store_true", help="Validate the OrthoXML file")
    parser.add_argument("--completeness", type=float, help="Completeness score threshold for filtering")
    
    # Global argument for the file path
    parser.add_argument("file", help="Path to the OrthoXML file")
    
    subparsers = parser.add_subparsers(title="subcommands", dest="command", required=True)
    
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
    
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
