import argparse
import sys
from p_hacking import detect_kinnex, infer_adapters, inspect_adapters, infer_adapters_from_concatemers, add_bx_tag, add_zm_tag, build_adapters_fasta

def main():
    parser = argparse.ArgumentParser(description="PacBio Tools Suite")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # ---------------------------------------------------------
    # Command: detect-kinnex
    # ---------------------------------------------------------
    parser_kinnex = subparsers.add_parser(
        "detect-kinnex", 
        help="Detect Kinnex concatemer signatures",
        description=detect_kinnex.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # Register arguments from the module
    detect_kinnex.add_args(parser_kinnex)
    # Set the function to call when this command is chosen
    parser_kinnex.set_defaults(func=detect_kinnex.main)

    # ---------------------------------------------------------
    # Command: infer-adapters
    # ---------------------------------------------------------
    parser_infer = subparsers.add_parser(
        "infer-adapters", 
        help="Infer multiplexed Iso-Seq primer sequences",
        description=infer_adapters.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    infer_adapters.add_args(parser_infer)
    parser_infer.set_defaults(func=infer_adapters.main)

    # ---------------------------------------------------------
    # Command: inspect-adapters
    # ---------------------------------------------------------
    parser_inspect = subparsers.add_parser(
        "inspect-adapters", 
        help="Highlight primer and polyA/polyT signatures inside reads",
        description=inspect_adapters.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    inspect_adapters.add_args(parser_inspect)
    parser_inspect.set_defaults(func=inspect_adapters.main)

    # ---------------------------------------------------------
    # Command: infer-adapters-from-concatemers
    # ---------------------------------------------------------
    parser_infer_concat = subparsers.add_parser(
        "infer-adapters-from-concatemers", 
        help="Infer primer cores from concatemer reads polyA/polyT boundaries",
        description=infer_adapters_from_concatemers.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    infer_adapters_from_concatemers.add_args(parser_infer_concat)
    parser_infer_concat.set_defaults(func=infer_adapters_from_concatemers.main)

    # ---------------------------------------------------------
    # Command: add-bx-tag
    # ---------------------------------------------------------
    parser_add_bx = subparsers.add_parser(
        "add-bx-tag", 
        help="Add or replace the per-read `bx` tag with constant barcode lengths",
        description=add_bx_tag.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    add_bx_tag.add_args(parser_add_bx)
    parser_add_bx.set_defaults(func=add_bx_tag.main)

    # ---------------------------------------------------------
    # Command: add-zm-tag
    # ---------------------------------------------------------
    parser_add_zm = subparsers.add_parser(
        "add-zm-tag", 
        help="Add/replace the per-read `zm` tag in a BAM and enrich related header/tags",
        description=add_zm_tag.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    add_zm_tag.add_args(parser_add_zm)
    parser_add_zm.set_defaults(func=add_zm_tag.main)

    # ---------------------------------------------------------
    # Command: build-adapters-fasta
    # ---------------------------------------------------------
    parser_build_fasta = subparsers.add_parser(
        "build-adapters-fasta", 
        help="Extract barcode primers from PacBio FL BAM headers and write a FASTA",
        description=build_adapters_fasta.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    build_adapters_fasta.add_args(parser_build_fasta)
    parser_build_fasta.set_defaults(func=build_adapters_fasta.main)

    # Parse arguments
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    
    # Execute the selected function
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
