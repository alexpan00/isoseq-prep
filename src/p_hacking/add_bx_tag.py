"""Add or replace the per-read `bx` tag with constant barcode lengths.

Assigns a two-element integer array stored in the
SAM/BAM tag `bx`. The default payload is `bx:B:i,40,39`, corresponding to
barcode lengths of 40 and 39 bases. You can override those lengths via CLI
options. It also adds a `bc:B:S` tag with a dummy sequence value of `[0,12]`.
Probably need to understand better how isoseq refine used those tags to provide
more meaningful values.
"""

import argparse
import os
import sys
from typing import Optional
from array import array

import pysam

DEFAULT_BC1_LEN = 40
DEFAULT_BC2_LEN = 39
PROGRESS_INTERVAL = 100_000


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("input", help="Input BAM file (indexed or unindexed)")
    parser.add_argument(
        "-o",
        "--output",
        help=(
            "Output BAM path. If omitted, will create '<input>.bx.bam' when the input"
            " ends with '.bam', otherwise '<input>.bx.bam'."
        ),
    )
    parser.add_argument(
        "--bc1-len",
        type=int,
        default=DEFAULT_BC1_LEN,
        help=f"Length to store as the first bx value (default: {DEFAULT_BC1_LEN})",
    )
    parser.add_argument(
        "--bc2-len",
        type=int,
        default=DEFAULT_BC2_LEN,
        help=f"Length to store as the second bx value (default: {DEFAULT_BC2_LEN})",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Overwrite output if it exists",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help=f"Print progress every {PROGRESS_INTERVAL:,} reads",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Add/replace `bx` array tag in a BAM file."
    )
    add_args(parser)
    return parser.parse_args()


def derive_output_path(input_path: str, explicit_output: Optional[str]) -> str:
    if explicit_output:
        return explicit_output
    if input_path.lower().endswith(".bam"):
        return f"{input_path[:-4]}.bx.bam"
    return f"{input_path}.bx.bam"


def main(args: argparse.Namespace) -> None:
    in_path = args.input
    out_path = derive_output_path(in_path, args.output)

    if os.path.exists(out_path):
        if not args.force:
            print(
                f"Error: output file '{out_path}' already exists. Use --force to overwrite.",
                file=sys.stderr,
            )
            sys.exit(1)
        os.remove(out_path)

    bx_payload = [args.bc1_len, args.bc2_len]

    count = 0
    with pysam.AlignmentFile(in_path, "rb", check_sq=False) as inh, pysam.AlignmentFile(
        out_path, "wb", template=inh
    ) as outh:
        for aln in inh.fetch(until_eof=True):
            aln.set_tag("bx", array("i", bx_payload))
            all_tags = aln.get_tags()
            aln.set_tags(all_tags + [("bc", [0, 12], "S")])
            outh.write(aln)
            count += 1
            if args.verbose and count % PROGRESS_INTERVAL == 0:
                print(f"Processed {count:,} reads...", file=sys.stderr)

    if args.verbose:
        print(
            f"Finished: wrote {count:,} reads to {out_path} with bx:B:i,{bx_payload[0]},{bx_payload[1]}",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main(parse_args())
