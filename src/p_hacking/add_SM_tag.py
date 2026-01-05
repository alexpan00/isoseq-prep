#!/usr/bin/env python3
"""Replace SM tags in BAM header by the given sample names.

This script reads a BAM file and replaces the sample names in the SM tags
of the read group (RG) headers with new sample names provided by the user.
It outputs a new BAM file with the updated headers. This is useful because
the SM tag is used by isoseq collapse to provide per sample quantification.
"""

from __future__ import annotations

import argparse
import os
from typing import Dict, Iterable, List, Sequence, Tuple, Optional

import pysam



def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--bam",
        required=True,
        help="Input BAM file to modify",
    )
    parser.add_argument(
        "--sample-names",
        required=True,
        help="Comma-separated list of new sample names to replace in SM tags",
    )
    parser.add_argument(
        "--output-bam",
        default=None,
        help="Output BAM file with updated SM tags. Default:input bam base name with _SM.bam suffix",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Replace SM tags in BAM header with new sample names."
    )
    add_args(parser)
    return parser.parse_args()




def main(args: argparse.Namespace) -> None:
    # Prepare output BAM filename
    if args.output_bam is None:
        base, ext = os.path.splitext(args.bam)
        args.output_bam = f"{base}_SM{ext}"

    # Parse new sample names
    new_sample_names = args.sample_names.split(",")
    
    # Open input BAM file
    with pysam.AlignmentFile(args.bam, "rb", check_sq=False) as in_bam:
        header = in_bam.header.to_dict()
        
        # Update SM tags in RG headers
        rg_headers = header.get("RG", [])
        if len(new_sample_names) != len(rg_headers):
            raise ValueError("Number of new sample names must match number of RG headers in BAM.")
        
        for rg, new_sm in zip(rg_headers, new_sample_names):
            rg["SM"] = new_sm
        
        # Create output BAM file with updated header
        with pysam.AlignmentFile(args.output_bam, "wb", header=header) as out_bam:
            for read in in_bam:
                out_bam.write(read)
    
    print(f"Updated BAM written to {args.output_bam}")


if __name__ == "__main__":
    main(parse_args())
