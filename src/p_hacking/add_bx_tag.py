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

try:
    from .utils import iterate_sequences
except ImportError:
    from utils import iterate_sequences
    
DEFAULT_BC1_LEN = 40
DEFAULT_BC2_LEN = 39
PROGRESS_INTERVAL = 100_000
CHECK_END_BASES = 11


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("input", help="Input BAM file")
    parser.add_argument(
        "-o",
        "--output",
        help=(
            "Output BAM path. If omitted, will create '<input>.bx.bam'."
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
        "--primers-fasta",
        help="Path to FASTA file with primer sequences",
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
    '''
    Derive output path from input if explicit output wasn't provided

    Args:
        input_path (str): path to input bam file
        explicit_output (Optional[str]): explicit output path, if provided

    Returns:
        str: derived output path
    '''
    if explicit_output:
        return explicit_output
    if input_path.lower().endswith(".bam"):
        return f"{input_path[:-4]}.bx.bam"
    return f"{input_path}.bx.bam"

def get_RG_from_header(header: dict) -> set[str]:
    '''
    Extract read group IDs from BAM header

    Args:
        header (dict): BAM file header
    Returns:
        set[str]: set of read group IDs or empty set if no RG present
    '''
    if "RG" not in header:
        return set()
    return {rg.get("ID") for rg in header["RG"] if "ID" in rg}


def get_informative_barcode_tags(aln: pysam.AlignedSegment) -> dict[str, str]:
    '''
    Extract informative barcode tags from a read alignment

    Args:
        aln (pysam.AlignedSegment): read alignment
    Returns:
        dict[str, str]: dictionary of informative barcode tags
        the informative tags are RG and bc if present
    '''
    tags = {"RG": None, "bc": None}
    rg_tag = aln.get_tag("RG") if aln.has_tag("RG") else None
    if rg_tag:
        tags["RG"] = rg_tag
    bc_tag = aln.get_tag("bc") if aln.has_tag("bc") else None
    if bc_tag:
        tags["bc"] = bc_tag
    return tags

def get_nt_content(seq: str, nucleotide: str) -> float:
    '''
    Count the frequency of a nucleotide in a sequence

    Args:
        seq (str): nucleotide sequence
        nucleotide (str): nucleotide to count

    Returns:
        float: frequency of the nucleotide in the sequence
    '''
    return seq.count(nucleotide)/len(seq)


def find_orientation(seq: str) -> str:
    '''
    Function to determine the orientation of the read. Uses A content of 3'
    vs T content of 5'

    Args:
        seq (str): read sequence

    Returns:
        str: description of the orientation: 'fwd' o 'rev'
    '''
    # get the fragment to check
    five_prime = seq[:CHECK_END_BASES]
    three_prime = seq[-CHECK_END_BASES:]
    
    # Compute the relevant frequencies, i.e. T in 5' and A in 3'
    five_t_content = get_nt_content(five_prime, "T")
    three_a_content = get_nt_content(three_prime, "A")
    if three_a_content > five_t_content:
        orientation =  "fwd"
    else:
        orientation = "rev"
    return orientation

def complete_missing_rgs(missing_rgs: set[str], header: pysam.AlignmentHeader, primers: list[str]) -> None:
    '''
    Add the missing RG found in the reads but not in the header back to the header

    Args:
        missing_rgs (set[str]): _description_
        header (pysam.AlignmentHeader): _description_
        primers (list[str]): _description_
    '''
    
    # if there is a RG in header, use it as template
    template_rg = None
    if "RG" in header:
        template_rg = header["RG"][0]
    else:
        header["RG"] = []
        template_rg = {"ID": "template", 
                       "PL":"PACBIO", 
                       "DS":"READTYPE=CCS;",
                       "SM": "sample1",
                       "LB": "lib1",
                       "BC": "N/A"}
    for rg in missing_rgs:
        indx = rg.split("/")[1]
        idx5, idx3 = indx.split("--")
        idx5, idx3 = int(idx5), int(idx3)
        primer5 = primers[idx5]
        primer3 = primers[idx3]
        new_rg = template_rg.copy()
        new_rg["ID"] = rg
        new_rg["BC"] = f"{primer5}-{primer3}"
        header["RG"].append(new_rg)

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
    
    # if priemers were supplied use them to derive bx lengths
    primers = []
    if args.primers_fasta is not None:
        primers = [seq for seq in iterate_sequences(args.primers_fasta, 1e6)]

    count = 0
    rg_ids_reads= set()
    with pysam.AlignmentFile(in_path, "rb", check_sq=False) as inh, pysam.AlignmentFile(
        out_path, "wb", template=inh
    ) as outh:
        # Get the IDs of read groups from the header
        input_header = inh.header.to_dict()
        rg_ids_header = get_RG_from_header(input_header)

        for aln in inh.fetch(until_eof=True):
            # if not priemrs provide there is not much we can do but add constant bx
            # add record RG ids in case some is missing from header
            if not primers:
                if aln.has_tag("RG"):
                    rg_ids_reads.add(aln.get_tag("RG"))
                aln.set_tag("bx", array("i", bx_payload))
                all_tags = aln.get_tags()
                # get read orientation
                orientation = find_orientation(aln.query_sequence)
                # add bc tag with dummy values but propepr orientation
                if orientation == "fwd":
                    aln.set_tags(all_tags + [("bc", [0, 12], "S")])
                else:
                    aln.set_tags(all_tags + [("bc", [12, 0], "S")])
                outh.write(aln)
            else:
                # Get tag related barcodes
                barcode_tags = get_informative_barcode_tags(aln)
                
                # BC is the most straightforward way to get the barcode index, don't
                # need to consider orientation here as the primers are stored in the
                # same order as in the read
                if barcode_tags["bc"] is not None:
                    barcodes_index = barcode_tags["bc"][1]
                    bx_len = [len(primers[barcodes_index][0]), len(primers[barcodes_index][1])]
                    aln.set_tag("bx", array("i", bx_len))
                    outh.write(aln)
                # if RG is present save be and in case that bc was missing
                # use it to get the primers
                if barcode_tags["RG"] is not None:
                    rg_ids_reads.add(barcode_tags["RG"])
                    if barcode_tags["bc"] is None:
                        # get the index from RG
                        indx = barcode_tags["RG"].split("/")[1]
                        idx5, idx3 = indx.split("--")
                        idx5, idx3 = int(idx5), int(idx3)
                        # Reverse primers if necessary
                        orientation = find_orientation(aln.query_sequence)
                        if orientation == "rev":
                            idx5, idx3 = idx3, idx5
                        primer5 = primers[idx5]
                        primer3 = primers[idx3]
                        
                        # Priemrs lengths go into bx
                        bx_len = [len(primer5), len(primer3)]
                        aln.set_tag("bx", array("i", bx_len))
                        all_tags = aln.get_tags()
                        
                        # Primers indices go into bc
                        aln.set_tags(all_tags + [("bc", [idx5, idx3], "S")])
                        outh.write(aln)
            count += 1
            if args.verbose and count % PROGRESS_INTERVAL == 0:
                print(f"Processed {count:,} reads...", file=sys.stderr)

    if args.verbose:
        print(
            f"Finished: wrote {count:,} reads to {out_path} with bx:B:i,{bx_payload[0]},{bx_payload[1]}",
            file=sys.stderr,
        )
    
    # we compare the RG from the reads and the header to see if any is missing
    missing_rgs = rg_ids_reads - rg_ids_header
    if missing_rgs:
        if args.verbose:
            print(f"Adding {len(missing_rgs)} missing RGs to header...", file=sys.stderr)
        # Add missing RGs to header and write updated bam
        complete_missing_rgs(missing_rgs, input_header, primers)
        with pysam.AlignmentFile(out_path, "rb", check_sq=False) as inh, pysam.AlignmentFile(
            f"{out_path}.tmp", "wb", header=input_header
        ) as outh:
            for aln in inh.fetch(until_eof=True):
                outh.write(aln)
        os.replace(f"{out_path}.tmp", out_path)
        if args.verbose:
            print(f"Updated header written to {out_path}", file=sys.stderr)
        


if __name__ == "__main__":
    main(parse_args())
