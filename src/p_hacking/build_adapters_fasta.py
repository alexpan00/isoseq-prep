"""Compose a primer FASTA by harvesting barcode metadata from FL BAM headers.

Each input BAM is expected to contain PacBio-style read-group (RG) entries where
barcode information can be recovered from the RG fields:

- `DS` includes `BarcodeCount=<N>` indicating the total number of primer entries.
- `BC` (or alternatively the last colon-delimited segment of `PM`) stores the
  5' and 3' barcode sequences separated by `-`.
- `ID` follows the pattern `<hash>/<five_idx>--<three_idx>` providing the
  zero-based indices each barcode occupies in the original primer FASTA.

The script verifies that all supplied BAMs agree on the barcode count and that
no index is assigned conflicting sequences. Missing indices are emitted as a
placeholder 25nt `N` sequence to preserve primer ordering.
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, Optional, Sequence, Tuple

import pysam

PLACEHOLDER_SEQ = "N" * 25
DEFAULT_OUTPUT = "primers.fasta"


class PrimerExtractionError(RuntimeError):
    """Raised when the header lacks the information required to build primers."""


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "inputs",
        nargs="+",
        help="One or more FL BAM files containing barcode metadata in RG headers.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=DEFAULT_OUTPUT,
        help=f"Output FASTA path (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Overwrite the output file if it already exists.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Log progress information to stderr.",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract barcode primers from PacBio FL BAM headers and write a FASTA."
    )
    add_args(parser)
    return parser.parse_args()


def extract_barcode_count(ds_value: str) -> Optional[int]:
    '''
    Extract abrcode count from bam header

    Args:
        ds_value (str): DS field from RG entry

    Raises:
        PrimerExtractionError: raised if BarcodeCount value is invalid

    Returns:
        Optional[int]: Number of barcodes, or None if not found
    '''
    for part in ds_value.split(";"):
        part = part.strip()
        if part.startswith("BarcodeCount="):
            try:
                return int(part.split("=", 1)[1])
            except ValueError:
                raise PrimerExtractionError(f"Invalid BarcodeCount value in DS field: '{part}'")
    return None


def extract_sequences(rg_entry: Dict[str, str]) -> Tuple[str, str]:
    bc_value = rg_entry.get("BC")
    if not bc_value or "-" not in bc_value:
        raise PrimerExtractionError(
            "RG entry lacks barcode sequence information (expected BC tag with '5-3' sequence)."
        )
    left, right = bc_value.split("-", 1)
    return left.upper(), right.upper()


def extract_indices(rg_entry: Dict[str, str]) -> Tuple[int, int]:
    rg_id = rg_entry.get("ID")
    if not rg_id:
        raise PrimerExtractionError("RG entry missing ID field; cannot determine primer indices.")
    indx = rg_id.split("/")[1]
    try:
        idx5, idx3 = indx.split("--")
    except ValueError:
        raise PrimerExtractionError(
            f"RG ID '{rg_id}' does not match expected '<hash>/<five>--<three>' format."
        )
    
    return int(idx5), int(idx3)


def harvest_from_bam(path: str, verbose: bool = False) -> Tuple[int, Dict[int, str], Dict[int, str]]:
    '''
    Extract the number of barcodes and the sequences and index in the fasta
    file from the bam header

    Args:
        path (str): path to bam file
        verbose (bool, optional): whether to print verbose output. Defaults to False.

    Raises:
        PrimerExtractionError: if the BAM header is malformed or missing required fields
        PrimerExtractionError: if the RG entry is missing the BC tag
        PrimerExtractionError: if the RG entry is missing the PM tag
        PrimerExtractionError: if the RG entry is missing the ID tag
        PrimerExtractionError: if the RG entry is missing the DS tag
        PrimerExtractionError: if the RG entry is missing the 5' or 3' primer sequences

    Returns:
        Tuple[int, Dict[int, str], Dict[int, str]]: number of primers, 
            5' primer dict, 3' primer dict. The dicts map index to sequence.
    '''
    if verbose:
        print(f"Scanning header of {path}...", file=sys.stderr)
    with pysam.AlignmentFile(path, "rb", check_sq=False) as bam:
        header_obj = bam.header
        header_dict = header_obj.to_dict() if hasattr(header_obj, "to_dict") else dict(header_obj)
        rg_entries: Sequence[Dict[str, str]] = header_dict.get("RG", [])
        if not rg_entries:
            raise PrimerExtractionError(f"No @RG entries found in header of {path}")

        barcode_count: Optional[int] = None
        five: Dict[int, str] = {}
        three: Dict[int, str] = {}

        for rg in rg_entries:
            ds_value = rg.get("DS", "")
            bc_count = extract_barcode_count(ds_value)
            if bc_count is None:
                raise PrimerExtractionError(
                    f"RG ID '{rg.get('ID', '<unknown>')}' missing BarcodeCount in DS field."
                )
            # check that all the RG entries agree on the barcode count
            if barcode_count is None:
                barcode_count = bc_count
            elif bc_count != barcode_count:
                raise PrimerExtractionError(
                    f"Conflicting BarcodeCount values within header of {path}: {barcode_count} vs {bc_count}."
                )

            idx5, idx3 = extract_indices(rg)
            seq5, seq3 = extract_sequences(rg)

            if idx5 in five and five[idx5] != seq5:
                raise PrimerExtractionError(
                    f"5' primer index {idx5} has conflicting sequences ('{five[idx5]}' vs '{seq5}') in {path}."
                )
            if idx3 in three and three[idx3] != seq3:
                raise PrimerExtractionError(
                    f"3' primer index {idx3} has conflicting sequences ('{three[idx3]}' vs '{seq3}') in {path}."
                )

            five[idx5] = seq5
            three[idx3] = seq3

        if barcode_count is None:
            raise PrimerExtractionError(f"Unable to determine BarcodeCount for {path}")

    return barcode_count, five, three


def write_fasta(output_path: str, barcode_count: int, five: Dict[int, str], three: Dict[int, str]) -> None:
    pad = 2
    with open(output_path, "w", encoding="ascii") as fh:
        for idx in range(barcode_count):
            idx_str = f"{idx:0{pad}d}"
            if idx in five:
                end = 5
                seq = five[idx]
            elif idx in three:
                end = 3
                seq = three[idx]
            else:
                end = "unknown"
                seq = PLACEHOLDER_SEQ

            fh.write(f">{end}p_bc10{idx_str}\n{seq}\n")


def main(args: argparse.Namespace) -> None:
    if os.path.exists(args.output):
        if not args.force:
            print(
                f"Error: output file '{args.output}' already exists. Use --force to overwrite.",
                file=sys.stderr,
            )
            sys.exit(1)
        os.remove(args.output)

    expected_count: Optional[int] = None
    combined_five: Dict[int, str] = {}
    combined_three: Dict[int, str] = {}

    for path in args.inputs:
        bc_count, five, three = harvest_from_bam(path, verbose=args.verbose)
        
        # check if all the samples have the same number of barcodes
        if expected_count is None:
            expected_count = bc_count
        elif bc_count != expected_count:
            print(
                f"Error: BarcodeCount mismatch ({bc_count} vs {expected_count}) detected in '{path}'.",
                file=sys.stderr,
            )
            sys.exit(1)

        # check that all the barcodes from different samples agree
        for idx, seq in five.items():
            if idx in combined_five and combined_five[idx] != seq:
                print(
                    f"Error: 5' primer index {idx} already assigned to a different sequence.",
                    file=sys.stderr,
                )
                sys.exit(1)
            combined_five[idx] = seq

        for idx, seq in three.items():
            if idx in combined_three and combined_three[idx] != seq:
                print(
                    f"Error: 3' primer index {idx} already assigned to a different sequence.",
                    file=sys.stderr,
                )
                sys.exit(1)
            combined_three[idx] = seq

    if expected_count is None:
        print("Error: no barcode information extracted from inputs.", file=sys.stderr)
        sys.exit(1)

    write_fasta(args.output, expected_count, combined_five, combined_three)

    if args.verbose:
        print(
            f"Finished: wrote primers for {expected_count} barcode positions to {args.output}.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main(parse_args())
