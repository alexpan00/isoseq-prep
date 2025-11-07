#!/usr/bin/env python3
"""Inspect reads for residual primer signatures and polyA/polyT tracts.

This utility scans reads within a BAM file and
attempts to highlight any primer sequences that remain, along with long
polyA/polyT segments. Matches are identified with the `edlib` approximate
string-matching library so small sequencing errors do not prevent detection.
"""

from __future__ import annotations

import argparse
import os
from typing import Dict, Iterable, List, Sequence, Tuple, Optional

import edlib  # type: ignore
import pysam

RESET = "\033[0m"
COLORS = {
    "primer_5p": "\033[91m",  # bright red
    "primer_5p_rc": "\033[31m",  # red (darker) for reverse complement
    "primer_3p": "\033[94m",  # bright blue
    "primer_3p_rc": "\033[34m",  # blue (darker) for reverse complement
    "polyA": "\033[93m",  # yellow
    "polyT": "\033[96m",  # cyan
}

Annotation = Tuple[int, int, str, int, str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Highlight primer and polyA/polyT signatures inside reads."
    )
    parser.add_argument(
        "--bam",
        default="discarded.bam",
        help="Concatemer BAM file to inspect (default: discarded.bam)",
    )
    parser.add_argument(
        "--primers",
        default="test/multiplexed_barcodes.fasta",
        help="FASTA file listing primer sequences (default: test/multiplexed_barcodes.fasta)",
    )
    parser.add_argument(
        "--max-reads",
        type=int,
        default=10,
        help="Process at most this many reads (default: 10)",
    )
    parser.add_argument(
        "--primer-error-frac",
        type=float,
        default=0.1,
        help="Maximum mismatch fraction allowed when aligning primers (default: 0.1)",
    )
    parser.add_argument(
        "--poly-min-length",
        type=int,
        default=20,
        help="Minimum length for polyA/polyT detection (default: 20)",
    )
    parser.add_argument(
        "--poly-max-mismatches",
        type=int,
        default=2,
        help="Maximum mismatches permitted in a polyA/polyT stretch (default: 2)",
    )
    parser.add_argument(
        "--wrap-width",
        type=int,
        default=80,
        help="Wrap highlighted sequence every N bases (default: 80)",
    )
    return parser.parse_args()


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def load_primers(fasta_path: str) -> List[Dict[str, str]]:
    primers: List[Dict[str, str]] = []
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Primer FASTA not found: {fasta_path}")

    with open(fasta_path, "r", encoding="utf-8") as fh:
        name: Optional[str] = None
        seq_chunks: List[str] = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name and seq_chunks:
                    sequence = "".join(seq_chunks).upper()
                    primers.append({"name": name, "sequence": sequence})
                name = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if name and seq_chunks:
            sequence = "".join(seq_chunks).upper()
            primers.append({"name": name, "sequence": sequence})
    if not primers:
        raise ValueError(f"No primers found in FASTA: {fasta_path}")
    augmented: List[Dict[str, str]] = []
    for primer in primers:
        name = primer["name"]
        sequence = primer["sequence"]
        augmented.append(primer)
        augmented.append({"name": f"{name}_rc", "sequence": reverse_complement(sequence)})
    return augmented


def align_primer(
    primer_seq: str,
    read_seq: str,
    max_errors: int,
) -> Tuple[int, List[Tuple[int, int]]]:
    result = edlib.align(primer_seq, read_seq, mode="HW", task="locations", k=max_errors)
    edit_distance = result.get("editDistance", -1)
    locations = result.get("locations", [])
    return edit_distance, locations


def collect_primer_matches(
    read_seq: str,
    primers: List[Dict[str, str]],
    error_fraction: float,
) -> Tuple[List[Annotation], List[str]]:
    annotations: List[Annotation] = []
    descriptions: List[str] = []

    for primer in primers:
        name = primer["name"]
        sequence = primer["sequence"]
        direction = "rc" if name.endswith("_rc") else "fwd"
        base_name = name[:-3] if direction == "rc" else name
        orientation = "3p" if base_name.endswith("_3p") else "5p"
        color_key = f"primer_{orientation}_{'rc' if direction == 'rc' else ''}".rstrip("_")

        max_errors = max(1, int(len(sequence) * error_fraction))

        edit_distance, locations = align_primer(sequence, read_seq, max_errors=max_errors)
        if edit_distance != -1:
            for start, end in locations:
                annotations.append((start, end + 1, COLORS[color_key], 3, f"{name}") )
                descriptions.append(
                    f"Primer {name} at {start}-{end + 1} (edit distance {edit_distance})"
                )

    return annotations, descriptions


def collect_poly_matches(
    read_seq: str,
    min_len: int,
    max_mismatches: int,
) -> Tuple[List[Annotation], List[str]]:
    annotations: List[Annotation] = []
    descriptions: List[str] = []

    for base, label in (("A", "polyA"), ("T", "polyT")):
        pattern = base * min_len
        result = edlib.align(pattern, read_seq, mode="HW", task="locations", k=max_mismatches)
        locations = result.get("locations", [])
        if not locations:
            continue

        merged: List[Tuple[int, int]] = []
        for start, end in sorted(locations):
            seg_start = start
            seg_end = end + 1
            # extend exact matches surrounding the seed
            while seg_start > 0 and read_seq[seg_start - 1] == base:
                seg_start -= 1
            while seg_end < len(read_seq) and read_seq[seg_end] == base:
                seg_end += 1
            if merged and seg_start <= merged[-1][1]:
                merged[-1] = (min(merged[-1][0], seg_start), max(merged[-1][1], seg_end))
            else:
                merged.append((seg_start, seg_end))

        for seg_start, seg_end in merged:
            color = COLORS[label]
            annotations.append((seg_start, seg_end, color, 1, label))
            descriptions.append(
                f"{label} at {seg_start}-{seg_end} (length {seg_end - seg_start})"
            )

    return annotations, descriptions


def build_color_map(length: int, annotations: Sequence[Annotation]) -> List[str | None]:
    color_priority: List[Tuple[str | None, int]] = [[None, 0] for _ in range(length)]  # type: ignore
    for start, end, color, priority, _ in annotations:
        clamped_start = max(0, start)
        clamped_end = min(length, end)
        for idx in range(clamped_start, clamped_end):
            if priority >= color_priority[idx][1]:
                color_priority[idx] = [color, priority]
    return [entry[0] for entry in color_priority]


def render_with_wrap(seq: str, colors: Sequence[str | None], width: int) -> str:
    lines: List[str] = []
    for start in range(0, len(seq), width):
        end = min(len(seq), start + width)
        line_parts: List[str] = []
        current_color: str | None = None
        for idx in range(start, end):
            color = colors[idx]
            if color != current_color:
                if current_color is not None:
                    line_parts.append(RESET)
                if color is not None:
                    line_parts.append(color)
                current_color = color
            line_parts.append(seq[idx])
        if current_color is not None:
            line_parts.append(RESET)
        lines.append("".join(line_parts))
    return "\n".join(lines)


def iterate_reads(bam_path: str, limit: int) -> Iterable[pysam.AlignedSegment]:
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        count = 0
        for aln in bam.fetch(until_eof=True):
            if aln.is_secondary or aln.is_supplementary:
                continue
            if not aln.query_sequence:
                continue
            yield aln
            count += 1
            if 0 < limit <= count:
                break


def main() -> None:
    args = parse_args()

    primers = load_primers(args.primers)

    for aln in iterate_reads(args.bam, args.max_reads):
        sequence = aln.query_sequence.upper()
        read_len = len(sequence)

        primer_annotations, primer_descriptions = collect_primer_matches(
            sequence, primers, error_fraction=args.primer_error_frac
        )
        poly_annotations, poly_descriptions = collect_poly_matches(
            sequence,
            min_len=args.poly_min_length,
            max_mismatches=args.poly_max_mismatches,
        )

        all_annotations = primer_annotations + poly_annotations
        colors = build_color_map(read_len, all_annotations)

        print(f"Read: {aln.query_name} | length {read_len}")
        if primer_descriptions:
            print("  Primers:")
            for desc in primer_descriptions:
                print(f"    - {desc}")
        else:
            print("  Primers: none detected")

        if poly_descriptions:
            print("  PolyA/T segments:")
            for desc in poly_descriptions:
                print(f"    - {desc}")
        else:
            print("  PolyA/T segments: none detected")

        print()
        print(render_with_wrap(sequence, colors, width=args.wrap_width))
        print(RESET)
        print("=" * 80)


if __name__ == "__main__":
    main()
