#!/usr/bin/env python3
"""Infer primer cores from concatemer readspolyA/polyT boundaries.

The script samples FL reads, looks for terminal polyA/polyT signatures that
indicate residual adapters, collects the sequences between those signals, and
builds a weighted consensus to recover the shared primer core. An optional
FASTA export captures the detected core for downstream refinement steps.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
from collections import Counter
from typing import Iterable, List, Optional, Tuple
from utils import build_consensus

try:
    import pysam
except Exception:  # pragma: no cover - optional dependency
    pysam = None


MERGE_GAP = 1
N_THRESHOLD = 0.1
DEFAULT_SUPPORT_FRAC = 0.8
END_WINDOW = 50  # bases from read end to search for polyA/polyT


def revcomp_seq(seq: str) -> str:
    """Return the reverse complement of ``seq`` (A/T/C/G/N only)."""
    table = str.maketrans("ATCGN", "TAGCN")
    return seq.translate(table)[::-1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Infer multiplexed primers from polyA/polyT signatures"
    )
    parser.add_argument("input", help="Input BAM or FASTQ file (optionally gzipped)")
    parser.add_argument(
        "--sample",
        type=int,
        default=2000,
        help="Number of reads to sample (default: 2000)",
    )
    parser.add_argument(
        "--min-poly",
        type=int,
        default=10,
        help="Minimum consecutive A/T bases to call polyA/polyT (default: 10)",
    )
    parser.add_argument(
        "--search-window",
        type=int,
        default=END_WINDOW,
        help="Window length from read end/start to search for polyA/polyT (default: 50)",
    )
    parser.add_argument(
        "--min-support",
        type=float,
        default=DEFAULT_SUPPORT_FRAC,
    help="Minimum support fraction for consensus columns (default: 0.3)",
    )
    parser.add_argument(
        "--n-threshold",
        type=float,
        default=N_THRESHOLD,
        help="Maximum fraction of 'N' votes allowed per consensus column (default: 0.8)",
    )
    parser.add_argument(
        "--no-fasta",
        action="store_true",
        help="Disable writing adapter FASTA sequences.",
    )
    parser.add_argument(
        "--fasta-prefix",
        default=None,
        help=(
            "Base filename for FASTA output. Defaults to <input_basename>_primers.fasta "
            "if omitted."
        ),
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print progress every 1000 reads",
    )
    return parser.parse_args()


def read_fastq_seqs(path: str, limit: int) -> Iterable[str]:
    opener = gzip.open if path.lower().endswith(".gz") else open
    with opener(path, "rt") as fh:
        for i in range(limit):
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            if not seq:
                break
            yield seq.strip()
            fh.readline()  # '+'
            fh.readline()  # quality


def read_bam_seqs(path: str, limit: int) -> Iterable[str]:
    if pysam is None:
        raise RuntimeError("pysam is required for BAM input; install pysam or provide FASTQ")
    with pysam.AlignmentFile(path, "rb", check_sq=False) as bam:
        count = 0
        for aln in bam.fetch(until_eof=True):
            seq = aln.query_sequence
            if seq:
                yield seq
                count += 1
                if count >= limit:
                    break


def iterate_sequences(path: str, sample: int) -> Iterable[str]:
    path_lower = path.lower()
    if path_lower.endswith(".bam"):
        return read_bam_seqs(path, sample)
    if path_lower.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
        return read_fastq_seqs(path, sample)
    # fallback attempts
    if pysam is not None and os.path.exists(path):
        try:
            return read_bam_seqs(path, sample)
        except Exception:
            return read_fastq_seqs(path, sample)
    return read_fastq_seqs(path, sample)

'''
Potential improvements for this functions:
    - Right now the merge is unnecessary because if there is no regex hit
        there is no polyA to merge, and if there is a hit the trailing polyA
        is leading polyT is already selected.
    - the logic of finding the polyA/T and extracting the adapter could be
        split into two functions for clarity.
'''
def find_polyA_candidates(
    seq: str, min_poly: int, search_window: int
) -> Tuple[Optional[str], Optional[str]]:
    """Return the adapter sequence (if any) and hit type ('A' or 'T')."""
    seq = seq.upper()
    
    def merge_regions(regions: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
        merged: List[Tuple[int, int]] = []
        start, end = regions[0]
        for s, e in regions[1:]:
            if s - end <= MERGE_GAP:
                end = max(end, e)
            else:
                merged.append((start, end))
                start, end = s, e
        merged.append((start, end))
        return merged

    tail = seq[-search_window:] if len(seq) >= search_window else seq
    matches = list(re.finditer(r"A{%d,}" % min_poly, tail))

    if matches:
        # Define read start position based on tail offset
        base_offset = len(seq) - len(tail)
        regions: List[Tuple[int, int]] = []
        for m in matches:
            regions.append((base_offset + m.start(), base_offset + m.end()))
        merged: List[Tuple[int, int]] = merge_regions(regions)
        
        # use the last merged region to extract adapter. This avoids picking up
        # internal polyA stretches.
        start, end = merged[-1]
        if end < len(seq):
            # the adapter is everything downstream of the polyA stretch
            adapter = seq[end :]
            return adapter, "A"
        elif end == len(seq):
            # no adapter found, polyA is at the end of the read
            return "", "A"

    head = seq[:search_window] if len(seq) >= search_window else seq
    matches_t = list(re.finditer(r"T{%d,}" % min_poly, head))
    if matches_t:
        regions_t: List[Tuple[int, int]] = [(m.start(), m.end()) for m in matches_t]
        merged_t: List[Tuple[int, int]] = merge_regions(regions_t)
        start, end = merged_t[0]
        if start > 0:
            # the adapter is everything upstream of the polyT stretch
            adapter = seq[:start]
            return revcomp_seq(adapter), "T"
        elif start == 0:
            # no adapter found, polyT is at the start of the read
            return "", "T"
    return None, None


def resolve_output_path(
    input_path: str,
    prefix: Optional[str],
    *,
    default_suffix: str,
    default_extension: str,
    allowed_extensions: Tuple[str, ...],
) -> str:
    """Build an output path using an optional prefix and sensible default stem."""

    if prefix:
        base = prefix
    else:
        stem = os.path.splitext(os.path.basename(input_path))[0] or "primer"
        base = f"{stem}_{default_suffix}"

    if base.lower().endswith(allowed_extensions):
        return base
    return f"{base}{default_extension}"

def concatemer_primer(seq: str) -> str:
    """Return a putative primer sequence if polyA/polyT signals flank it."""
    adapter: str = ""
    
    # strip the ends of the sequence to search for polyA/T
    region_5p = seq[:END_WINDOW].upper()
    region_3p = seq[-END_WINDOW:].upper()

    # search for polyA and polyT candidates independently in both ends
    adapter_5p, hit_5p = find_polyA_candidates(region_5p, min_poly=10, search_window=len(region_5p))
    adapter_3p, hit_3p = find_polyA_candidates(region_3p, min_poly=10, search_window=len(region_3p))

    # Concatemers usually have a polyA-inset-polyT structure, so both ends
    # should have hits.
    if (hit_5p is not None and hit_3p is not None) and (hit_5p != hit_3p):
        adapter = adapter_5p if len(adapter_5p) > len(adapter_3p) else adapter_3p
    
    return adapter


def find_primer_in_concatemers(seq_iter: Iterable[str]) -> Counter[str]:
    '''
    Try to find 3' primers in concatemer sequences

    Args:
        seq_iter (Iterable[str]): iterable of reads

    Returns:
        Counter[str]: putative primer counter
    '''
    adapter_counts: Counter[str] = Counter()
    for seq in seq_iter:
        adapter = concatemer_primer(seq)
        if adapter:
            adapter_counts[adapter] += 1
    return adapter_counts

def write_adapter_fasta(results: dict, output_path: str) -> None:
    """Write adapter FASTA sequences up to the detected knee."""

    if not output_path:
        raise ValueError("output_path must be provided for FASTA export")

    with open(output_path, "w") as fh:
        for label in ("5p", "3p"):
            data = results.get(label, {})
            core = data.get("core", "")
            header = f">{label}_adapter_core"
            fh.write(header + "\n")
            fh.write(core + "\n")
        
        
def main() -> None:
    args = parse_args()

    # Generate sequence iterator
    seq_iter = iterate_sequences(args.input, args.sample)

    # find putative primers in concatemer reads
    adapter_counts = find_primer_in_concatemers(seq_iter).most_common()

    # build consensus from the adapters found
    core = build_consensus(adapter_counts, 
                           args.min_support, 
                           args.n_threshold, 
                           "left").strip("N")
    # prepare results
    primers_res = {
        "5p": {
            "core": core[:-1],
        },
        "3p": {
            "core": core,
        },
    }

    print(f"Detected primer core length: {len(core)}")
    if core:
        print(f"Detected primer core: {core}")
        if not args.no_fasta:
            fasta_path = resolve_output_path(
                args.input,
                args.fasta_prefix,
                default_suffix="primers",
                default_extension=".fasta",
                allowed_extensions=(".fa", ".fasta", ".fna"),
            )
            print(f"Writing adapter sequences to: {fasta_path}")
            write_adapter_fasta(primers_res, fasta_path)
    else:
        print("No primer core detected.")

if __name__ == "__main__":
    main()
