#!/usr/bin/env python3
"""Detect Kinnex-style concatemer reads by counting multiple polyA tracts.

For each sampled read the script detects the number of polyA hits and, for
each hit, a tuple with the 100 bp after the polyA hit and the first 100 bp of
the read. That information is then used to infer if the data contains Kinnex
concatemers and, in that case, to extract the Kinnex linkers sequences.
"""

from __future__ import annotations

import argparse
import re
from typing import List, Tuple
from collections import Counter
import edlib
from string import ascii_uppercase
try:
    from .utils import compute_core_and_barcodes, revcomp_seq, iterate_sequences
except ImportError:
    from utils import compute_core_and_barcodes, revcomp_seq, iterate_sequences
try:
    import pysam
except Exception:  # pragma: no cover - optional dependency
    pysam = None

MERGE_GAP = 1
MIN_POLYA_LEN = 15
MIN_SUPPORT = 0.6
N_THRESHOLD = 0.2


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("input", help="Input BAM or FASTQ (gz) file")
    parser.add_argument("--sample", type=int, default=1000, help="Number of reads to sample (default: 1000)")
    parser.add_argument("--min-poly", type=int, default=MIN_POLYA_LEN, help="Minimum consecutive A bases to call a polyA hit (default: 15)")
    parser.add_argument("--post-len", type=int, default=100, help="Number of bases to extract after polyA (default: 100)")
    parser.add_argument("--head-len", type=int, default=100, help="Number of bases to extract from read start (default: 100)")
    parser.add_argument("--output", default=None, help="Write JSONL output to file (default: stdout)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Detect Kinnex concatemer signatures")
    add_args(p)
    return p.parse_args()


def detect_polyA(seq: str, min_poly: int) -> list[tuple[int, int]]:
    """Return list of (start, end) positions for polyA runs in the sequence."""
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
    seq = seq.upper()
    pattern = rf"A{{{min_poly},}}"
    polyA_canditates = [(m.start(), m.end()) for m in re.finditer(pattern, seq)]
    if not polyA_canditates:
        return []
    polyA_merged = merge_regions(polyA_canditates)
    return polyA_merged
    
def process_sequences(seq_iter, min_poly: int, head_len: int, post_len: int) -> List[dict]:
    """Process sequences to detect polyA hits and extract relevant info.
    
    Args:
        seq_iter: Iterable of (name, sequence) tuples
        min_poly: Minimum length of polyA to detect
        head_len: Length of 5' head sequence to extract
        post_len: Length of sequence to extract after polyA hit
    Returns:
        List[dict]: List of processed read information with polyA hits
    """
    processed_reads: List[dict] = []
    count = 0
    for name, seq in seq_iter:
        count += 1
        poly_hits = detect_polyA(seq, min_poly)
        
        # Check reverse complement for more polyA hits if few found
        if len(poly_hits) <= 2:
            seq_rev = revcomp_seq(seq)
            polyT_hits = detect_polyA(seq_rev, min_poly)
            if len(polyT_hits) > len(poly_hits):
                poly_hits = polyT_hits
                seq = seq_rev
        head = seq[: head_len]
        hits = []
        for start, end in poly_hits:
            post = seq[end : end + post_len]
            hits.append({"pos": start, "seq": post})

        record = {"read": name, "num_polyA": len(poly_hits), "hits": hits, "head": head}
        processed_reads.append(record)

    return processed_reads

def build_common_ends(seq1: str, seq2: str, end: str) -> str:
    """Build common sequence from the two ends"""
    if end == "3p":
        seq1 = seq1[::-1]
        seq2 = seq2[::-1]
    i = 0
    while i < min(len(seq1), len(seq2)) and seq1[i] == seq2[i]:
        i += 1
    common_seq = seq1[:i]
    if end == "3p":
        common_seq = common_seq[::-1]
    return common_seq
    
def count_n_blocks(seq: str) -> list[tuple[int, int]]:
    """Count the number of N blocks in a sequence and their starts and ends"""
    n_blocks = []
    pattern = r"N+"
    for m in re.finditer(pattern, seq):
        n_blocks.append((m.start(), m.end()))
    return n_blocks


def infer_kinnex_from_reads(processed_reads: List[dict], verbose: bool) -> list[str]:
    '''
    Try to infer it we have kinnex data. In case it is try to identify the
    primers.

    Args:
        processed_reads (List[dict]): relevant information to compute primers:
            - num_polyA: number of polyA hits
            - head: 5' head sequence
            - hits: list of dicts with 'seq' key for 3' adapter sequences
        verbose (bool): whether to print verbose output
    Returns:
        List[str]: list of kinnex linkers sequences including the first adapter
    '''
    # Counter for the different features
    segments_counter = Counter() # of number of polyA segments per read
    heads_counter = Counter() # of 5' head sequences
    adapter_3p_counter = Counter() # of 3' adapter sequences
    
    # this gives us the length of the multiplexing array
    for rec in processed_reads:
        segments_counter[rec["num_polyA"]] +=1
    most_common_num = segments_counter.most_common(1)[0][0]
    
    # check if we have concatemers
    if most_common_num < 2:
        print("No concatemers detected. The most common number of polyA segments is less than 2.")
        return []

    # Recover sequences and their frequencies to build consensuses
    for rec in processed_reads:
        # only from full concatemers
        if rec["num_polyA"] == most_common_num:
            heads_counter[rec["head"]] +=1
            for hit in rec["hits"]:
                adapter_3p_counter[hit["seq"]] +=1
    
    # Build initial consensuses
    seq_3p, _  = compute_core_and_barcodes(adapter_3p_counter, "3p", MIN_SUPPORT, N_THRESHOLD, False)
    seq_5p, _ = compute_core_and_barcodes(heads_counter, "5p", MIN_SUPPORT, N_THRESHOLD, False)
    if verbose:
        print(f"Consensus 5' head sequence: {seq_5p}")
        print(f"Consensus 3' adapter sequence: {seq_3p}")
    '''
    Build common sequence between 5' head and 3' sequence
    This works because the 3' sequence contains at the end the 5'
    adapter sequence:
    5' adapter: Kinnex - 5' Adapter - Insert
    3' adapter: Insert - polyA - 3' Adapter - Kinnex - 5' Adapter
    '''
    adapter_5p = build_common_ends(seq_5p, seq_3p, "3p")
    if verbose:
        print(f"Common sequence between 5' head and 3' adapter: {adapter_5p}")
    
    # Extract the first kinnex adapter sequence
    a_sequence = seq_5p[:-len(adapter_5p)]
    a_len = len(a_sequence)
        
    # Extract the 3' adapter
    adapter_3p = seq_3p[:-len(adapter_5p)]

    # Extract kinnex linkers and add the a_sequence at the start
    linkers_seqs = extract_kinnex_linkers(processed_reads, most_common_num, adapter_3p)
    linkers_seqs = [a_sequence] + linkers_seqs

    return linkers_seqs
    
def extract_kinnex_linkers(processed_reads: List[dict], most_common_num: int, adapter_3p: str) -> List[str]:
    '''
    Given the 3' adapter, extract the kinnex linkers that are immediately
    downstream. In order to find the adapter map it using edlib.

    Args:
        processed_reads (List[dict]): reads with polyA hits and sequences
        most_common_num (int): most common number of polyA segments
        adapter_3p (str): 3' adapter consensus sequence

    Returns:
        List[str]: list of extracted kinnex linkers
    '''
    def extract_adapter_bait(adapter_3p: str) -> str:
        """Extract adapter bait sequence from 3' adapter consensus
        
        Args:
            adapter_3p (str): 3' adapter consensus sequence

        Returns:
            str: Extracted adapter bait sequence
        """
        if len(n_bloccks_3p) == 1: # just kinnex adapater, no multiplexing
            n_block_start = n_bloccks_3p[0][0]
            adapter_bait = adapter_3p[(n_block_start-10):n_block_start]
        else:
            # take sequence before last N block
            n_block_start = n_bloccks_3p[-1][0]
            adapater_bait_start = max(n_bloccks_3p[-2][1], n_block_start - 10)
            adapter_bait = adapter_3p[adapater_bait_start:n_block_start]
        return adapter_bait
    
    # Get the possition and length of kinnex linkers based on N blocks
    n_bloccks_3p = count_n_blocks(adapter_3p)
    n_len = n_bloccks_3p[-1][1] - n_bloccks_3p[-1][0] if n_bloccks_3p else 0
    
    # Extract smaller sequence that will be used for alignment
    adapter_bait = extract_adapter_bait(adapter_3p)

    # we leverage the fact that kinnex linkers are in the same order in
    # all of the reads
    l_linker_counter: list[Counter[str]] = [Counter() for _ in range(most_common_num)]
    for rec in processed_reads:
        if rec["num_polyA"] == most_common_num:
            for idx, hit in enumerate(rec["hits"]):
                seq = hit["seq"] # extract the sequences downstream of polyA
                # Find the adapter position using edlib
                aligment = edlib.align(adapter_bait, seq, mode="HW", task="locations")
                aligment_end = aligment["locations"][0][1] if aligment["locations"] else None
                if aligment_end is not None:
                    # Aligment end is inclusive, so we add 1
                    aligment_end +=1
                    seq_kinnex = seq[aligment_end:aligment_end+n_len]
                    l_linker_counter[idx][seq_kinnex] +=1
    
    # if errors are low enough the most common sequence should be the
    # correct kinnex linker
    l_kinnex_linkers = [count.most_common(1)[0][0] for count in l_linker_counter]
    
    return l_kinnex_linkers


def write_linkers_fasta(linkers: List[str], outfile: str)-> None:
    """Write the detected kinnex linkers to a FASTA file.

    Args:
        linkers (List[str]): List of kinnex linker sequences
        outfile (str): Output FASTA file path
    """
    with open(outfile, "w") as fh:
        for idx, seq in enumerate(linkers):
            fh.write(f">{ascii_uppercase[idx]}\n{seq}\n")

def main(args: argparse.Namespace) -> None:
    outfile = args.output if args.output else "kinnex_linkers.fasta"
    # Opens the file and builds an iterator over sequences
    seq_iter = iterate_sequences(args.input, args.sample, return_name=True)

    # preprocess sequences to detect polyA hits
    processed_reads = process_sequences(seq_iter, args.min_poly, args.head_len, args.post_len)

    # extract the linkers
    linkers = infer_kinnex_from_reads(processed_reads, args.verbose)
    
    if linkers:
        write_linkers_fasta(linkers, outfile)
    else:
        print("No kinnex linkers detected.")

if __name__ == "__main__":
    args = parse_args()
    main(args)
