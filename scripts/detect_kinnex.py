#!/usr/bin/env python3
"""Detect Kinnex-style concatemer reads by counting multiple polyA tracts.

For each sampled read the script reports the number of polyA hits and, for
each hit, a tuple with the 100 bp after the polyA hit and the first 100 bp of
the read. Output is JSON Lines (one JSON object per read). Works with BAM
input (requires pysam) or FASTQ/FASTQ.gz.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
from typing import Iterable, List, Tuple
from collections import Counter
import edlib
from infer_adapters import compute_core_and_barcodes
try:
    import pysam
except Exception:  # pragma: no cover - optional dependency
    pysam = None

MERGE_GAP = 1
def _build_weighted_consensus(
    seq_counts: List[Tuple[str, int]],
    min_support_frac: float,
    n_threshold: float,
    align: str,
) -> str:
    '''
    Function to generate a consensus sequences from counter

    Args:
        seq_counts (List[Tuple[str, int]]): Counter of adapter sequences
        min_support_frac (float): minimum support fraction for consensus
        n_threshold (float): N threshold for consensus
        align (str): alignment direction ("left" or "right")

    Returns:
        str: consensus sequence
    '''
    if not seq_counts:
        return ""
    # max length of sequences (needed for padding)
    maxlen = max(len(s) for s, _ in seq_counts)
    padded = []
    
    # Padding of the reads, right for 3' consensus and left for 5' consensus
    # ex. AAAAAGCT, AAAAAGGCT (with G duplication) -> GCT, GGCT -> NGCT, GGCT
    # this way the presence of insertions in 3' adaptors lead to positions in 5'
    # consensus with lots of Ns that get discarded
    for seq, weight in seq_counts:
        if align == "right":
            pad = "N" * (maxlen - len(seq)) + seq
        else:
            pad = seq + "N" * (maxlen - len(seq))
        padded.append((pad, weight))

    # list to store the consensus bases
    collected: List[str] = []
    for idx in range(maxlen):
        counter = Counter()
        total = 0
        n_weight = 0
        # Build the base counts for this position
        for seq, weight in padded:
            base = seq[idx]
            counter[base] += weight
            total += weight
            if base == "N":
                n_weight += weight

        # Compute statistics and decide consensus base
        non_n = {base: w for base, w in counter.items() if base != "N"}
        base, count = Counter(non_n).most_common(1)[0]
        support = total - n_weight
        # Ratio of N votes too high
        if n_weight / total > n_threshold:
            collected.append("N")
        # Ratio of top base too low
        elif count / support < min_support_frac:
            collected.append("N")
        else:
            collected.append(base)

    return "".join(collected)

def parse_args():
    p = argparse.ArgumentParser(description="Detect Kinnex concatemer signatures")
    p.add_argument("input", help="Input BAM or FASTQ (gz) file")
    p.add_argument("--sample", type=int, default=1000, help="Number of reads to sample (default: 1000)")
    p.add_argument("--min-poly", type=int, default=15, help="Minimum consecutive A bases to call a polyA hit (default: 10)")
    p.add_argument("--post-len", type=int, default=100, help="Number of bases to extract after polyA (default: 100)")
    p.add_argument("--head-len", type=int, default=100, help="Number of bases to extract from read start (default: 100)")
    p.add_argument("--output", default=None, help="Write JSONL output to file (default: stdout)")
    return p.parse_args()

def revcomp_seq(seq: str) -> str:
    """Return the reverse complement of ``seq`` (A/T/C/G/N only)."""
    table = str.maketrans("ATCGN", "TAGCN")
    return seq.translate(table)[::-1]

def read_fastq_seqs(path: str, limit: int) -> Iterable[tuple[str, str]]:
    opener = gzip.open if path.lower().endswith(".gz") else open
    with opener(path, "rt") as fh:
        count = 0
        while count < limit:
            hdr = fh.readline()
            if not hdr:
                break
            seq = fh.readline().strip()
            fh.readline()  # +
            fh.readline()  # qual
            name = hdr.strip().split()[0][1:] if hdr.startswith("@") else f"read_{count}"
            yield name, seq
            count += 1


def read_bam_seqs(path: str, limit: int) -> Iterable[tuple[str, str]]:
    if pysam is None:
        raise RuntimeError("pysam is required for BAM input; install pysam or provide FASTQ")
    with pysam.AlignmentFile(path, "rb", check_sq=False) as bam:
        count = 0
        for aln in bam.fetch(until_eof=True):
            if aln.is_secondary or aln.is_supplementary:
                continue
            seq = aln.query_sequence
            if not seq:
                continue
            yield aln.query_name, seq
            count += 1
            if count >= limit:
                break


def iterate_sequences(path: str, sample: int) -> Iterable[tuple[str, str]]:
    path_lower = path.lower()
    if path_lower.endswith(".bam"):
        return read_bam_seqs(path, sample)
    if path_lower.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
        return read_fastq_seqs(path, sample)
    # fallback: try BAM first if pysam present
    if pysam is not None and os.path.exists(path):
        try:
            return read_bam_seqs(path, sample)
        except Exception:
            return read_fastq_seqs(path, sample)
    return read_fastq_seqs(path, sample)


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
    
def process_sequences(seq_iter, args, out_fh=None):
    processed_reads: List[dict] = []
    count = 0
    for name, seq in seq_iter:
        count += 1
        poly_hits = detect_polyA(seq, args.min_poly)
        if len(poly_hits) <= 2:
            seq_rev = revcomp_seq(seq)
            polyT_hits = detect_polyA(seq_rev, args.min_poly)
            if len(polyT_hits) > len(poly_hits):
                poly_hits = polyT_hits
                seq = seq_rev
        head = seq[: args.head_len]
        hits = []
        for start, end in poly_hits:
            post = seq[end : end + args.post_len]
            hits.append({"pos": start, "seq": post})

        record = {"read": name, "num_polyA": len(poly_hits), "hits": hits, "head": head}
        processed_reads.append(record)
        line = json.dumps(record)
        if out_fh:
            out_fh.write(line + "\n")

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


def infer_kinnex_from_reads(processed_reads: List[dict]) -> None:
    '''
    Try to infer it we have kinnex data. In case it is try to identify the
    primers.

    Args:
        processed_reads (List[dict]): relevant information to compute primers:
            - num_polyA: number of polyA hits
            - head: 5' head sequence
            - hits: list of dicts with 'seq' key for 3' adapter sequences
    '''
    # Counter for the different features
    segments_counter = Counter() # of number of polyA segments per read
    heads_counter = Counter() # of 5' head sequences
    adapter_3p_counter = Counter() # of 3' adapter sequences
    
    # this gives us the length of the multiplexing array
    for rec in processed_reads:
        segments_counter[rec["num_polyA"]] +=1
    most_common_num = segments_counter.most_common(1)[0][0]

    # Recover sequences and their frequencies to build consensuses
    for rec in processed_reads:
        heads_counter[rec["head"]] +=1
        if rec["num_polyA"] == most_common_num:
            for hit in rec["hits"]:
                adapter_3p_counter[hit["seq"]] +=1
    
    # Build initial consensuses
    seq_3p, _  = compute_core_and_barcodes(adapter_3p_counter, "3p", 0.6,0.2, False)
    print(f"Consensus 3' adapter sequence: {seq_3p}")
    seq_5p, _ = compute_core_and_barcodes(heads_counter, "5p", 0.6,0.2, False)

    print(f"Consensus 5' head sequence: {seq_5p}")
    adapter_5p = build_common_ends(seq_5p, seq_3p, "3p")
    print(f"Common sequence between 5' head and 3' adapter: {adapter_5p}")
    a_sequence = seq_5p[:-len(adapter_5p)]
    print(a_sequence)
    print(len(a_sequence)+ len(adapter_5p))
    print(len(seq_5p))
    # count the number of N blocks in adapter_5p
    n_blocks_5p_adapter = count_n_blocks(adapter_5p)
    print(f"N blocks in adapter_5p: {n_blocks_5p_adapter}")
    # count the number of N blocks in seq_3p before adapter_5p
    n_blocks = count_n_blocks(seq_3p[:-len(adapter_5p)])
    print(f"N blocks before adapter_5p: {n_blocks}")
    
def main() -> None:
    args = parse_args()

    out_fh = open(args.output, "w") if args.output else None

    seq_iter = iterate_sequences(args.input, args.sample)
    processed_reads = process_sequences(seq_iter, args, out_fh)


    if out_fh:
        out_fh.close()
    infer_kinnex_from_reads(processed_reads)

if __name__ == "__main__":
    main()
