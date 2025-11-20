"""Shared utilities for scripts in this directory.

Currently contains consensus-building helpers extracted from multiple scripts.
"""
from __future__ import annotations

from collections import Counter
import os
from typing import List, Tuple, Iterable

import pysam


def _build_weighted_consensus(
    seq_counts: List[Tuple[str, int]],
    min_support_frac: float,
    n_threshold: float,
    align: str,
) -> str:
    """Build a weighted consensus sequence from a list of (sequence, weight).

    Args:
        seq_counts: list of (sequence, weight) tuples (usually Counter.most_common())
        min_support_frac: minimum fraction of non-N votes required for a column
        n_threshold: maximum fraction of N votes allowed per column
        align: "left" or "right" alignment for padding

    Returns:
        consensus sequence (may contain 'N')
    """
    if not seq_counts:
        return ""
    # max length of sequences (needed for padding)
    maxlen = max(len(s) for s, _ in seq_counts)
    padded = []

    # Padding of the reads, right for 3' consensus and left for 5' consensus
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
        if not non_n:
            collected.append("N")
            continue
        base, count = Counter(non_n).most_common(1)[0]
        support = total - n_weight
        # Ratio of N votes too high
        if total == 0:
            collected.append("N")
        elif n_weight / total > n_threshold:
            collected.append("N")
        # Ratio of top base too low
        elif support == 0 or count / support < min_support_frac:
            collected.append("N")
        else:
            collected.append(base)

    return "".join(collected)


def build_consensus(
    seq_counts: List[Tuple[str, int]],
    min_support_frac: float,
    n_threshold: float,
    align: str,
) -> str:
    """Wrapper around _build_weighted_consensus for convenience."""
    return _build_weighted_consensus(seq_counts, min_support_frac, n_threshold, align)


def revcomp_seq(seq: str) -> str:
    """Return the reverse complement of ``seq`` (A/T/C/G/N)."""
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


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


def read_bam_seqs(path: str, limit: int, return_name: bool) -> Iterable[str]:
    if pysam is None:
        raise RuntimeError("pysam is required for BAM input; install pysam or provide FASTQ")
    with pysam.AlignmentFile(path, "rb", check_sq=False) as bam:
        count = 0
        for aln in bam.fetch(until_eof=True):
            seq = aln.query_sequence
            if seq:
                if return_name:
                    yield aln.query_name, seq
                else:
                    yield seq
                count += 1
                if count >= limit:
                    break


def iterate_sequences(path: str, sample: int, return_name: bool = False) -> Iterable[str]:
    path_lower = path.lower()
    if path_lower.endswith(".bam"):
        return read_bam_seqs(path, sample, return_name=return_name)
    if path_lower.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
        return read_fastq_seqs(path, sample)
    # fallback attempts
    if pysam is not None and os.path.exists(path):
        try:
            return read_bam_seqs(path, sample)
        except Exception:
            return read_fastq_seqs(path, sample)
    return read_fastq_seqs(path, sample)