"""Shared utilities for scripts in this directory.

Currently contains consensus-building helpers extracted from multiple scripts.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import os
import gzip
from typing import List, Tuple, Iterable

import pysam
import edlib

################
#     MISC     #
################

def revcomp_seq(seq: str) -> str:
    """Return the reverse complement of ``seq`` (A/T/C/G/N)."""
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]

    
#################################
#     Primer reconstruction     #
#################################

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

def compute_core_and_barcodes(
    counter: Counter[str],
    orientation: str,
    min_support_frac: float,
    n_threshold: float,
    multiplexed: bool,
) -> Tuple[str, Counter[str]]:
    '''
    Function to call the core primer sequence and potentially
    identify the sample barcodes

    Args:
        counter (Counter[str]): counter of adapter sequences
        orientation (str): 5' or 3' end
        min_support_frac (float): minimum support fraction for consensus
        n_threshold (float): N threshold for consensus
        multiplexed (bool): whether the data is multiplexed

    Returns:
        Tuple[str, Counter[str]]: core primer sequence and barcode counts, if not
            multiplexed the second element is empty
    '''
    if not counter:
        return "", Counter()
    seq_counts = counter.most_common()
    align_map = {
        "3p": "left",   # build 3' consensus from right to left
        "5p": "left",  # build 5' consensus from left to right
    }
    align = align_map.get(orientation, "left")
    core = build_consensus(seq_counts, min_support_frac, n_threshold, align)

    if not core:
        return "", Counter()
    else:
        core = core.strip("N")

    barcode_counts: Counter[str] = Counter()
    if multiplexed:
        if "N" in core:
            idx_start = core.find("N")
            idx_end = core.rfind("N") + 1
            for seq, weight in counter.items():
                if len(seq) < len(core):
                    continue
                barcode = seq[idx_start: idx_end]
                barcode_counts[barcode] += weight
        else:
            core_len = len(core)
            for seq, weight in counter.items():
                if len(seq) < core_len:
                    continue
                if orientation == "5p":
                    idx = seq.rfind(core[:10])
                    if idx == -1:
                        barcode = seq[:-core_len]
                    else:
                        barcode = seq[:idx]
                else:
                    barcode = seq[core_len:]

                barcode_counts[barcode] += weight

        barcode_counts = collapse_barcodes(barcode_counts)
    return core, barcode_counts

def collapse_barcodes(
    barcode_counts: Counter[str],
    *,
    top_n: int = 15,
    selected_target_fraction: float = 0.9,
    min_candidate_fraction: float = 0.005,
) -> Counter[str]:
    """Collapse shorter/noisier barcodes onto high-confidence representatives.

    The heuristic follows the user-specified procedure: determine the dominant
    barcode length from the top ``top_n`` sequences, keep those as anchors, and
    iteratively absorb other barcodes whose nearest anchor (by edit distance)
    is unique while either overall anchor coverage remains below
    ``selected_target_fraction`` or the candidate itself exceeds
    ``min_candidate_fraction`` of observations.
    """

    if not barcode_counts:
        return barcode_counts
    if edlib is None:
        raise RuntimeError(
            "edlib is required for barcode collapsing; install edlib (e.g. pip install edlib)."
        )

    total = sum(barcode_counts.values())
    if total == 0:
        return barcode_counts

    top_items = barcode_counts.most_common(top_n)
    if not top_items:
        return barcode_counts

    length_scores: defaultdict[int, int] = defaultdict(int)
    for barcode, count in top_items:
        length_scores[len(barcode)] += count

    if not length_scores:
        return barcode_counts

    # Prefer the most supported length; break ties in favour of longer barcodes.
    target_length = max(length_scores.items(), key=lambda kv: (kv[1], kv[0]))[0]

    anchors = [barcode for barcode, _ in top_items if len(barcode) == target_length]
    if not anchors:
        anchors = [barcode for barcode, _ in top_items]

    anchors_set = set(anchors)
    result = Counter()
    selected_total = 0
    for barcode in anchors:
        count = barcode_counts[barcode]
        result[barcode] += count
        selected_total += count

    # Process remaining barcodes from most to least abundant.
    for barcode, count in barcode_counts.most_common():
        if barcode in anchors_set:
            continue

        candidate_fraction = count / total
        distances = []
        for anchor in anchors:
            alignment = edlib.align(barcode, anchor, mode="NW", task="distance")
            dist = alignment.get("editDistance")
            if dist == -1:
                continue
            distances.append((dist, anchor))

        if not distances:
            result[barcode] += count
            continue

        distances.sort(key=lambda pair: (pair[0], pair[1]))
        min_dist, best_anchor = distances[0]

        # Require a unique best anchor (no ties on distance).
        if len(distances) > 1 and distances[0][0] == distances[1][0]:
            result[barcode] += count
            continue

        selected_fraction = selected_total / total if total else 0.0
        if selected_fraction < selected_target_fraction or candidate_fraction > min_candidate_fraction:
            result[best_anchor] += count
            selected_total += count
        else:
            result[barcode] += count

    return result

#########################
#     I/O utilities     #
#########################
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