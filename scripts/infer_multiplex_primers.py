#!/usr/bin/env python3
"""Infer multiplexed Iso-Seq primer sequences.

This script extends the heuristics from ``infer_adapters.py`` to handle
multiplexed primers where several barcodes share a common primer core. It:

* samples reads from a BAM/FASTQ file
* finds polyA/polyT signals to locate adapter/primer segments
* aggregates both 3' (downstream of polyA) and 5' (upstream of polyT) segments
* computes weighted consensus sequences for the shared primer core
* reports barcode sequences (prefix/suffix surrounding the core) with counts

The goal is to distinguish the common primer portion from barcode-specific
sequence, supporting cases where one or both ends are multiplexed.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
from collections import Counter
from dataclasses import dataclass
from typing import Iterable, List, Optional, Tuple

try:
    import pysam
except Exception:  # pragma: no cover - optional dependency
    pysam = None

MERGE_GAP = 1
N_THRESHOLD = 0.1
DEFAULT_SUPPORT_FRAC = 0.6


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
        default=200,
        help="Window length from read end/start to search for polyA/polyT (default: 200)",
    )
    parser.add_argument(
        "--max-primer-len",
        type=int,
        default=45,
        help="Maximum primer length to record from each end (default: 45)",
    )
    parser.add_argument(
        "--min-support",
        type=float,
        default=DEFAULT_SUPPORT_FRAC,
        help="Minimum support fraction for consensus columns (default: 0.6)",
    )
    parser.add_argument(
        "--n-threshold",
        type=float,
        default=N_THRESHOLD,
        help="Maximum fraction of 'N' votes allowed per consensus column (default: 0.1)",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=20,
        help="Report at most this many barcode sequences per end (default: 20)",
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


def find_polyA_candidates(
    seq: str, min_poly: int, search_window: int, max_adapter_len: int
) -> Tuple[Optional[str], Optional[str]]:
    """Return the adapter sequence (if any) and hit type ('A' or 'T')."""
    seq = seq.upper()

    def revcomp(s: str) -> str:
        table = str.maketrans("ATCGN", "TAGCN")
        return s.translate(table)[::-1]

    tail = seq[-search_window:] if len(seq) >= search_window else seq
    matches = list(re.finditer(r"A{%d,}" % min_poly, tail))
    if matches:
        base_offset = len(seq) - len(tail)
        regions: List[Tuple[int, int]] = []
        for m in matches:
            regions.append((base_offset + m.start(), base_offset + m.end()))
        merged: List[Tuple[int, int]] = []
        start, end = regions[0]
        for s, e in regions[1:]:
            if s - end <= MERGE_GAP:
                end = max(end, e)
            else:
                merged.append((start, end))
                start, end = s, e
        merged.append((start, end))
        start, end = merged[-1]
        if end < len(seq):
            adapter = seq[end : end + max_adapter_len]
            if adapter:
                return adapter, "A"

    head = seq[:search_window] if len(seq) >= search_window else seq
    matches_t = list(re.finditer(r"T{%d,}" % min_poly, head))
    if matches_t:
        regions_t: List[Tuple[int, int]] = [(m.start(), m.end()) for m in matches_t]
        merged_t: List[Tuple[int, int]] = []
        start, end = regions_t[0]
        for s, e in regions_t[1:]:
            if s - end <= MERGE_GAP:
                end = max(end, e)
            else:
                merged_t.append((start, end))
                start, end = s, e
        merged_t.append((start, end))
        start, end = merged_t[0]
        if start > 0:
            s_idx = max(0, start - max_adapter_len)
            adapter = seq[s_idx:start]
            if adapter:
                return revcomp(adapter), "T"
    return None, None


def _build_weighted_consensus(
    seq_counts: List[Tuple[str, int]],
    min_support_frac: float,
    n_threshold: float,
    align: str,
) -> str:
    if not seq_counts:
        return ""
    maxlen = max(len(s) for s, _ in seq_counts)
    padded = []
    for seq, weight in seq_counts:
        if align == "right":
            pad = "N" * (maxlen - len(seq)) + seq
        else:
            pad = seq + "N" * (maxlen - len(seq))
        padded.append((pad, weight))

    if align == "right":
        indices = range(maxlen - 1, -1, -1)
        collect_rev = True
    else:
        indices = range(maxlen)
        collect_rev = False

    collected: List[str] = []
    for idx in indices:
        counter = Counter()
        total = 0
        n_weight = 0
        for seq, weight in padded:
            base = seq[idx]
            counter[base] += weight
            total += weight
            if base == "N":
                n_weight += weight
        if total == 0:
            collected.append("N")
            continue
        if n_weight / total > n_threshold:
            collected.append("N")
            continue
        non_n = {base: w for base, w in counter.items() if base != "N"}
        if not non_n:
            collected.append("N")
            continue
        base, count = Counter(non_n).most_common(1)[0]
        support = total - n_weight
        if support <= 0:
            collected.append("N")
            continue
        if count / support < min_support_frac:
            collected.append("N")
            continue
        collected.append(base)

    return "".join(reversed(collected)) if collect_rev else "".join(collected)


def build_consensus(
    seq_counts: List[Tuple[str, int]],
    min_support_frac: float,
    n_threshold: float,
    align: str,
) -> str:
    return _build_weighted_consensus(seq_counts, min_support_frac, n_threshold, align)


def compute_core_and_barcodes(
    counter: Counter[str],
    orientation: str,
    min_support_frac: float,
    n_threshold: float,
) -> Tuple[str, Counter[str]]:
    if not counter:
        return "", Counter()
    seq_counts = counter.most_common()
    align_map = {
        "3p": "left",   # build 3' consensus from left to right
        "5p": "right",  # build 5' consensus from right to left
    }
    align = align_map.get(orientation, "left")
    core = build_consensus(seq_counts, min_support_frac, n_threshold, align)
    if not core:
        core = "N"
    else:
        core = core.strip("N")
    if not core:
        return "", Counter()
    core_len = len(core)
    barcode_counts: Counter[str] = Counter()
    for seq, weight in counter.items():
        if len(seq) < core_len:
            continue
        if orientation == "5p":
            idx = seq.rfind(core)
            if idx == -1:
                barcode = seq[:-core_len]
            else:
                barcode = seq[:idx]
        else:
            idx = seq.find(core)
            if idx == -1:
                barcode = seq[core_len:]
            else:
                barcode = seq[idx + core_len :]
        barcode_counts[barcode] += weight
    return core, barcode_counts


@dataclass
class PrimerStats:
    core: str
    barcode_counts: Counter[str]

    def report(self, label: str, limit: int) -> None:
        print(f"\n{label} primer core length: {len(self.core)}")
        if self.core:
            print(f"{label} primer core: {self.core}")
        else:
            print(f"No reliable {label.lower()} primer core detected.")
            return
        if not self.barcode_counts:
            print(f"No barcode diversity detected for {label.lower()} primer.")
            return
        total = sum(self.barcode_counts.values())
        print(f"Top barcodes for {label.lower()} primer (n={total} observations):")
        for barcode, count in self.barcode_counts.most_common(limit):
            name = barcode if barcode else "<none>"
            frac = count / total
            print(f"  {name}\t{count}\t{frac:.2%}")


def analyze_primers(
    three_prime_counts: Counter[str],
    five_prime_counts: Counter[str],
    args: argparse.Namespace,
) -> None:
    three_core, three_barcodes = compute_core_and_barcodes(
        three_prime_counts, "3p", args.min_support, args.n_threshold
    )
    five_core, five_barcodes = compute_core_and_barcodes(
        five_prime_counts, "5p", args.min_support, args.n_threshold
    )

    PrimerStats(three_core, three_barcodes).report("3'", args.top)
    PrimerStats(five_core, five_barcodes).report("5'", args.top)


def main() -> None:
    args = parse_args()

    seq_iter = iterate_sequences(args.input, args.sample)

    adapter_counts: Counter[str] = Counter()
    five_prime_counts: Counter[str] = Counter()
    sampled = 0
    polyA_hits = 0
    polyT_hits = 0

    for seq in seq_iter:
        seq = seq.strip().upper()
        if not seq:
            continue
        sampled += 1
        adapter, hit = find_polyA_candidates(seq, args.min_poly, args.search_window, args.max_primer_len)
        if adapter:
            trimmed = adapter.rstrip("N")
            if trimmed:
                adapter_counts[trimmed] += 1
        if hit == "A":
            polyA_hits += 1
            if len(seq) >= args.max_primer_len:
                cand5 = seq[: args.max_primer_len]
            else:
                cand5 = seq
            cand5 = cand5.rstrip("N")
            if cand5:
                five_prime_counts[cand5] += 1
        elif hit == "T":
            polyT_hits += 1
            if len(seq) >= args.max_primer_len:
                tail = seq[-args.max_primer_len :]
            else:
                tail = seq
            cand5 = revcomp_seq(tail)
            cand5 = cand5.rstrip("N")
            if cand5:
                five_prime_counts[cand5] += 1
        if args.verbose and sampled % 1000 == 0:
            print(
                f"Processed {sampled} reads | 3' adapters: {len(adapter_counts)} | 5' candidates: {len(five_prime_counts)}"
            )

    print(f"Reads sampled: {sampled}")
    print(f"PolyA tail hits (3'): {polyA_hits}")
    print(f"PolyT head hits (5'): {polyT_hits}")
    print(five_prime_counts)
    if not adapter_counts and not five_prime_counts:
        print("No primer candidates detected with the current parameters.")
        return

    analyze_primers(adapter_counts, five_prime_counts, args)


if __name__ == "__main__":
    main()
