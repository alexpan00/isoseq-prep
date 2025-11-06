#!/usr/bin/env python3
"""Infer multiplexed Iso-Seq primer sequences.

This script extends the heuristics from ``infer_adapters.py`` to handle
multiplexed primers where several barcodes share a common primer core. It:

* samples reads from a BAM/FASTQ file
* finds polyA/polyT signals to locate adapter/primer segments
* aggregates both 3' (downstream of polyA) and 5' (upstream of polyT) segments
* computes weighted consensus sequences for the shared primer core
* reports barcode sequences (prefix/suffix surrounding the core) with counts
* optionally renders a knee plot summarising barcode frequency distributions

The goal is to distinguish the common primer portion from barcode-specific
sequence, supporting cases where one or both ends are multiplexed.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import accumulate
from typing import Iterable, List, Optional, Tuple

try:
    import pysam
except Exception:  # pragma: no cover - optional dependency
    pysam = None

try:
    import edlib  # type: ignore
except Exception:  # pragma: no cover - optional dependency
    edlib = None

try:
    import matplotlib.pyplot as plt  # type: ignore
except Exception:  # pragma: no cover - optional dependency
    plt = None

MERGE_GAP = 1
N_THRESHOLD = 0.1
DEFAULT_SUPPORT_FRAC = 0.6
MIN_POLYA_HIT_FRAC = 0.8  # minimum fraction of reads with polyA/polyT hits


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
        "--multiplexed",
        action="store_true",
        help="Enable multiplexed primer detection",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=20,
        help="Report at most this many barcode sequences per end (default: 20)",
    )
    parser.add_argument(
        "--plot-prefix",
        default=None,
        help=(
            "Base filename (without extension) for the barcode knee plot. "
            "Defaults to <input_basename>_primer_knee.png if omitted."
        ),
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Enable barcode knee plot generation.",
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
            if adapter:
                return adapter, "A"

    head = seq[:search_window] if len(seq) >= search_window else seq
    matches_t = list(re.finditer(r"T{%d,}" % min_poly, head))
    if matches_t:
        regions_t: List[Tuple[int, int]] = [(m.start(), m.end()) for m in matches_t]
        merged_t: List[Tuple[int, int]] = merge_regions(regions_t)
        start, end = merged_t[0]
        if start > 0:
            # the adapter is everything upstream of the polyT stretch
            adapter = seq[:start]
            if adapter:
                return revcomp_seq(adapter), "T"
    return None, None


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


def build_consensus(
    seq_counts: List[Tuple[str, int]],
    min_support_frac: float,
    n_threshold: float,
    align: str,
) -> str:
    '''
    Wrapper to _build_weighted_consensus

    Args:
        seq_counts (List[Tuple[str, int]]): Counter of adapter sequences
        min_support_frac (float): minimum support fraction for consensus
        n_threshold (float): N threshold for consensus
        align (str): alignment direction ("left" or "right")

    Returns:
        str: consensus sequence
    '''
    return _build_weighted_consensus(seq_counts, min_support_frac, n_threshold, align)


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
        Tuple[str, Counter[str]]: core primer sequence and barcode counts
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
) -> dict:
    """
    Compute primer cores and barcode counts for 3' and 5' ends and return
    a dictionary with the results. This function does not perform any
    printing; use `show_primer_results` to display the returned data.
    """
    three_core, three_barcodes = compute_core_and_barcodes(
        three_prime_counts, "3p", args.min_support, args.n_threshold, args.multiplexed
    )
    five_core, five_barcodes = compute_core_and_barcodes(
        five_prime_counts, "5p", args.min_support, args.n_threshold, args.multiplexed
    )

    return {
        "3p": {"core": three_core, "barcodes": three_barcodes},
        "5p": {"core": five_core, "barcodes": five_barcodes},
    }


def show_primer_results(results: dict, args: argparse.Namespace) -> None:
    """
    Pretty-print primer analysis results produced by `analyze_primers`.
    """
    three = results.get("3p", {})
    five = results.get("5p", {})

    PrimerStats(three.get("core", ""), three.get("barcodes", Counter())).report("3'", args.top)
    PrimerStats(five.get("core", ""), five.get("barcodes", Counter())).report("5'", args.top)


def save_knee_plot(results: dict, args: argparse.Namespace) -> None:
    """Render a barcode knee plot if data and dependencies are available."""

    series = []
    for label in ("3p", "5p"):
        barcodes: Counter[str] = results.get(label, {}).get("barcodes", Counter())
        if barcodes:
            counts = sorted(barcodes.values(), reverse=True)
            if counts:
                ranks = list(range(1, len(counts) + 1))
                cumulative = list(accumulate(counts))
                total = cumulative[-1]
                cumulative = [value / total for value in cumulative]
                series.append((label, ranks, counts, cumulative))

    if not series:
        print("No barcode counts available for knee plot.", file=sys.stderr)
        return

    prefix = args.plot_prefix
    if not prefix:
        stem = os.path.splitext(os.path.basename(args.input))[0] or "primer"
        prefix = f"{stem}_primer_knee"
    output_path = prefix if prefix.endswith(".png") else f"{prefix}.png"

    fig, ax_counts = plt.subplots()
    ax_cum = ax_counts.twinx()

    count_lines = []
    cumulative_lines = []
    for label, ranks, counts, cumulative in series:
        marker = "o" if len(ranks) <= 20 else None
        (count_line,) = ax_counts.plot(ranks, counts, marker=marker, label=f"{label} counts")
        (cum_line,) = ax_cum.plot(
            ranks,
            cumulative,
            linestyle="--",
            color=count_line.get_color(),
            label=f"{label} cumulative",
        )
        count_lines.append(count_line)
        cumulative_lines.append(cum_line)

    ax_counts.set_xscale("log")
    ax_counts.set_yscale("log")
    ax_counts.set_xlim(1, max(max(ranks) for _, ranks, _, _ in series))
    ax_counts.set_xlabel("Barcode rank (log scale)")
    ax_counts.set_ylabel("Barcode counts (log scale)")
    ax_counts.grid(True, which="both", alpha=0.3)

    ax_cum.set_ylim(0, 1.05)
    ax_cum.set_ylabel("Cumulative fraction of counts")
    ax_cum.axhline(0.9, color="grey", linestyle=":", linewidth=1, alpha=0.7)

    lines = count_lines + cumulative_lines
    labels = [line.get_label() for line in lines]
    ax_counts.legend(lines, labels, loc="best")
    ax_counts.set_title("Primer barcode knee plot")

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    print(f"Knee plot saved to {output_path}")


def main() -> None:
    args = parse_args()

    # Generate sequence iterator
    seq_iter = iterate_sequences(args.input, args.sample)

    # string counter for 5' and 3' end sequences
    adapter_counts: Counter[str] = Counter()
    five_prime_counts: Counter[str] = Counter()
    sampled = 0
    polyA_hits = 0
    polyT_hits = 0

    # sample sequences and extract candidate primers (could be refactored)
    for seq in seq_iter:
        seq = seq.strip().upper()
        if not seq:
            continue
        sampled += 1
        # Find polyA/polyT and extract candidate primers from 3' end
        adapter, hit = find_polyA_candidates(seq, args.min_poly, 
                                             args.search_window, 
                                             args.max_primer_len)
        if adapter:
            trimmed = adapter.rstrip("N")
            if trimmed:
                adapter_counts[trimmed] += 1
        # Adapter found at the 3' end
        if hit == "A":
            polyA_hits += 1
            if len(seq) >= args.max_primer_len:
                cand5 = seq[: args.max_primer_len]
            else:
                cand5 = seq
            cand5 = cand5.rstrip("N")
            if cand5:
                five_prime_counts[cand5] += 1
        # Adapter found at the 5' end
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

    # Initial sampling results, in case potential primers are found end here
    # Could mean that the reads end in polyA or have already been fully-trimmed
    print(f"Reads sampled: {sampled}")
    print(f"PolyA tail hits (3'): {polyA_hits}")
    print(f"PolyT head hits (5'): {polyT_hits}")
    
    # if the proportion of polyA/polyT hits is too low, likely that reads are
    # fully-trimmed (isoseq refine)
    if polyA_hits + polyT_hits < sampled*MIN_POLYA_HIT_FRAC:
        print("PolyA detection rate too low; reads may be fully-trimmed.")
        print("EXITING.")
    # if the sum of adapter counts is very low, likely no primers detected
    elif sum(adapter_counts.values()) + sum(five_prime_counts.values()) < sampled*MIN_POLYA_HIT_FRAC:
        print("No primer candidates detected downstream the polyA")
        print("EXITING.")
    else:
        primers_res = analyze_primers(adapter_counts, five_prime_counts, args)
        show_primer_results(primers_res, args)
        if args.plot:
            save_knee_plot(primers_res, args)

if __name__ == "__main__":
    main()
