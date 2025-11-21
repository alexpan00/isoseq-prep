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
* automatically detects the knee to recommend adapter counts and FASTA output

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
from typing import Iterable, List, Optional, Sequence, Tuple
from utils import build_consensus, revcomp_seq, iterate_sequences, compute_core_and_barcodes

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


def find_knee_point(counts: Sequence[int]) -> Tuple[int, float]:
    """Identify the knee point using a simple Kneedle-style heuristic.

    Returns the rank (1-based) of the knee and the cumulative fraction covered
    by adapters up to that knee. If no counts are provided, returns (0, 0.0).
    """

    if not counts:
        return 0, 0.0

    total = sum(counts)
    if total == 0:
        return 0, 0.0

    if len(counts) == 1:
        return 1, 1.0

    cumulative = list(accumulate(counts))
    n = len(counts)
    diffs = []
    for idx, cum in enumerate(cumulative):
        x_norm = idx / (n - 1)
        y_norm = cum / total
        diffs.append(y_norm - x_norm)

    knee_idx = max(range(n), key=lambda i: diffs[i])
    if diffs[knee_idx] <= 0:
        knee_idx = 0

    knee_rank = knee_idx + 1
    knee_fraction = cumulative[knee_idx] / total
    return knee_rank, knee_fraction


def adjust_knee_for_empty(
    sorted_counts: Sequence[Tuple[str, int]],
    cumulative: Sequence[int],
    knee_rank: int,
) -> Tuple[int, float]:
    """Nudge the knee rank left if it selects an empty adapter ("" barcode)."""

    if not sorted_counts or knee_rank <= 0:
        return 0, 0.0

    knee_rank = min(knee_rank, len(sorted_counts))
    total = cumulative[-1] if cumulative else 0
    if total == 0:
        return 0, 0.0

    while knee_rank > 1 and sorted_counts[knee_rank - 1][0] == "":
        knee_rank -= 1

    knee_fraction = cumulative[knee_rank - 1] / total
    return knee_rank, knee_fraction


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

@dataclass
class PrimerStats:
    core: str
    barcode_counts: Counter[str]
    sorted_barcodes: List[Tuple[str, int]]
    knee_rank: int = 0
    knee_fraction: float = 0.0

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
        if self.knee_rank > 0 and total > 0:
            print(
                f"Knee suggests {self.knee_rank} barcodes covering {self.knee_fraction:.2%} of reads."
            )
        barcodes_iter = self.sorted_barcodes if self.sorted_barcodes else self.barcode_counts.most_common()
        print(f"Top barcodes for {label.lower()} primer (n={total} observations):")
        for barcode, count in barcodes_iter[:limit]:
            name = barcode if barcode else "<none>"
            frac = count / total if total else 0.0
            print(f"  {name}\t{count}\t{frac:.2%}")


def analyze_primers(
    three_prime_counts: Counter[str],
    five_prime_counts: Counter[str],
    *,
    min_support: float,
    n_threshold: float,
    multiplexed: bool,
) -> dict:
    """Compute primer summaries, optionally inferring barcode knees."""

    three_core, three_barcodes = compute_core_and_barcodes(
        three_prime_counts, "3p", min_support, n_threshold, multiplexed
    )
    five_core, five_barcodes = compute_core_and_barcodes(
        five_prime_counts, "5p", min_support, n_threshold, multiplexed
    )

    # find the number of barcodes to report based on knee point only if multiplexed
    if multiplexed:
        # 3' end analysis
        three_sorted = three_barcodes.most_common()
        three_counts = [count for _, count in three_sorted]
        three_cumulative = list(accumulate(three_counts)) if three_counts else []
        three_knee_rank, three_knee_fraction = find_knee_point(three_counts)
        three_knee_rank, three_knee_fraction = adjust_knee_for_empty(
            three_sorted, three_cumulative, three_knee_rank
        )

        # 5' end analysis
        five_sorted = five_barcodes.most_common()
        five_counts = [count for _, count in five_sorted]
        five_cumulative = list(accumulate(five_counts)) if five_counts else []
        five_knee_rank, five_knee_fraction = find_knee_point(five_counts)
        five_knee_rank, five_knee_fraction = adjust_knee_for_empty(
            five_sorted, five_cumulative, five_knee_rank
        )

    else:
        three_sorted = three_cumulative = []
        three_knee_rank = 0
        three_knee_fraction = 0.0
        five_sorted = five_cumulative = []
        five_knee_rank = 0
        five_knee_fraction = 0.0

    return {
        "3p": {
            "core": three_core,
            "barcodes": three_barcodes,
            "sorted": three_sorted,
            "knee_rank": three_knee_rank,
            "knee_fraction": three_knee_fraction,
            "total": three_cumulative[-1] if three_cumulative else 0,
        },
        "5p": {
            "core": five_core,
            "barcodes": five_barcodes,
            "sorted": five_sorted,
            "knee_rank": five_knee_rank,
            "knee_fraction": five_knee_fraction,
            "total": five_cumulative[-1] if five_cumulative else 0,
        },
    }


def show_primer_results(results: dict, top_n: int) -> None:
    """
    Pretty-print primer analysis results produced by `analyze_primers`.
    """
    three = results.get("3p", {})
    five = results.get("5p", {})

    PrimerStats(
        three.get("core", ""),
        three.get("barcodes", Counter()),
        three.get("sorted", []),
        three.get("knee_rank", 0),
        three.get("knee_fraction", 0.0),
    ).report("3'", top_n)
    PrimerStats(
        five.get("core", ""),
        five.get("barcodes", Counter()),
        five.get("sorted", []),
        five.get("knee_rank", 0),
        five.get("knee_fraction", 0.0),
    ).report("5'", top_n)


def save_knee_plot(results: dict, output_path: str) -> None:
    """Render a barcode knee plot if data and dependencies are available."""

    if plt is None:
        print("matplotlib is not available; skipping barcode knee plot.", file=sys.stderr)
        return

    series = []
    for label in ("3p", "5p"):
        data = results.get(label, {})
        sorted_counts: List[Tuple[str, int]] = data.get("sorted", []) or []
        if not sorted_counts:
            continue
        counts = [count for _, count in sorted_counts]
        ranks = list(range(1, len(counts) + 1))
        cumulative_raw = list(accumulate(counts))
        total = cumulative_raw[-1]
        cumulative = [value / total for value in cumulative_raw]
        series.append(
            {
                "label": label,
                "ranks": ranks,
                "counts": counts,
                "cumulative": cumulative,
                "knee_rank": data.get("knee_rank", 0),
                "knee_fraction": data.get("knee_fraction", 0.0),
            }
        )

    if not series:
        print("No barcode counts available for knee plot.", file=sys.stderr)
        return

    fig, ax_counts = plt.subplots()
    ax_cum = ax_counts.twinx()

    count_lines = []
    cumulative_lines = []
    for entry in series:
        label = entry["label"]
        ranks = entry["ranks"]
        counts = entry["counts"]
        cumulative = entry["cumulative"]
        knee_rank = entry["knee_rank"]
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

        if knee_rank and 1 <= knee_rank <= len(ranks):
            ax_counts.axvline(
                knee_rank,
                color=count_line.get_color(),
                linestyle=":",
                linewidth=1,
                alpha=0.7,
            )
            ax_cum.scatter(
                [knee_rank],
                [cumulative[knee_rank - 1]],
                color=count_line.get_color(),
                marker="x",
                zorder=5,
            )

    ax_counts.set_xscale("log")
    ax_counts.set_yscale("log")
    ax_counts.set_xlim(1, max(max(entry["ranks"]) for entry in series))
    ax_counts.set_xlabel("Barcode rank (log scale)")
    ax_counts.set_ylabel("Barcode counts (log scale)")
    ax_counts.grid(True, which="both", alpha=0.3)

    ax_cum.set_ylim(0, 1.05)
    ax_cum.set_ylabel("Cumulative fraction of counts")
    ax_cum.axhline(0.9, color="grey", linestyle=":", linewidth=1, alpha=0.7)

    lines = count_lines + cumulative_lines
    labels = [line.get_label() for line in lines]
    ax_counts.legend(
        lines,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.2),
        ncol=2,
        frameon=False,
    )
    ax_counts.set_title("Primer barcode knee plot")

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    print(f"Knee plot saved to {output_path}")


def write_adapter_fasta(results: dict, output_path: str) -> None:
    """Write adapter FASTA sequences up to the detected knee."""

    if not output_path:
        raise ValueError("output_path must be provided for FASTA export")

    entries: List[Tuple[str, str]] = []

    for label in ("5p", "3p"):
        data = results.get(label, {})
        core = data.get("core", "")
        if not core:
            continue
        sorted_counts: List[Tuple[str, int]] = data.get("sorted", []) or []
        total = data.get("total", 0) or sum(count for _, count in sorted_counts)

        if not sorted_counts:
            header = f">{label}_adapter_core"
            sequence = core
            entries.append((header, sequence))
            continue

        knee_rank = data.get("knee_rank", len(sorted_counts)) or len(sorted_counts)
        knee_rank = max(1, min(len(sorted_counts), knee_rank))

        cumulative = 0
        for idx, (barcode, count) in enumerate(sorted_counts[:knee_rank], 1):
            barcode_seq = barcode or ""
            if label == "5p":
                if not "N" in core:
                    sequence = barcode_seq + core
                else:
                    # replace N by barcode
                    seq_list = list(core)
                    idx_start = core.find("N")
                    idx_end = core.rfind("N") + 1
                    seq_list[idx_start: idx_end] = list(barcode_seq)
                    sequence = "".join(seq_list)
            else:
                if not "N" in core:
                    sequence = core + barcode_seq
                else:
                    # replace N by barcode
                    seq_list = list(core)
                    idx_start = core.find("N")
                    idx_end = core.rfind("N") + 1
                    seq_list[idx_start: idx_end] = list(barcode_seq)
                    sequence = "".join(seq_list)
            cumulative += count
            frac = count / total if total else 0.0
            cum_frac = cumulative / total if total else 0.0
            header = (
                f">{label}_adapter_{idx}|count={count}|freq={frac:.4f}|cum={cum_frac:.4f}"
            )
            entries.append((header, sequence))

    if not entries:
        print("No adapter sequences available for FASTA output.", file=sys.stderr)
        return

    with open(output_path, "w", encoding="ascii") as fh:
        for header, sequence in entries:
            fh.write(header + "\n")
            fh.write(sequence + "\n")

    print(f"Adapter FASTA saved to {output_path}")


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
        primers_res = analyze_primers(
            adapter_counts,
            five_prime_counts,
            min_support=args.min_support,
            n_threshold=args.n_threshold,
            multiplexed=args.multiplexed,
        )
        show_primer_results(primers_res, args.top)

        if args.plot and args.multiplexed:
            plot_path = resolve_output_path(
                args.input,
                args.plot_prefix,
                default_suffix="primer_knee",
                default_extension=".png",
                allowed_extensions=(".png",),
            )
            save_knee_plot(primers_res, plot_path)

        if not args.no_fasta:
            fasta_path = resolve_output_path(
                args.input,
                args.fasta_prefix,
                default_suffix="primers",
                default_extension=".fasta",
                allowed_extensions=(".fa", ".fasta", ".fna"),
            )
            write_adapter_fasta(primers_res, fasta_path)

if __name__ == "__main__":
    main()
