#!/usr/bin/env python3
"""
Infer likely adapter sequences by searching for polyA tails (or polyT heads)
in reads and collecting the sequence region downstream (after polyA) or
upstream (before polyT) as adapter candidates. Produces a simple consensus
adapter and top candidate sequences seen in the sample.

Supports input as BAM (.bam) or FASTQ (.fastq/.fq, optionally gzipped).

Usage examples:
    python scripts/infer_adapters.py input.bam --sample 5000 --min-poly 10
    python scripts/infer_adapters.py input.bam --multiplexed --max-adapter-len 45

Note: this is a heuristic tool intended to find obvious polyA-adapter
signatures; adjust parameters (sample size, poly length, search window) to
match your data.

When ``--multiplexed`` is supplied the script additionally reports shared
primer cores and barcode counts for each end.
"""

import argparse
import gzip
import os
import re
from collections import Counter
from dataclasses import dataclass
from typing import Iterable, List, Optional, Tuple

try:
    import pysam
except Exception:
    pysam = None

# Maximum gap (in bases) between separate polyA matches to merge into one
# e.g. AAAAAAATAAAAA -> with MERGE_GAP=1 these will be merged into a single tail
MERGE_GAP = 1
# Fraction of weighted 'N' at a column above which the position is excluded
# from the consensus (e.g., 0.1 == 10%)
N_THRESHOLD = 0.1
DEFAULT_SUPPORT_FRAC = 0.6


def revcomp_seq(s: str) -> str:
    tr = str.maketrans('ATCGN', 'TAGCN')
    return s.translate(tr)[::-1]


def parse_args():
    p = argparse.ArgumentParser(description="Infer adapter from polyA/polyT signatures in reads")
    p.add_argument("input", help="Input BAM or FASTQ file")
    p.add_argument("--sample", type=int, default=2000, help="Number of reads to sample (default: 2000)")
    p.add_argument("--min-poly", type=int, default=10, help="Minimum consecutive A/Ts to call polyA/polyT (default: 10)")
    p.add_argument("--search-window", type=int, default=200, help="Search window from end/start for polyA/polyT (default: 200 nt)")
    p.add_argument("--max-adapter-len", "--max-primer-len", dest="max_adapter_len", type=int, default=30, help="Maximum adapter/primer length to record per read (default: 30)")
    p.add_argument("--multiplexed", action="store_true", help="Enable multiplex primer mode to infer shared cores and barcodes")
    p.add_argument("--min-support", type=float, default=DEFAULT_SUPPORT_FRAC, help="Minimum support fraction for consensus columns (default: 0.6)")
    p.add_argument("--n-threshold", type=float, default=N_THRESHOLD, help="Maximum 'N' fraction allowed per consensus column (default: 0.1)")
    p.add_argument("--top", type=int, default=20, help="Report at most this many barcode sequences per end when multiplexed (default: 20)")
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()


def read_fastq_seqs(path: str, limit: int) -> Iterable[str]:
    opener = gzip.open if path.lower().endswith('.gz') else open
    with opener(path, 'rt') as fh:
        i = 0
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().strip()
            fh.readline()  # +
            fh.readline()  # qual
            yield seq
            i += 1
            if i >= limit:
                break


def read_bam_seqs(path: str, limit: int) -> Iterable[str]:
    if pysam is None:
        raise RuntimeError("pysam is required to read BAM files; install pysam or provide FASTQ input")
    with pysam.AlignmentFile(path, "rb", check_sq=False) as bam:
        i = 0
        for aln in bam.fetch(until_eof=True):
            seq = aln.query_sequence
            if seq:
                yield seq
                i += 1
            if i >= limit:
                break


def find_polyA_candidates(seq: str, min_poly: int, search_window: int, max_adapter_len: int) -> Tuple[Optional[str], Optional[str]]:
    """Return the adapter sequence (if any) and hit type ('A' or 'T')."""
    seq = seq.upper()
    # helper: reverse complement
    def revcomp(s: str) -> str:
        tr = str.maketrans('ATCGN', 'TAGCN')
        return s.translate(tr)[::-1]

    # polyA: look in trailing window first. If found, do NOT search polyT.
    tail = seq[-search_window:] if len(seq) >= search_window else seq
    # find all polyA runs within the tail
    matches = [m for m in re.finditer(r'A{%d,}' % min_poly, tail)]
    if matches:
        # convert to full-seq coordinates
        regions = []
        base_off = len(seq) - len(tail)
        for m in matches:
            regions.append([base_off + m.start(), base_off + m.end()])
        # merge nearby regions using MERGE_GAP
        merged = []
        cur_start, cur_end = regions[0]
        for s, e in regions[1:]:
            if s - cur_end <= MERGE_GAP:
                # extend current region
                cur_end = max(cur_end, e)
            else:
                merged.append((cur_start, cur_end))
                cur_start, cur_end = s, e
        merged.append((cur_start, cur_end))
        # choose the last merged region (closest to read end)
        start, end = merged[-1]
        # adapter region is downstream of the polyA (within the read)
        if end < len(seq):
            adapter = seq[end:end + max_adapter_len]
            if adapter:
                return adapter, 'A'

    # If no polyA found, look for polyT in the leading window and apply the
    # same merging logic as for polyA. If multiple non-mergeable polyT regions
    # exist, report the first (most 5' end) region and take the upstream
    # sequence as the adapter (reverse-complemented to match polyA orientation).
    head = seq[:search_window] if len(seq) >= search_window else seq
    matches_t = [m for m in re.finditer(r'T{%d,}' % min_poly, head)]
    if matches_t:
        regions_t = []
        for m in matches_t:
            regions_t.append([m.start(), m.end()])
        # merge nearby regions using MERGE_GAP
        merged_t = []
        cur_start, cur_end = regions_t[0]
        for s, e in regions_t[1:]:
            if s - cur_end <= MERGE_GAP:
                cur_end = max(cur_end, e)
            else:
                merged_t.append((cur_start, cur_end))
                cur_start, cur_end = s, e
        merged_t.append((cur_start, cur_end))
        # choose the first merged region (closest to 5' end)
        start, end = merged_t[0]
        # adapter region is upstream/before the polyT (i.e., prefix immediately before polyT)
        if start > 0:
            s = max(0, start - max_adapter_len)
            adapter = seq[s:start]
            if adapter:
                # reverse-complement the upstream adapter so it is in the same
                # orientation as adapters found downstream of polyA
                return revcomp(adapter), 'T'
    return None, None


def _build_weighted_consensus(seq_counts: List[Tuple[str, int]], min_support_frac: float = 0.5,
                              n_threshold: float = 0.1, align: str = 'right') -> str:
    """Generalized weighted consensus builder.

    seq_counts: list of (sequence, count) pairs.
    align: 'right' for 3' consensus (pad left with N), 'left' for 5' consensus
           (pad right with N).
    The algorithm walks from the anchor end (right for 'right', left for
    'left') and collects bases until a column fails the N-threshold or
    the min_support_frac among non-N weighted votes.
    """
    if not seq_counts:
        return ""
    maxlen = max(len(s) for s, _ in seq_counts)
    padded = []
    for s, cnt in seq_counts:
        if align == 'right':
            pad = 'N' * (maxlen - len(s)) + s
        else:
            pad = s + 'N' * (maxlen - len(s))
        padded.append((pad, cnt))

    # choose iteration order depending on alignment
    if align == 'right':
        indices = range(maxlen - 1, -1, -1)
        collect_rev = True
    else:
        indices = range(0, maxlen)
        collect_rev = False

    collected = []
    for i in indices:
        col_counter = Counter()
        total = 0
        n_weight = 0
        for seq, cnt in padded:
            base = seq[i]
            col_counter[base] += cnt
            total += cnt
            if base == 'N':
                n_weight += cnt
        if total == 0:
            break
        if (n_weight / total) > n_threshold:
            break
        # exclude N for majority calculation
        col_counter_noN = Counter({b: w for b, w in col_counter.items() if b != 'N'})
        if not col_counter_noN:
            break
        base, count = col_counter_noN.most_common(1)[0]
        non_n_total = total - n_weight
        if non_n_total <= 0:
            break
        if (count / non_n_total) < min_support_frac:
            break
        collected.append(base)

    if collect_rev:
        return ''.join(reversed(collected))
    return ''.join(collected)


def build_consensus(seq_counts: List[Tuple[str, int]], min_support_frac: float = 0.5, n_threshold: float = 0.1) -> str:
    return _build_weighted_consensus(seq_counts, min_support_frac=min_support_frac, n_threshold=n_threshold, align='right')


def build_consensus_left(seq_counts: List[Tuple[str, int]], min_support_frac: float = 0.5, n_threshold: float = 0.1) -> str:
    return _build_weighted_consensus(seq_counts, min_support_frac=min_support_frac, n_threshold=n_threshold, align='left')


def build_consensus_padded(
    seq_counts: List[Tuple[str, int]],
    min_support_frac: float,
    n_threshold: float,
    align: str,
) -> str:
    """Consensus builder that keeps walking after weak columns by emitting 'N'."""
    if not seq_counts:
        return ""
    maxlen = max(len(s) for s, _ in seq_counts)
    padded = []
    for seq, weight in seq_counts:
        if align == 'right':
            pad = 'N' * (maxlen - len(seq)) + seq
        else:
            pad = seq + 'N' * (maxlen - len(seq))
        padded.append((pad, weight))

    if align == 'right':
        indices = range(maxlen - 1, -1, -1)
        reverse_out = True
    else:
        indices = range(maxlen)
        reverse_out = False

    collected: List[str] = []
    for idx in indices:
        column = Counter()
        total = 0
        n_weight = 0
        for seq, weight in padded:
            base = seq[idx]
            column[base] += weight
            total += weight
            if base == 'N':
                n_weight += weight
        if total == 0:
            collected.append('N')
            continue
        if (n_weight / total) > n_threshold:
            collected.append('N')
            continue
        non_n = {base: weight for base, weight in column.items() if base != 'N'}
        if not non_n:
            collected.append('N')
            continue
        base, count = Counter(non_n).most_common(1)[0]
        support = total - n_weight
        if support <= 0:
            collected.append('N')
            continue
        if (count / support) < min_support_frac:
            collected.append('N')
            continue
        collected.append(base)

    result = ''.join(reversed(collected)) if reverse_out else ''.join(collected)
    return result


def compute_core_and_barcodes(
    counter: Counter[str],
    orientation: str,
    min_support_frac: float,
    n_threshold: float,
) -> Tuple[str, Counter[str]]:
    if not counter:
        return "", Counter()
    seq_counts = counter.most_common()
    align = 'left' if orientation == '3p' else 'right'
    core = build_consensus_padded(seq_counts, min_support_frac, n_threshold, align)
    core = core.strip('N') if core else ""
    if not core:
        return "", Counter()
    core_len = len(core)
    barcode_counts: Counter[str] = Counter()
    for seq, weight in counter.items():
        if len(seq) < core_len:
            continue
        if orientation == '5p':
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
        if not self.core:
            print(f"No reliable {label.lower()} primer core detected.")
            return
        print(f"{label} primer core: {self.core}")
        if not self.barcode_counts:
            print(f"No barcode diversity detected for {label.lower()} primer.")
            return
        total = sum(self.barcode_counts.values())
        print(f"Top barcodes for {label.lower()} primer (n={total} observations):")
        for barcode, count in self.barcode_counts.most_common(limit):
            name = barcode if barcode else "<none>"
            frac = count / total
            print(f"  {name}\t{count}\t{frac:.2%}")


def analyze_multiplexed_primers(
    three_prime_counts: Counter[str],
    five_prime_counts: Counter[str],
    min_support_frac: float,
    n_threshold: float,
    report_limit: int,
) -> None:
    three_core, three_barcodes = compute_core_and_barcodes(three_prime_counts, '3p', min_support_frac, n_threshold)
    five_core, five_barcodes = compute_core_and_barcodes(five_prime_counts, '5p', min_support_frac, n_threshold)

    PrimerStats(three_core, three_barcodes).report("3'", report_limit)
    PrimerStats(five_core, five_barcodes).report("5'", report_limit)


def main():
    args = parse_args()

    path = args.input
    sample = args.sample
    min_poly = args.min_poly
    search_window = args.search_window
    max_adapter_len = args.max_adapter_len

    seq_iter = None
    if path.lower().endswith('.bam'):
        if pysam is None:
            raise SystemExit("pysam not available; install pysam to read BAM files or provide FASTQ input")
        seq_iter = read_bam_seqs(path, sample)
    elif any(path.lower().endswith(ext) for ext in ('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
        seq_iter = read_fastq_seqs(path, sample)
    else:
        # try BAM first, fallback to FASTQ
        if pysam is not None and os.path.exists(path):
            try:
                seq_iter = read_bam_seqs(path, sample)
            except Exception:
                seq_iter = read_fastq_seqs(path, sample)
        else:
            seq_iter = read_fastq_seqs(path, sample)

    adapter_counts = Counter()
    five_prime_counts = Counter()
    total_reads = 0
    polyA_hits = 0
    polyT_hits = 0
    sampled_seqs = []

    for seq in seq_iter:
        seq = seq.strip()
        sampled_seqs.append(seq)
        total_reads += 1
        adapter, hit = find_polyA_candidates(seq, min_poly, search_window, max_adapter_len)
        if adapter:
            c2 = adapter.rstrip('N')
            if c2:
                adapter_counts[c2] += 1
            if hit == 'A':
                polyA_hits += 1
            cand5 = seq[:max_adapter_len] if len(seq) >= max_adapter_len else seq
            cand5 = cand5.rstrip('N')
            if cand5:
                five_prime_counts[cand5] += 1
            elif hit == 'T':
                polyT_hits += 1
            tail = seq[-max_adapter_len:] if len(seq) >= max_adapter_len else seq
            cand5 = revcomp_seq(tail)
            cand5 = cand5.rstrip('N')
            if cand5:
                five_prime_counts[cand5] += 1
        if args.verbose and total_reads % 1000 == 0:
            print(f"Processed {total_reads} reads, adapter candidates so far: {len(adapter_counts)}")

    # summarize
    print(f"Reads sampled: {total_reads}")
    print(f"PolyA-like hits (tail): {polyA_hits}")
    print(f"PolyT-like hits (head): {polyT_hits}")

    if not adapter_counts:
        print("No adapter candidates found with the current parameters.")
        return

    # show top candidates
    print("\nTop adapter candidate sequences:")
    top = adapter_counts.most_common(10)
    for s, cnt in top:
        print(f"  {s}\t{cnt}")

    # pass (sequence, count) pairs to build_consensus so consensus is weighted
    consensus = build_consensus(top, min_support_frac=args.min_support, n_threshold=args.n_threshold)
    if consensus:
        print(f"\nConsensus adapter (3' end, from top candidates): {consensus}")
    else:
        print("\nNo strong column-wise consensus could be built from top candidates.")

    # Now attempt to infer the 5' adapter consensus of the same length
    L = len(consensus)
    if L == 0:
        return

    adapter5_counts = Counter()
    for seq in sampled_seqs:
        # detect if this read had a polyA or polyT (reuse detection)
        _, hit = find_polyA_candidates(seq, min_poly, search_window, max_adapter_len)
        if not hit:
            continue
        # for reads with polyA at the tail, take the first L bases
        if hit == 'A':
            if len(seq) >= L:
                cand = seq[:L].upper()
                adapter5_counts[cand] += 1
        elif hit == 'T':
            # for polyT reads, extract last L bases and reverse-complement them
            if len(seq) >= L:
                tail = seq[-L:].upper()
                cand = revcomp_seq(tail)
                adapter5_counts[cand] += 1

    if not adapter5_counts:
        print("\nNo 5' adapter candidates found (insufficient data).")
        return

    print("\nTop 5' adapter candidate sequences:")
    top5 = adapter5_counts.most_common(1000)
    i = 0
    while i < len(top5) and i < 5:
        s, cnt = top5[i]
        print(f"  {s}\t{cnt}")
        i += 1

    consensus5 = build_consensus_left(top5, min_support_frac=args.min_support, n_threshold=args.n_threshold)
    if consensus5:
        print(f"\nConsensus adapter (5' end): {consensus5}")
    else:
        print("\nNo strong 5' consensus could be built from top candidates.")

    if args.multiplexed:
        analyze_multiplexed_primers(
            adapter_counts,
            five_prime_counts,
            args.min_support,
            args.n_threshold,
            args.top,
        )


if __name__ == '__main__':
    main()
