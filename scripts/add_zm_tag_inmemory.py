"""Single-pass variant of add_zm_tag that reheaders once reads are written.

This script streams the input BAM exactly once, tagging reads and collecting the
set of PU values as it writes the modified alignments to a temporary BAM. After
the pass completes it computes the updated header and rewrites only the header
section of the temporary BAM, yielding the final output file. This lets you
compare performance with the baseline two-pass implementation while avoiding a
second read of the source BAM.
"""

import argparse
import os
import re
import sys
import tempfile
from typing import Dict, Optional, Set

import pysam


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Add/replace `zm` integer tag in a BAM based on read name (single-pass variant)."
    )
    parser.add_argument("input", help="Input BAM file (indexed or unindexed)")
    parser.add_argument(
        "-o",
        "--output",
        help="Output BAM path. If omitted, will create '<input>.inmem.zm.bam'",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Overwrite output if it exists",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print progress every 100k reads",
    )
    return parser.parse_args()


def extract_zm_from_qname(qname: str) -> Optional[int]:
    """Return the integer after the first '/' in the read name or None."""
    if not qname:
        return None
    match = re.match(r"^[^/]+/(\d+)", qname)
    if not match:
        return None
    try:
        return int(match.group(1))
    except ValueError:
        return None


def derive_pu_from_qname(qname: str) -> str:
    """Return the prefix before the first '/' in the read name or 'unknown'."""
    if not qname:
        return "unknown"
    if "/" in qname:
        prefix = qname.split("/", 1)[0]
        return prefix or "unknown"
    return qname


def prepare_rg_state(header_dict: Dict) -> Dict[str, str]:
    """Normalize existing RG entries and build the PU→RG map."""
    rg_list = list(header_dict.get("RG", []) or [])
    header_dict["RG"] = rg_list
    pu_to_rg: Dict[str, str] = {}
    used_ids: Set[str] = set()
    next_index = 0

    def next_rg_id() -> str:
        nonlocal next_index
        while True:
            candidate = f"RG{next_index}"
            next_index += 1
            if candidate not in used_ids:
                used_ids.add(candidate)
                return candidate

    for rg in rg_list:
        rg_id = rg.get("ID")
        if not rg_id:
            rg_id = next_rg_id()
            rg["ID"] = rg_id
        used_ids.add(rg_id)
        pu_value = rg.get("PU") or "unknown"
        rg["PU"] = pu_value
        pu_to_rg.setdefault(pu_value, rg_id)
        ds_value = rg.get("DS", "") or ""
        if "READTYPE=CCS" not in ds_value:
            separator = "" if not ds_value or ds_value.endswith(";") else ";"
            rg["DS"] = f"{ds_value}{separator}READTYPE=CCS" if ds_value else "READTYPE=CCS"
        rg.setdefault("SM", pu_value)

    header_dict.setdefault("RG", rg_list)
    header_dict.setdefault("PG", header_dict.get("PG", []))
    header_dict.setdefault("SQ", header_dict.get("SQ", []))

    header_dict["__used_rg_ids"] = used_ids
    header_dict["__next_rg_index"] = next_index
    return pu_to_rg


def assign_rg(header_dict: Dict, pu_to_rg: Dict[str, str], pu_value: str) -> str:
    """Return an RG ID for the given PU value, creating a new one as needed."""
    if pu_value in pu_to_rg:
        return pu_to_rg[pu_value]

    used_ids: Set[str] = header_dict["__used_rg_ids"]
    next_index: int = header_dict["__next_rg_index"]

    while True:
        rg_id = f"RG{next_index}"
        next_index += 1
        if rg_id not in used_ids:
            used_ids.add(rg_id)
            break

    header_dict["__next_rg_index"] = next_index

    rg_entry = {"ID": rg_id, "PU": pu_value, "DS": "READTYPE=CCS", "SM": pu_value}
    header_dict["RG"].append(rg_entry)
    pu_to_rg[pu_value] = rg_id
    return rg_id


def main() -> None:
    args = parse_args()

    in_path = args.input
    if args.output:
        out_path = args.output
    else:
        base = in_path[:-4] if in_path.lower().endswith(".bam") else in_path
        out_path = f"{base}.inmem.zm.bam"

    if os.path.exists(out_path):
        if not args.force:
            print(
                f"Error: output file '{out_path}' already exists. Use --force to overwrite.",
                file=sys.stderr,
            )
            sys.exit(1)
        os.remove(out_path)

    tmp_dir = os.path.dirname(os.path.abspath(out_path)) or "."
    tmp_handle = tempfile.NamedTemporaryFile(prefix="add_zm_", suffix=".bam", dir=tmp_dir, delete=False)
    tmp_path = tmp_handle.name
    tmp_handle.close()

    total_reads = 0
    missing_zm_counter = 0
    unique_pus: Set[str] = set()

    with pysam.AlignmentFile(in_path, "rb", check_sq=False) as inh, pysam.AlignmentFile(
        tmp_path, "wb", template=inh
    ) as outh:
        header_obj = inh.header
        header_dict = header_obj.to_dict() if hasattr(header_obj, "to_dict") else dict(header_obj)
        pu_to_rg = prepare_rg_state(header_dict)

        for total_reads, aln in enumerate(inh.fetch(until_eof=True), start=1):
            qname = aln.query_name or ""
            pu = derive_pu_from_qname(qname)
            unique_pus.add(pu)

            rg_id = assign_rg(header_dict, pu_to_rg, pu)
            aln.set_tag("RG", rg_id, value_type="Z")

            zm = extract_zm_from_qname(qname)
            if zm is not None:
                aln.set_tag("zm", zm, value_type="i")
            else:
                aln.set_tag("zm", missing_zm_counter, value_type="i")
                missing_zm_counter += 1

            outh.write(aln)
            if args.verbose and total_reads % 100000 == 0:
                print(f"Processed {total_reads} reads...", file=sys.stderr)

    if not unique_pus:
        unique_pus.add("unknown")
        assign_rg(header_dict, pu_to_rg, "unknown")

    header_dict.pop("__used_rg_ids", None)
    header_dict.pop("__next_rg_index", None)

    new_header = pysam.AlignmentHeader.from_dict(header_dict)

    if args.verbose:
        print(f"Identified {len(unique_pus)} PU value(s)", file=sys.stderr)

    with pysam.AlignmentFile(tmp_path, "rb", check_sq=False) as interim, pysam.AlignmentFile(
        out_path, "wb", header=new_header
    ) as final_out:
        for idx, aln in enumerate(interim.fetch(until_eof=True), start=1):
            final_out.write(aln)
            if args.verbose and idx % 100000 == 0:
                print(f"Copied {idx} reads during header rewrite...", file=sys.stderr)

    try:
        os.unlink(tmp_path)
    except OSError:
        pass

    if args.verbose:
        print(f"Finished: wrote {out_path} ({total_reads} reads)")


if __name__ == "__main__":
    main()
