"""Add/replace the per-read `zm` tag in a BAM and enrich related header/tags.

This utility reads a BAM (or a FASTQ which it will convert to BAM using
`samtools import`) and writes a new BAM where every alignment has a numeric
`zm:i:<N>` tag derived from the read name. The `zm` value is taken from the
integer immediately following the first ``/`` in the read name (for example
``movie/12345/ccs`` -> ``zm:i:12345``). When the pattern is not present the
script assigns a monotonically increasing integer as a fallback.

Additional behavior implemented by the script:
- Enriches the BAM header `@RG` (read-group) lines to ensure each observed
	platform-unit (PU) is present. New `RG` records are created as needed and
	receive `DS` metadata containing `READTYPE=CCS` to satisfy `isoseq refine`.
- Optionally computes an `rq` tag (read quality / accuracy) from base
	qualities with `--compute_quality` and writes it as a float tag when
	`rq` is missing. Otherwise a default `rq` value is written (`--default-quality`).
- Optionally add `np`tag when number of passes is missing. This is useful for
	sc data, where this tag is used in the groupedup step.
- Supports `--force` to overwrite existing output and `--verbose` for progress
	messages.

Dependencies and notes:
- Requires `pysam` to read/write BAM files. If the input is FASTQ the script
	uses the `samtools` binary to perform a temporary conversion to BAM; ensure
	`samtools` is available in `PATH` when supplying FASTQ input.
- The script preserves secondary/supplementary alignments behavior and writes
	updated RG/zm/rq tags per alignment.

Example:
		python scripts/add_zm_tag.py movie.bam -o movie.zm.bam --compute_quality --verbose

"""

import pysam
import argparse
import os
import re
import sys
import tempfile

DEFAULT_QUALITY = 0.99
DEFAULT_NP = 10


def add_args(p: argparse.ArgumentParser) -> None:
	p.add_argument("input", help="Input BAM file")
	p.add_argument("-o", 
                "--output", 
                help="Output BAM path. If omitted, will create '<input>.zm.bam'")
	p.add_argument("-f", 
                "--force", 
                action="store_true", 
                help="Overwrite output if it exists")
	group = p.add_mutually_exclusive_group()
	group.add_argument("--default-quality", 
                type=float, default=DEFAULT_QUALITY, 
                help="Default RQ value to assign if missing (default: 0.99) or --compute_quality flag is not used")
	group.add_argument("--compute_quality", 
                action="store_true", 
                help="Compute basewise accuracy from quality scores to assign rq tag if missing")
	p.add_argument("--np-tag",
                type=int, 
                default=DEFAULT_NP, 
                help="If provided, add np tag with this value when missing. (default: 10)")
	p.add_argument("-v", 
                "--verbose", 
                action="store_true", 
                help="Print progress every 100k reads")


def parse_args():
	p = argparse.ArgumentParser(
		description="Add/replace `zm` integer tag in a BAM based on read name."
	)
	add_args(p)
	return p.parse_args()


def extract_zm_from_qname(qname):
	"""Return the integer after the first '/' in the read name.

	If the pattern is missing or not numeric, return None.
	Examples:
	  'm64011_180916_202355/12345/ccs' -> 12345
	  'movie/42' -> 42
	  'no_slash' -> None
	"""
	if not qname:
		return None
	m = re.match(r'^[^/]+/(\d+)', qname)
	if m:
		try:
			return int(m.group(1))
		except ValueError:
			return None
	return None


def derive_pu_from_qname(qname):
	"""Return the prefix before the first '/' in the read name or 'unknown'."""
	qn = qname or ''
	pu = qn.split('/', 1)[0] if '/' in qn else qn
	return pu if pu else 'unknown'


def collect_header_and_pus(in_path):
	"""Scan the BAM once to collect header info and observed PU values."""
	pu_values = set()
	with pysam.AlignmentFile(in_path, "rb", check_sq=False) as inh:
		header_obj = inh.header
		if hasattr(header_obj, 'to_dict'):
			header_dict = header_obj.to_dict()
		else:
			header_dict = dict(header_obj)
		for aln in inh.fetch(until_eof=True):
			pu = derive_pu_from_qname(aln.query_name)
			pu_values.add(pu)
	return header_dict, pu_values


def ensure_bam_input(in_path, verbose=False):
	"""Return a BAM path, converting FASTQ inputs via samtools import when needed."""
	if in_path.lower().endswith('.bam'):
		return in_path, None

	tmp_fd, tmp_path = tempfile.mkstemp(prefix="add_zm_", suffix=".bam")
	os.close(tmp_fd)

	if verbose:
		print(f"Converting FASTQ to BAM", file=sys.stderr)
	try:
		pysam.samtools.fqimport("-o", tmp_path, in_path, catch_stdout=False)
	except Exception as exc:
		print(f"Error: samtools import failed with exit code {exc.returncode}.", file=sys.stderr)
		if os.path.exists(tmp_path):
			os.remove(tmp_path)
		sys.exit(1)

	return tmp_path, tmp_path


def prepare_header(header_dict: dict, pu_values: set) -> tuple[dict, dict]:
	"""Add RG lines for new PU values while keeping existing entries simple.
	Args:
		header_dict: BAM header as a dictionary.
		pu_values: Set of observed PU values.
	Returns:
		Tuple of (updated_header_dict, pu_to_rg_dict)
  	"""
	rg_list = list(header_dict.get('RG', []) or [])
	pu_to_rg = {}
	used_ids = set()
	next_index = len(rg_list)
 
	# check if the existing RG entries have IDs and PUs
	for rg in rg_list:
		rg_id = rg.get('ID')
		# if there is not id create one
		if not rg_id:
			rg_id = f"RG{next_index}"
			while rg_id in used_ids:
				next_index += 1
				rg_id = f"RG{next_index}"
			next_index += 1
			rg['ID'] = rg_id
		if rg_id:
			used_ids.add(rg_id)

		# if there is not PU, assign one from the set
		if rg.get('PU') is None:
			pu_value = pu_values.pop()
			rg['PU'] = pu_value
		pu_to_rg.setdefault(pu_value, rg_id)
		# Ensure DS tag includes READTYPE=CCS, otherwise isoseq refine will complain
		ds_value = rg.get('DS', '') or ''
		if 'READTYPE=CCS' not in ds_value:
			separator = '' if not ds_value or ds_value.endswith(';') else ';'
			rg['DS'] = f"{ds_value}{separator}READTYPE=CCS" if ds_value else 'READTYPE=CCS'
		rg['SM'] = pu_value

	# Now add RG entries for any remaining PU values
	while pu_values:
		pu = pu_values.pop()
		rg_id = f"RG{next_index}"
		while rg_id in used_ids:
			next_index += 1
			rg_id = f"RG{next_index}"
		rg_list.append({'ID': rg_id, 'SM':pu, 'PU': pu, 'DS': 'READTYPE=CCS'})
		pu_to_rg[pu] = rg_id
		used_ids.add(rg_id)
		next_index += 1

	header_dict['RG'] = rg_list
	return header_dict, pu_to_rg


def compute_quality_metrics(aln: pysam.AlignedSegment, default_rq: float = DEFAULT_QUALITY)-> float:
	"""Return accuracy_basewise or placeholder if no quals.
	
	To calculate basewise accuracy from quality values we need to convert
	the Phred scores to probabilities by using the formula:
		P(error) = 10^(-Q/10)
	Then we can compute accuracy as:
	accuracy = 1 - (sum of error probabilities / length of sequence)
	
	Args:
		aln: pysam AlignedSegment object
		default_rq: default accuracy value to return if no quality scores are present
	Returns:
		accuracy_basewise: float
	"""
	quals = aln.query_qualities  # list of ints, or None
	if quals is None:
		return default_rq
	probs = [10 ** (-q / 10.0) for q in quals]
	L = len(probs)
	expected_errors = sum(probs)
	accuracy_basewise = 1.0 - (expected_errors / L) if L > 0 else 0.0
	return accuracy_basewise


def tag_getter_setter(aln: pysam.AlignedSegment, tag: str, value, value_type: str) -> None:
	"""Utility to get/set tags in pysam AlignedSegment objects."""
	try:
		aln.get_tag(tag)
	except KeyError:
		aln.set_tag(tag, value, value_type=value_type)

def main(args: argparse.Namespace) -> None:
	# I/O paths
	in_path = args.input
	out_path = args.output or (in_path + ".zm.bam" if not in_path.lower().endswith('.bam') else in_path[:-4] + ".zm.bam")
	if os.path.exists(out_path) and not args.force:
		print(f"Error: output file '{out_path}' already exists. Use --force to overwrite.", file=sys.stderr)
		sys.exit(1)
	if os.path.exists(out_path) and args.force:
		os.remove(out_path)

	converted_path = None
	try:
		bam_path, converted_path = ensure_bam_input(in_path, verbose=args.verbose)

		# First pass: collect header and PU values
		header_dict, pu_values = collect_header_and_pus(bam_path)
		header_dict, pu_to_rg = prepare_header(header_dict, pu_values)

		if args.verbose and pu_values:
			unique_pus = len(pu_values)
			print(f"Identified {unique_pus} PU value(s)", file=sys.stderr)

		count = 0
		# Second pass: rewrite alignments with the enriched header and tags.
		with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as inh, pysam.AlignmentFile(out_path, "wb", header=header_dict) as outh:
			for aln in inh.fetch(until_eof=True):
				qname = aln.query_name or ''
				pu = derive_pu_from_qname(qname)
				rg_id = pu_to_rg.get(pu)
				if rg_id:
					aln.set_tag('RG', rg_id, value_type='Z')
				# ZM tag
				if not (zm := extract_zm_from_qname(qname)):
					zm = count
				tag_getter_setter(aln, 'zm', zm, 'i')

				# Read quality tag
				accuracy_basewise = compute_quality_metrics(aln, default_rq=args.default_quality) if args.compute_quality else args.default_quality
				tag_getter_setter(aln, "rq", accuracy_basewise, "f")

				# Number of passes tag
				tag_getter_setter(aln, "np", args.np_tag, "i")
				outh.write(aln)
				count += 1
				if args.verbose and count % 100000 == 0:
					print(f"Processed {count} reads...", file=sys.stderr)

	finally:
		if converted_path and os.path.exists(converted_path):
			os.remove(converted_path)

	if args.verbose:
		print(f"Finished: wrote {out_path} ({count} reads)")


if __name__ == '__main__':
	main(parse_args())

