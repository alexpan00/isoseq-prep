'''
Script that reads a bam file and adds the zm tag to each read based on its
read name. The zm is the number after the first slash in the read name.
'''

import pysam
import argparse
import os
import re
import sys
import shutil
import subprocess
import tempfile


def parse_args():
	p = argparse.ArgumentParser(
		description="Add/replace `zm` integer tag in a BAM based on read name."
	)
	p.add_argument("input", help="Input BAM file (indexed or unindexed)")
	p.add_argument("-o", "--output", help="Output BAM path. If omitted, will create '<input>.zm.bam'")
	p.add_argument("-f", "--force", action="store_true", help="Overwrite output if it exists")
	p.add_argument("-v", "--verbose", action="store_true", help="Print progress every 100k reads")
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

	if shutil.which('samtools') is None:
		print("Error: samtools executable not found in PATH; required to convert FASTQ inputs.", file=sys.stderr)
		sys.exit(1)

	tmp_fd, tmp_path = tempfile.mkstemp(prefix="add_zm_", suffix=".bam")
	os.close(tmp_fd)
	cmd = [
		'samtools',
		'import',
		'-o', tmp_path,
		in_path,
	]
	if verbose:
		print(f"Converting FASTQ to BAM via: {' '.join(cmd)}", file=sys.stderr)
	try:
		subprocess.run(cmd, check=True)
	except subprocess.CalledProcessError as exc:
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


def main():
	args = parse_args()

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

				zm = extract_zm_from_qname(qname)
				if zm is not None:
					aln.set_tag('zm', zm, value_type='i')
				else:
					aln.set_tag('zm', count, value_type='i')

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
	main()

