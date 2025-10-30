'''
Script that reads a bam file and adds the zm tag to each read based on its
read name. The zm is the number after the first slash in the read name.
'''

import pysam
import argparse
import os
import re
import sys


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

	If the pattern is missing or not numeric, return 0.
	Examples:
	  'm64011_180916_202355/12345/ccs' -> 12345
	  'movie/42' -> 42
	  'no_slash' -> 0
	"""
	m = re.match(r'^[^/]+/(\d+)', qname)
	if m:
		try:
			return int(m.group(1))
		except ValueError:
			return 0
	return 0


def main():
	args = parse_args()

	in_path = args.input
	out_path = args.output or (in_path + ".zm.bam" if not in_path.lower().endswith('.bam') else in_path[:-4] + ".zm.bam")

	if os.path.exists(out_path) and not args.force:
		print(f"Error: output file '{out_path}' already exists. Use --force to overwrite.", file=sys.stderr)
		sys.exit(1)

	# Determine PU from first read (the part before the first '/')
	pu = 'unknown'
	first_aln = None
	with pysam.AlignmentFile(in_path, "rb", check_sq=False) as inh_probe:
		try:
			first_aln = next(inh_probe.fetch(until_eof=True))
		except (ValueError, StopIteration):
			# empty file / no alignments / fetch raised error
			first_aln = None

	if first_aln is not None:
		qn = first_aln.query_name
		pu = qn.split('/', 1)[0] if '/' in qn else qn

	# Prepare/modify header to ensure RG entries have PU
	with pysam.AlignmentFile(in_path, "rb", check_sq=False) as inh_header:
		if hasattr(inh_header.header, 'to_dict'):
			header_dict = inh_header.header.to_dict()
		else:
			# fallback: try converting to a plain dict
			header_dict = dict(inh_header.header)

	# Ensure 'RG' list exists and each entry has 'PU'
	rg_list = header_dict.get('RG', [])
	if rg_list:
		for rg in rg_list:
			if 'PU' not in rg or not rg['PU']:
				rg['PU'] = pu
	else:
		# add a minimal RG entry
		header_dict['RG'] = [{'ID': 'RG0', 'PU': pu}]

	# Now process alignments and write output with modified header
	count = 0
	with pysam.AlignmentFile(in_path, "rb", check_sq=False) as inh:
		with pysam.AlignmentFile(out_path, "wb", header=header_dict) as outh:
			for aln in inh.fetch(until_eof=True):
				qname = aln.query_name
				zm = extract_zm_from_qname(qname)
				# set or replace tag 'zm' as integer
				try:
					aln.set_tag('zm', zm, value_type='i')
				except TypeError:
					# older pysam versions accept just set_tag(tag, value)
					aln.set_tag('zm', zm)
				outh.write(aln)
				count += 1
				if args.verbose and count % 100000 == 0:
					print(f"Processed {count} reads...", file=sys.stderr)

	if args.verbose:
		print(f"Finished: wrote {out_path} ({count} reads)")


if __name__ == '__main__':
	main()

