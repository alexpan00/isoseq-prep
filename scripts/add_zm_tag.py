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

	If the pattern is missing or not numeric, return None.
	Examples:
	  'm64011_180916_202355/12345/ccs' -> 12345
	  'movie/42' -> 42
	  'no_slash' -> None
	"""
	m = re.match(r'^[^/]+/(\d+)', qname)
	if m:
		try:
			return int(m.group(1))
		except ValueError:
			return None
	return None


def prepare_header(header_dict, pu_value):
	"""Ensure every RG has PU and a READTYPE=CCS descriptor."""
	rg_list = header_dict.get('RG', [])
	if rg_list:
		for rg in rg_list:
			if 'PU' not in rg or not rg['PU']:
				rg['PU'] = pu_value
			ds_value = rg.get('DS', '') or ''
			if 'READTYPE=CCS' not in ds_value:
				separator = '' if not ds_value or ds_value.endswith(';') else ';'
				rg['DS'] = f"{ds_value}{separator}READTYPE=CCS" if ds_value else 'READTYPE=CCS'
	else:
		header_dict['RG'] = [{'ID': 'RG0', 'PU': pu_value, 'DS': 'READTYPE=CCS'}]
	return header_dict


def main():
	args = parse_args()

	in_path = args.input
	out_path = args.output or (in_path + ".zm.bam" if not in_path.lower().endswith('.bam') else in_path[:-4] + ".zm.bam")

	if os.path.exists(out_path) and not args.force:
		print(f"Error: output file '{out_path}' already exists. Use --force to overwrite.", file=sys.stderr)
		sys.exit(1)

	# Single pass through BAM: open once, inspect first read, then stream through
	count = 0
	with pysam.AlignmentFile(in_path, "rb", check_sq=False) as inh:
		header_obj = inh.header
		if hasattr(header_obj, 'to_dict'):
			header_dict = header_obj.to_dict()
		else:
			header_dict = dict(header_obj)

		stream = inh.fetch(until_eof=True)
		first_aln = next(stream, None)

		pu = 'unknown'
		if first_aln is not None:
			qn = first_aln.query_name or ''
			pu = qn.split('/', 1)[0] if '/' in qn else qn

		header_dict = prepare_header(header_dict, pu)

		with pysam.AlignmentFile(out_path, "wb", header=header_dict) as outh:
			if first_aln is not None:
				qname = first_aln.query_name
				zm = extract_zm_from_qname(qname)
				if zm is not None:
					first_aln.set_tag('zm', zm, value_type='i')
				else:
					first_aln.set_tag('zm', count, value_type='i')
				outh.write(first_aln)
				count += 1

			for aln in stream:
				qname = aln.query_name
				zm = extract_zm_from_qname(qname)
				if zm is not None:
					aln.set_tag('zm', zm, value_type='i')
				else:
					aln.set_tag('zm', count, value_type='i')
				outh.write(aln)
				count += 1
				if args.verbose and count % 100000 == 0:
					print(f"Processed {count} reads...", file=sys.stderr)

	if args.verbose:
		print(f"Finished: wrote {out_path} ({count} reads)")


if __name__ == '__main__':
	main()

