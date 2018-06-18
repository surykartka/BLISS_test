from __future__ import print_function
import sys

infile, outfile = sys.argv[1], sys.argv[2]
with open(outfile, 'w') as f:
	for line in open(infile):
		if line[:1] == '>':
			print(line.strip(), file=f)
		else:
			print(line.strip()[:10], file=f)
