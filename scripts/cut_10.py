from __future__ import print_function
import sys

l = 10
infile, outfile = sys.argv[1], sys.argv[2]
with open(outfile, 'w') as f:
	for line in open(infile):
		if line[:1] == '>':
			print(line.strip(), file=f)
		else:
			seq = line.strip()[:10]
			n = len(seq)
			if n < l:
				seq = seq + ''.join((l-n) * ['-'])
			print(seq, file=f)
