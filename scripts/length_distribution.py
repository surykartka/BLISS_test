"""
Script to generate statistics of number of alignments for reads of different lengths.
"""

from __future__ import print_function
import sys, os

def get_hist(d):
	vs = d.values()
	for v in range(min(vs), max(vs)+1):
		n = len([x for x in d if d[x] == v])
		print(v, n)

fastq1 = sys.argv[1]
fastq2 = sys.argv[2]
sam = sys.argv[3]

lengths = {}
for i, line in enumerate(open(fastq1)):
	n = i % 4
	if n == 1:
		seq = line.strip()
		length = len(seq)
		lengths[seqid] = length
	elif n == 0:
		seqid = ':'.join(line.split()[0].split(':')[-2:]) + ' 1'
get_hist(lengths)
for i, line in enumerate(open(fastq2)):
	n = i % 4
	if n == 1:
		seq = line.strip()
		length = len(seq)
		lengths[seqid] = length
	elif n == 0:
		seqid = ':'.join(line.split()[0].split(':')[-2:]) + ' 2'
print(''); get_hist(lengths)

alns = {}
for line in open(sam):
	if line.startswith('@'):
		continue
	tab = line.split()
	seqid = ':'.join(tab[0].split(':')[-2:])
	if int(tab[1]) > 127:
		seqid += ' 2'
	else:
		seqid += ' 1'
	alns[seqid] = alns.get(seqid, 0) + 1


print([x for x in alns.items() if x[1] > 1 and x[0].endswith(' 1')])


