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
zeros = set()
mappings = {}
for line in open(sam):
	if line.startswith('@'):
		continue
	tab = line.split()
	seqid = ':'.join(tab[0].split(':')[-2:])
	sam_score = int(tab[3])
	if sam_score == 0:
		zeros.add(seqid)
		continue
	start, end = int(tab[3]), int(tab[7])
	chrom = tab[2]
	if start > end:
		seqid += ' 2'
	else:
		seqid += ' 1'
	alns[seqid] = alns.get(seqid, 0) + 1
	mappings[seqid] = (chrom, start, end)

print('')
print("How many paired-end reads do not map at all?", len(zeros))

one_1 = [x for x in alns.items() if x[1] == 1 and x[0].endswith(' 1')]
print("How many map from 5' only once?", len(one_1))
print("How many map from 5' to more than one place?", len([x for x in alns.items() if x[1] > 1 and x[0].endswith(' 1')]))

one_2 = [x for x in alns.items() if x[1] == 1 and x[0].endswith(' 2')]
print("How many map from 3' only once?", len(one_2))
print("How many map from 3' to more than one place?", len([x for x in alns.items() if x[1] > 1 and x[0].endswith(' 2')]))

both_once = set([x[0].split()[0] for x in one_1]) & set([x[0].split()[0] for x in one_2])
print("How many reads map from both sides only once?", len(both_once))

print("Where do they map?")
chr2length = {} # dictionary mapping chromosomes to their length
for line in open(sam):
	if line.startswith('@SQ') and 'LN:' in line:
		chrom = line.split('SN:')[1].split()[0]
		length = int(line.split('LN:')[1].strip())
		chr2length[chrom] = length
	if line.startswith('@PG'):
		break

nparts = 100
chr2parts = {}
for chrom, length in chr2length.items():
	dec_length = 1.*length/nparts
	chr2parts[chrom] = {}
	for i in range(nparts):
		chr2parts[chrom][(dec_length*i, dec_length*(i+1))] = set()

chr2mil = {}
for seqid_prefix in both_once:
	seqid1 = seqid_prefix + ' 1'
	seqid2 = seqid_prefix + ' 2'
	chrom, start, end = mappings[seqid1]
	#print(lengths[seqid1], lengths[seqid2], mappings[seqid1], mappings[seqid2])
	for (dec_start, dec_end) in chr2parts[chrom]:
		if start >= dec_start and start <= dec_end:
			chr2parts[chrom][(dec_start, dec_end)].add(seqid1)

"""
print('')
for chrom in pos2reads:
	for pos in sorted(pos2reads[chrom]):
		if pos2reads[chrom][pos] < 10:
			continue
	print(chrom, pos, pos2reads[chrom][pos])
"""

pos2reads = {}
for chrom in sorted(chr2parts):
	for dec_bounds in sorted(chr2parts[chrom]):
		if len(chr2parts[chrom][dec_bounds]) > 9:
			for seqid in chr2parts[chrom][dec_bounds]:
				chrom, start, end = mappings[seqid]
				pos2reads[(chrom, start, end)] = pos2reads.get((chrom, start, end), 0) + 1
				#for i in range(start, end+1):
				#	pos2reads[i] = pos2reads.get(i, 0) + 1
			#for i in sorted(pos2reads):
			#	if pos2reads[i] > 9:
			#		print(chrom, i, pos2reads[i])
		#print(chrom, dec_bounds, chr2parts[chrom][dec_bounds])

print('')
for x,n in sorted(pos2reads.items(), key=lambda y: y[1], reverse=True):
	if n > 9:
		print('* [%s:%s-%s](https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=%s%%3A%s-%s)' % (x+x), n, sep=': ')
print('')
#chr18 6258510 429


