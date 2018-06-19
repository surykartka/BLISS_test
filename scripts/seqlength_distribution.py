from __future__ import print_function

import sys
from Bio import SeqIO

fasta = sys.argv[1]
lengths = [len(rec.seq) for rec in SeqIO.parse(open(fasta), 'fasta')]
print(*lengths, sep='\n')
