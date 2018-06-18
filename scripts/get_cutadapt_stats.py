from __future__ import print_function

"""This scripts reads cutadapt std output and writes the basic statistics to a tab-delimited table.
The table is then visualized in Numbers and plotted to cutadapt_stats.png.
"""

cutadapt_output = 'analysis/cutadapt.out' # std output of cutadapt
output = 'analysis/cutadapt.stats' # output with final table

with open(output, 'w') as f:
	print('File', 'Total reads', 'Reads at least 6 nt long after adapter trimming', 'Reads with barcode of non-zero length (one mismatch allowed)', 'Reads with barcode of non-zero length', sep='\t', file=f)
	for line in open(cutadapt_output):
		if line.startswith('Command line parameters'):
			infile = line.split()[-1].split('/')[-1]
			outfile = line.split()[-2].split('/')[-1]
			if outfile == '--discard-untrimmed':
				outfile = line.split()[-3].split('/')[-1]
			print(infile, outfile)
		if line.startswith('Total reads processed'):
			if 'cutadapt' not in infile:
				total = line.split()[-1].replace(',','')
		if line.startswith('Reads written (passing filters)'):
			n = line.split()[-2].replace(',','')
			if 'cutadapt' not in infile:
				after_trim = n
			elif 'cutadapt_2' in outfile:
				barcode_7 = n
			else:
				barcode_6 = n
				print(infile, total, after_trim, barcode_6, barcode_7, sep='\t', file=f)

print('Stats written to', output)