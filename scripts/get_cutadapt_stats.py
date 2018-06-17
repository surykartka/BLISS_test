from __future__ import print_function

"""This scripts reads cutadapt std output and writes the basic statistics to a tab-delimited table.
The table is then visualized in Numbers and plotted to cutadapt_stats.png.
"""

cutadapt_output = 'analysis/cutadapt.out' # std output of cutadapt
output = 'analysis/cutadapt.stats' # output with final table

with open(output, 'w') as f:
	print('File', 'Total reads', 'Reads at least 6 nt long after adapter trimming', 'Reads with the barcode of non-zero length', sep='\t', file=f)
	for line in open(cutadapt_output):
		if line.startswith('Command line parameters'):
			file = line.split()[-1].split('/')[-1]
		if line.startswith('Total reads processed'):
			if 'cutadapt_1' not in file:
				total = line.split()[-1].replace(',','')
		if line.startswith('Reads written (passing filters)'):
			if 'cutadapt_1' not in file:
				after_trim = line.split()[-2].replace(',','')
			else:
				barcode = line.split()[-2].replace(',','')
				print(file, total, after_trim, barcode, sep='\t', file=f)

print('Stats written to', output)