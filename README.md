# BLISS test

Analysis of BLISS test runs performed in June 2018. 

The commands used can be found in [`cmd.sh`](cmd.sh).

## Basic questions

1. How many reads start with the barcode ('AGACTCT'/'AGAGTCT')?
2. How many reads can be mapped to a reference human genome?
3. How many reads can be mapped to telomeres?

## Input data

The input fastq files were downloaded from [here](http://bio4.cent.uw.edu.pl/BCHLT) to the `fastq` directory, but are also accessible on the server (`bigram:/ngs/Fastq/new_BCHLT/`). They correspond to sequencing results from two samples.

## Quality control

* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) was run to check quality of the reads (outputs written to [`FastQC`](FastQC)). Overrepresented sequences suggest contamination with 'RNA PCR Primer, Index 6' (R1 C1_S1), 'RNA PCR Primer, Index 12' (R1 NB_S2), and 'Illumina RNA PCR Primer' (R2). Sequences of these adapter (their reverse complement) were found [here](https://github.com/csf-ngs/fastqc/blob/master/Contaminants/contaminant_list.txt) and written to [`primers.fa`](trim_adapters/primers.fa). (*As the first 'T' could overlap with the last 'T' of the barcode sequence, it was stripped from the adapter sequences.* –– *actually, it wasn't*)

* The adapters were trimmed with [cutadapt](https://cutadapt.readthedocs.io/) (allowing at minimum three bases match between adapter and match, and discarding proccesed reads that are shorter than 6 nt), e.g.: `cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG --minimum-length 6 -o trim_adapters/B_SC-BLESS_C1_S1_L001_R1_BCHLT_cutadapt_1.fastq.gz fastq/B_SC-BLESS_C1_S1_L001_R1_BCHLT.fastq.gz > analysis/cutadapt.out`:

		=== Summary ===

		Total reads processed:               5,253,177
		Reads with adapters:                 5,038,644 (95.9%)
		Reads that were too short:             737,275 (14.0%)
		Reads written (passing filters):     4,515,902 (86.0%)

		Total basepairs processed:   372,975,567 bp
		Total written (filtered):     56,381,610 bp (15.1%)

* Then, only reads with `AGACTCT` at 5' end were left, e.g.: `cutadapt -g ^AGACTCT -e 0 --overlap 7 --minimum-length 1 -o trim_adapters/B_SC-BLESS_C1_S1_L001_R1_BCHLT_cutadapt_2.fastq.gz --discard-untrimmed trim_adapters/B_SC-BLESS_C1_S1_L001_R1_BCHLT_cutadapt.fastq.gz >> analysis/cutadapt.out`:

		=== Summary ===

		Total reads processed:               4,515,902
		Reads with adapters:                   107,137 (2.4%)
		Reads that were too short:               4,412 (0.1%)
		Reads written (passing filters):       102,725 (2.3%)
		
		Total basepairs processed:    56,381,610 bp
		Total written (filtered):      4,685,925 bp (8.3%)

* Reads with barcodes with max. one mismatch were also filtered, e.g.: `cutadapt -g ^AGACTCT -e 0.15 --overlap 7 --minimum-length 1 -o trim_adapters/B_SC-BLESS_C1_S1_L001_R1_BCHLT_cutadapt_3.fastq.gz --discard-untrimmed trim_adapters/B_SC-BLESS_C1_S1_L001_R1_BCHLT_cutadapt.fastq.gz >> analysis/cutadapt.out`

* Similarly for the remaining files (or run `bash cmd.sh &`), write standard output to [`cutadapt.out`](analysis/cutadapt.out).

## How many reads start with the barcode?

Script [`get_cutadapt_stats.py`](scripts/get_cutadapt_stats.py) writes a table with trimming statistics: `python scripts/get_cutadapt_stats.py`. A plot below (from `analysis/cutadapt_stats.numbers`) shows numbers of reads after trimming and with barcodes for the respective fastq files.

![alt text](analysis/cutadapt_stats.png)
