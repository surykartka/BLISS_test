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

* The adapters were trimmed with [cutadapt](https://cutadapt.readthedocs.io/) (allowing at minimum three bases match between adapter and match, and discarding proccesed reads that are shorter than 6 nt), and only the reads without matching adapters were left, e.g.: `cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG --discard-trimmed --minimum-length 6 -o trim_adapters/B_SC-BLESS_C1_S1_L001_R1_BCHLT_cutadapt_1.fastq.gz fastq/B_SC-BLESS_C1_S1_L001_R1_BCHLT.fastq.gz > analysis/cutadapt.out`:

		=== Summary ===

		Total reads processed:               5,253,177
		Reads with adapters:                 5,038,644 (95.9%)
		Reads that were too short:             737,275 (14.0%)
		Reads written (passing filters):       214,533 (4.1%)
		
		Total basepairs processed:   372,975,567 bp
		Total written (filtered):     15,231,843 bp (4.1%)

* Then, only reads with `AGACTCT` at 5' end were left, e.g.: `cutadapt -g ^AGACTCT -e 0 --overlap 7 --minimum-length 1 -o trim_adapters/B_SC-BLESS_C1_S1_L001_R1_BCHLT_cutadapt_2.fastq.gz --discard-untrimmed trim_adapters/B_SC-BLESS_C1_S1_L001_R1_BCHLT_cutadapt.fastq.gz >> analysis/cutadapt.out`:

		=== Summary ===

		Total reads processed:                 214,533
		Reads with adapters:                    44,853 (20.9%)
		Reads that were too short:                   0 (0.0%)
		Reads written (passing filters):        44,853 (20.9%)
		
		Total basepairs processed:    15,231,843 bp
		Total written (filtered):      2,870,592 bp (18.8%)

* Reads with barcodes with max. one mismatch were also filtered, e.g.: `cutadapt -g ^AGACTCT -e 0.15 --overlap 7 --minimum-length 1 -o trim_adapters/B_SC-BLESS_C1_S1_L001_R1_BCHLT_cutadapt_3.fastq.gz --discard-untrimmed trim_adapters/B_SC-BLESS_C1_S1_L001_R1_BCHLT_cutadapt.fastq.gz >> analysis/cutadapt.out`

* Similarly for the remaining files (or run `bash cmd.sh &`), write standard output to [`cutadapt.out`](analysis/cutadapt.out).

## How many reads start with the barcode?

Script [`get_cutadapt_stats.py`](scripts/get_cutadapt_stats.py) writes a table with trimming statistics: `python scripts/get_cutadapt_stats.py`. A plot below (from `analysis/cutadapt_stats.numbers`) shows numbers of reads after trimming and with barcodes for the respective fastq files.

![alt text](analysis/cutadapt_stats.png)

* Large contamination with adapters (over 90% of reads are affected).
* Small fraction of reads with perfect match to the barcode at 5' end.

File | Total reads | Without adapters (min 6 nt) | With barcode (min 1 nt, max one mismatch) | With barcode (min 1 nt)
-----|-------------|-----------------------------------|--------------------------------------------------------------|--------------------------------------
C1_S1_L001_R1 | 5,253,177 | 214,533 | 172,057 | 44,853
C1_S1_L001_R2 | 5,253,177 | 853,167 | 4,141 | 88
NB_S2_L001_R1 | 3,793,031 | 927,052 | 775,657 | 143,443
NB_S2_L001_R2 | 3,793,031 | 958,800 | 3,527 | 55


## How many of the barcoded reads can be mapped to human genome?

Mapping of the barcoded reads to the human genome was done as follows:

* The reference indexed GRCh38/hg38 genome was downloaded from [NCBI ftp](ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/) (`GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz`), as linked in [the Bowtie 2 website](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).


