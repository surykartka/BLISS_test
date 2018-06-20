gunzip -c trim_adapters/B_SC-BLESS_NB_S2_L001_R1_BCHLT_cutadapt_native_adapters_4.fastq.gz > mapping/B_SC-BLESS_NB_S2_L001_R1_BCHLT_cutadapt_native_adapters_4.fastq
gunzip -c trim_adapters/B_SC-BLESS_C1_S1_L001_R1_BCHLT_cutadapt_native_adapters_4.fastq.gz > mapping/B_SC-BLESS_C1_S1_L001_R1_BCHLT_cutadapt_native_adapters_4.fastq
gunzip -c trim_adapters/B_SC-BLESS_NB_S2_L001_R2_BCHLT_cutadapt_native_adapters_4.fastq.gz > mapping/B_SC-BLESS_NB_S2_L001_R2_BCHLT_cutadapt_native_adapters_4.fastq
gunzip -c trim_adapters/B_SC-BLESS_C1_S1_L001_R2_BCHLT_cutadapt_native_adapters_4.fastq.gz > mapping/B_SC-BLESS_C1_S1_L001_R2_BCHLT_cutadapt_native_adapters_4.fastq

fastq_to_fasta -Q33 -i mapping/B_SC-BLESS_NB_S2_L001_R1_BCHLT_cutadapt_native_adapters_4.fastq -o mapping/R1_5prime/C1_S1_R1.fa
fastq_to_fasta -Q33 -i mapping/B_SC-BLESS_C1_S1_L001_R1_BCHLT_cutadapt_native_adapters_4.fastq -o mapping/R1_5prime/NB_S2_R1.fa
fastq_to_fasta -Q33 -i mapping/B_SC-BLESS_NB_S2_L001_R2_BCHLT_cutadapt_native_adapters_4.fastq -o mapping/R2_5prime/C1_S1_R2.fa
fastq_to_fasta -Q33 -i mapping/B_SC-BLESS_C1_S1_L001_R2_BCHLT_cutadapt_native_adapters_4.fastq -o mapping/R2_5prime/NB_S2_R2.fa

python scripts/cut_10.py mapping/R1_5prime/C1_S1_R1.fa mapping/R1_5prime/C1_S1_R1_first10.fa
python scripts/cut_10.py mapping/R1_5prime/NB_S2_R1.fa mapping/R1_5prime/NB_S2_R1_first10.fa
python scripts/cut_10.py mapping/R2_5prime/C1_S1_R2.fa mapping/R2_5prime/C1_S1_R2_first10.fa
python scripts/cut_10.py mapping/R2_5prime/NB_S2_R2.fa mapping/R2_5prime/NB_S2_R2_first10.fa

## distribution of read lengths
python scripts/seqlength_distribution.py mapping/R1_5prime/C1_S1_R1.fa > mapping/R1_5prime/C1_seqlengths.txt
python scripts/seqlength_distribution.py mapping/R1_5prime/NB_S2_R1.fa > mapping/R1_5prime/NB_seqlengths.txt
python scripts/seqlength_distribution.py mapping/R2_5prime/C1_S1_R2.fa > mapping/R2_5prime/C1_seqlengths.txt
python scripts/seqlength_distribution.py mapping/R2_5prime/NB_S2_R2.fa > mapping/R2_5prime/NB_seqlengths.txt
R --no-save < scripts/plot_read_length_distribution.R

weblogo -f mapping/R1_5prime/C1_S1_R1_first10.fa -D fasta -o mapping/R1_5prime/C1_S1_R1_first10.pdf -F pdf -A dna --errorbars no -c classic
weblogo -f mapping/R1_5prime/NB_S2_R1_first10.fa -D fasta -o mapping/R1_5prime/NB_S2_R1_first10.pdf -F pdf -A dna --errorbars no -c classic
weblogo -f mapping/R2_5prime/C1_S1_R2_last10.fa -D fasta -o mapping/R2_5prime/C1_S1_R2_last10.pdf -F pdf -A dna --errorbars no -c classic
weblogo -f mapping/R2_5prime/NB_S2_R2_last10.fa -D fasta -o mapping/R2_5prime/NB_S2_R2_last10.pdf -F pdf -A dna --errorbars no -c classic

cutadapt -g ^AGACTC -e 0 --overlap 6 --minimum-length 1 --discard-untrimmed -o barcodes/B_SC-BLESS_C1_S1_L001_R1_BCHLT_barcodes.fastq.gz -p barcodes/B_SC-BLESS_C1_S1_L001_R2_BCHLT_barcodes.fastq.gz trim_adapters/B_SC-BLESS_C1_S1_L001_R1_BCHLT_cutadapt_native_adapters_4.fastq.gz trim_adapters/B_SC-BLESS_C1_S1_L001_R2_BCHLT_cutadapt_native_adapters_4.fastq.gz > analysis/cutadapt_4.out
cutadapt -g ^AGACTC -e 0 --overlap 6 --minimum-length 1 --discard-untrimmed -o barcodes/B_SC-BLESS_NB_S2_L001_R1_BCHLT_barcodes.fastq.gz -p barcodes/B_SC-BLESS_NB_S2_L001_R2_BCHLT_barcodes.fastq.gz trim_adapters/B_SC-BLESS_NB_S2_L001_R1_BCHLT_cutadapt_native_adapters_4.fastq.gz trim_adapters/B_SC-BLESS_NB_S2_L001_R2_BCHLT_cutadapt_native_adapters_4.fastq.gz>> analysis/cutadapt_4.out
cutadapt -a GAGTCT -e 0.2 --overlap 5 --minimum-length 1 --discard-untrimmed -o barcodes/B_SC-BLESS_C1_S1_L001_R2_BCHLT_barcodes_2.fastq.gz -p barcodes/B_SC-BLESS_C1_S1_L001_R1_BCHLT_barcodes_2.fastq.gz barcodes/B_SC-BLESS_C1_S1_L001_R2_BCHLT_barcodes.fastq.gz barcodes/B_SC-BLESS_C1_S1_L001_R1_BCHLT_barcodes.fastq.gz >> analysis/cutadapt_4.out
cutadapt -a GAGTCT -e 0.2 --overlap 5 --minimum-length 1 --discard-untrimmed -o barcodes/B_SC-BLESS_NB_S2_L001_R2_BCHLT_barcodes_2.fastq.gz -p barcodes/B_SC-BLESS_NB_S2_L001_R1_BCHLT_barcodes_2.fastq.gz barcodes/B_SC-BLESS_NB_S2_L001_R2_BCHLT_barcodes.fastq.gz barcodes/B_SC-BLESS_NB_S2_L001_R1_BCHLT_barcodes.fastq.gz >> analysis/cutadapt_4.out

