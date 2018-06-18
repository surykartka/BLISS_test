cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG --discard-trimmed --minimum-length 6 -o trim_adapters/B_SC-BLESS_C1_S1_L001_R1_BCHLT_cutadapt_paired_1.fastq.gz -p trim_adapters/B_SC-BLESS_C1_S1_L001_R1_BCHLT_cutadapt_paired_2.fastq.gz fastq/B_SC-BLESS_C1_S1_L001_R1_BCHLT.fastq.gz fastq/B_SC-BLESS_C1_S1_L001_R2_BCHLT.fastq.gz > analysis/cutadapt_2.out
cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG --discard-trimmed --minimum-length 6 -o trim_adapters/B_SC-BLESS_NB_S2_L001_R1_BCHLT_cutadapt_paired_1.fastq.gz -p trim_adapters/B_SC-BLESS_NB_S2_L001_R1_BCHLT_cutadapt_paired_2.fastq.gz fastq/B_SC-BLESS_NB_S2_L001_R1_BCHLT.fastq.gz fastq/B_SC-BLESS_NB_S2_L001_R2_BCHLT.fastq.gz >> analysis/cutadapt_2.out

fastq_to_fasta -Q33 -i mapping/B_SC-BLESS_C1_S1_L001_R1_BCHLT_cutadapt_paired_1.fastq -o mapping/R1_5prime/C1_S1_R1.fa # remember to extract fastq here!!
fastq_to_fasta -Q33 -i mapping/B_SC-BLESS_NB_S2_L001_R1_BCHLT_cutadapt_paired_1.fastq -o mapping/R1_5prime/NB_S2_R1.fa

python scripts/cut_10.py mapping/R1_5prime/C1_S1_R1.fa mapping/R1_5prime/C1_S1_R1_first10.fa
python scripts/cut_10.py mapping/R1_5prime/NB_S2_R1.fa mapping/R1_5prime/NB_S2_R1_first10.fa

weblogo -f mapping/R1_5prime/C1_S1_R1_first10.fa -D fasta -o mapping/R1_5prime/C1_S1_R1_first10.pdf -F pdf -A dna --errorbars no -c classic
weblogo -f mapping/R1_5prime/NB_S2_R1_first10.fa -D fasta -o mapping/R1_5prime/NB_S2_R1_first10.pdf -F pdf -A dna --errorbars no -c classic
