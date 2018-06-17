cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG --minimum-length 6 -o trim_adapters/B_SC-BLESS_C1_S1_L001_R1_BCHLT_cutadapt_1.fastq.gz fastq/B_SC-BLESS_C1_S1_L001_R1_BCHLT.fastq.gz > analysis/cutadapt.out
cutadapt -g ^AGACTCT -e 0 --overlap 7 --minimum-length 1 -o trim_adapters/B_SC-BLESS_C1_S1_L001_R1_BCHLT_cutadapt_2.fastq.gz --discard-untrimmed trim_adapters/B_SC-BLESS_C1_S1_L001_R1_BCHLT_cutadapt_1.fastq.gz >> analysis/cutadapt.out

cutadapt -a TCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT --minimum-length 6 -o trim_adapters/B_SC-BLESS_C1_S1_L001_R2_BCHLT_cutadapt_1.fastq.gz fastq/B_SC-BLESS_C1_S1_L001_R2_BCHLT.fastq.gz >> analysis/cutadapt.out
cutadapt -g ^AGACTCT -e 0 --overlap 7 --minimum-length 1 -o trim_adapters/B_SC-BLESS_C1_S1_L001_R2_BCHLT_cutadapt_2.fastq.gz --discard-untrimmed trim_adapters/B_SC-BLESS_C1_S1_L001_R2_BCHLT_cutadapt_1.fastq.gz >> analysis/cutadapt.out

cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG --minimum-length 6 -o trim_adapters/B_SC-BLESS_NB_S2_L001_R1_BCHLT_cutadapt_1.fastq.gz fastq/B_SC-BLESS_NB_S2_L001_R1_BCHLT.fastq.gz >> analysis/cutadapt.out
cutadapt -g ^AGACTCT -e 0 --overlap 7 --minimum-length 1 -o trim_adapters/B_SC-BLESS_NB_S2_L001_R1_BCHLT_cutadapt_2.fastq.gz --discard-untrimmed trim_adapters/B_SC-BLESS_NB_S2_L001_R1_BCHLT_cutadapt_1.fastq.gz >> analysis/cutadapt.out

cutadapt -a TCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT --minimum-length 6 -o trim_adapters/B_SC-BLESS_NB_S2_L001_R2_BCHLT_cutadapt_1.fastq.gz fastq/B_SC-BLESS_NB_S2_L001_R2_BCHLT.fastq.gz >> analysis/cutadapt.out
cutadapt -g ^AGACTCT -e 0 --overlap 7 --minimum-length 1 -o trim_adapters/B_SC-BLESS_NB_S2_L001_R2_BCHLT_cutadapt_2.fastq.gz --discard-untrimmed trim_adapters/B_SC-BLESS_NB_S2_L001_R2_BCHLT_cutadapt_1.fastq.gz >> analysis/cutadapt.out
