cutadapt -g ^AGACTC -e 0 --overlap 6 --minimum-length 30 --discard-untrimmed -o barcodes/B_SC-BLESS_C1_S1_L001_R1_BCHLT_barcodes_30.fastq.gz -p barcodes/B_SC-BLESS_C1_S1_L001_R2_BCHLT_barcodes_30.fastq.gz trim_adapters/B_SC-BLESS_C1_S1_L001_R1_BCHLT_cutadapt_native_adapters_4.fastq.gz trim_adapters/B_SC-BLESS_C1_S1_L001_R2_BCHLT_cutadapt_native_adapters_4.fastq.gz > analysis/cutadapt_5.out
cutadapt -g ^AGACTC -e 0 --overlap 6 --minimum-length 30 --discard-untrimmed -o barcodes/B_SC-BLESS_NB_S2_L001_R1_BCHLT_barcodes_30.fastq.gz -p barcodes/B_SC-BLESS_NB_S2_L001_R2_BCHLT_barcodes_30.fastq.gz trim_adapters/B_SC-BLESS_NB_S2_L001_R1_BCHLT_cutadapt_native_adapters_4.fastq.gz trim_adapters/B_SC-BLESS_NB_S2_L001_R2_BCHLT_cutadapt_native_adapters_4.fastq.gz>> analysis/cutadapt_5.out

bowtie2 -x mapping/bowtie_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 barcodes/B_SC-BLESS_NB_S2_L001_R1_BCHLT_barcodes_30.fastq.gz -2 barcodes/B_SC-BLESS_NB_S2_L001_R2_BCHLT_barcodes_30.fastq.gz --threads 6 -S mapping/B_SC-BLESS_NB_S2_L001_BCHLT_bowtie.sam
bowtie2 -x mapping/bowtie_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 barcodes/B_SC-BLESS_C1_S1_L001_R1_BCHLT_barcodes_30.fastq.gz -2 barcodes/B_SC-BLESS_C1_S1_L001_R2_BCHLT_barcodes_30.fastq.gz --threads 6 -S mapping/B_SC-BLESS_C1_S1_L001_BCHLT_bowtie.sam