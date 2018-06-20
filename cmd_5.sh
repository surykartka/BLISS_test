cutadapt -g ^AGACTC -e 0 --overlap 6 --minimum-length 30 --discard-untrimmed -o barcodes/B_SC-BLESS_C1_S1_L001_R1_BCHLT_barcodes_30.fastq.gz -p barcodes/B_SC-BLESS_C1_S1_L001_R2_BCHLT_barcodes_30.fastq.gz trim_adapters/B_SC-BLESS_C1_S1_L001_R1_BCHLT_cutadapt_native_adapters_4.fastq.gz trim_adapters/B_SC-BLESS_C1_S1_L001_R2_BCHLT_cutadapt_native_adapters_4.fastq.gz > analysis/cutadapt_5.out
cutadapt -g ^AGACTC -e 0 --overlap 6 --minimum-length 30 --discard-untrimmed -o barcodes/B_SC-BLESS_NB_S2_L001_R1_BCHLT_barcodes_30.fastq.gz -p barcodes/B_SC-BLESS_NB_S2_L001_R2_BCHLT_barcodes_30.fastq.gz trim_adapters/B_SC-BLESS_NB_S2_L001_R1_BCHLT_cutadapt_native_adapters_4.fastq.gz trim_adapters/B_SC-BLESS_NB_S2_L001_R2_BCHLT_cutadapt_native_adapters_4.fastq.gz>> analysis/cutadapt_5.out
cutadapt -a GAGTCT -e 0.2 --overlap 5 --minimum-length 1 --discard-untrimmed -o barcodes/B_SC-BLESS_C1_S1_L001_R2_BCHLT_barcodes_2_30.fastq.gz -p barcodes/B_SC-BLESS_C1_S1_L001_R1_BCHLT_barcodes_2_30.fastq.gz barcodes/B_SC-BLESS_C1_S1_L001_R2_BCHLT_barcodes_30.fastq.gz barcodes/B_SC-BLESS_C1_S1_L001_R1_BCHLT_barcodes_30.fastq.gz >> analysis/cutadapt_5.out
cutadapt -a GAGTCT -e 0.2 --overlap 5 --minimum-length 1 --discard-untrimmed -o barcodes/B_SC-BLESS_NB_S2_L001_R2_BCHLT_barcodes_2_30.fastq.gz -p barcodes/B_SC-BLESS_NB_S2_L001_R1_BCHLT_barcodes_2_30.fastq.gz barcodes/B_SC-BLESS_NB_S2_L001_R2_BCHLT_barcodes_30.fastq.gz barcodes/B_SC-BLESS_NB_S2_L001_R1_BCHLT_barcodes_30.fastq.gz >> analysis/cutadapt_5.out

cutadapt -g ^AGACTC -e 0 --overlap 6 --minimum-length 20 --discard-untrimmed -o barcodes/B_SC-BLESS_C1_S1_L001_R1_BCHLT_barcodes_20.fastq.gz -p barcodes/B_SC-BLESS_C1_S1_L001_R2_BCHLT_barcodes_20.fastq.gz trim_adapters/B_SC-BLESS_C1_S1_L001_R1_BCHLT_cutadapt_native_adapters_4.fastq.gz trim_adapters/B_SC-BLESS_C1_S1_L001_R2_BCHLT_cutadapt_native_adapters_4.fastq.gz > analysis/cutadapt_6.out
cutadapt -g ^AGACTC -e 0 --overlap 6 --minimum-length 20 --discard-untrimmed -o barcodes/B_SC-BLESS_NB_S2_L001_R1_BCHLT_barcodes_20.fastq.gz -p barcodes/B_SC-BLESS_NB_S2_L001_R2_BCHLT_barcodes_20.fastq.gz trim_adapters/B_SC-BLESS_NB_S2_L001_R1_BCHLT_cutadapt_native_adapters_4.fastq.gz trim_adapters/B_SC-BLESS_NB_S2_L001_R2_BCHLT_cutadapt_native_adapters_4.fastq.gz>> analysis/cutadapt_6.out
cutadapt -a GAGTCT -e 0.2 --overlap 5 --minimum-length 1 --discard-untrimmed -o barcodes/B_SC-BLESS_C1_S1_L001_R2_BCHLT_barcodes_2_20.fastq.gz -p barcodes/B_SC-BLESS_C1_S1_L001_R1_BCHLT_barcodes_2_20.fastq.gz barcodes/B_SC-BLESS_C1_S1_L001_R2_BCHLT_barcodes_30.fastq.gz barcodes/B_SC-BLESS_C1_S1_L001_R1_BCHLT_barcodes_30.fastq.gz >> analysis/cutadapt_6.out
cutadapt -a GAGTCT -e 0.2 --overlap 5 --minimum-length 1 --discard-untrimmed -o barcodes/B_SC-BLESS_NB_S2_L001_R2_BCHLT_barcodes_2_20.fastq.gz -p barcodes/B_SC-BLESS_NB_S2_L001_R1_BCHLT_barcodes_2_20.fastq.gz barcodes/B_SC-BLESS_NB_S2_L001_R2_BCHLT_barcodes_30.fastq.gz barcodes/B_SC-BLESS_NB_S2_L001_R1_BCHLT_barcodes_30.fastq.gz >> analysis/cutadapt_6.out

bowtie2 --fr -x hg_mapping/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 barcodes/B_SC-BLESS_C1_S1_L001_R1_BCHLT_barcodes_2_30.fastq.gz -2 barcodes/B_SC-BLESS_C1_S1_L001_R2_BCHLT_barcodes_2_30.fastq.gz --threads 6 -S hg_mapping/B_SC-BLESS_C1_S1_L001_BCHLT_bowtie.sam 2> analysis/bowtie.out
bowtie2 --fr -x hg_mapping/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 barcodes/B_SC-BLESS_NB_S2_L001_R1_BCHLT_barcodes_2_30.fastq.gz -2 barcodes/B_SC-BLESS_NB_S2_L001_R2_BCHLT_barcodes_2_30.fastq.gz --threads 6 -S hg_mapping/B_SC-BLESS_NB_S2_L001_BCHLT_bowtie.sam 2>> analysis/bowtie.out
bowtie2 --fr -N 1 -x hg_mapping/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 barcodes/B_SC-BLESS_C1_S1_L001_R1_BCHLT_barcodes_2_30.fastq.gz -2 barcodes/B_SC-BLESS_C1_S1_L001_R2_BCHLT_barcodes_2_30.fastq.gz --threads 6 -S hg_mapping/B_SC-BLESS_C1_S1_L001_BCHLT_bowtie_N1.sam 2>> analysis/bowtie.out
bowtie2 --fr -N 1 -x hg_mapping/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 barcodes/B_SC-BLESS_NB_S2_L001_R1_BCHLT_barcodes_2_30.fastq.gz -2 barcodes/B_SC-BLESS_NB_S2_L001_R2_BCHLT_barcodes_2_30.fastq.gz --threads 6 -S hg_mapping/B_SC-BLESS_NB_S2_L001_BCHLT_bowtie_N1.sam 2>> analysis/bowtie.out

## sort for igv
samtools sort -o hg_mapping/B_SC-BLESS_C1_S1_L001_BCHLT_bowtie_sorted.sam hg_mapping/B_SC-BLESS_C1_S1_L001_BCHLT_bowtie.sam
samtools sort -o hg_mapping/B_SC-BLESS_NB_S2_L001_BCHLT_bowtie_sorted.sam hg_mapping/B_SC-BLESS_NB_S2_L001_BCHLT_bowtie.sam
samtools sort -o hg_mapping/B_SC-BLESS_C1_S1_L001_BCHLT_bowtie_v1_m1_sorted.sam hg_mapping/B_SC-BLESS_C1_S1_L001_BCHLT_bowtie_v1_m1.sam
samtools sort -o hg_mapping/B_SC-BLESS_NB_S2_L001_BCHLT_bowtie_v1_m1_sorted.sam hg_mapping/B_SC-BLESS_NB_S2_L001_BCHLT_bowtie_v1_m1.sam

## convert to bam
samtools view -Sb hg_mapping/B_SC-BLESS_C1_S1_L001_BCHLT_bowtie.sam > hg_mapping/B_SC-BLESS_C1_S1_L001_BCHLT_bowtie.bam
samtools view -Sb hg_mapping/B_SC-BLESS_NB_S2_L001_BCHLT_bowtie.sam > hg_mapping/B_SC-BLESS_NB_S2_L001_BCHLT_bowtie.bam

## sort bam files
samtools sort hg_mapping/B_SC-BLESS_C1_S1_L001_BCHLT_bowtie.bam > hg_mapping/B_SC-BLESS_C1_S1_L001_BCHLT_bowtie_sorted.bam
samtools sort hg_mapping/B_SC-BLESS_NB_S2_L001_BCHLT_bowtie.bam > hg_mapping/B_SC-BLESS_NB_S2_L001_BCHLT_bowtie_sorted.bam

## and to bed
bedtools bamtobed -i hg_mapping/B_SC-BLESS_C1_S1_L001_BCHLT_bowtie.bam > hg_mapping/B_SC-BLESS_C1_S1_L001_BCHLT_bowtie.bed
bedtools bamtobed -i hg_mapping/B_SC-BLESS_NB_S2_L001_BCHLT_bowtie.bam > hg_mapping/B_SC-BLESS_NB_S2_L001_BCHLT_bowtie.bed

## sort bed
bedtools sort -i hg_mapping/B_SC-BLESS_C1_S1_L001_BCHLT_bowtie.bed > hg_mapping/B_SC-BLESS_C1_S1_L001_BCHLT_bowtie_sorted.bed
bedtools sort -i hg_mapping/B_SC-BLESS_NB_S2_L001_BCHLT_bowtie.bed > hg_mapping/B_SC-BLESS_NB_S2_L001_BCHLT_bowtie_sorted.bed

## create bedgraphs
bedtools genomecov -ibam hg_mapping/B_SC-BLESS_C1_S1_L001_BCHLT_bowtie_sorted.bam -bg > hg_mapping/B_SC-BLESS_C1_S1_L001_BCHLT_bowtie.bedgraph
bedtools genomecov -ibam hg_mapping/B_SC-BLESS_NB_S2_L001_BCHLT_bowtie_sorted.bam -bg > hg_mapping/B_SC-BLESS_NB_S2_L001_BCHLT_bowtie.bedgraph
