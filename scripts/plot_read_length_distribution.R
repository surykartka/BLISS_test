setwd('~/Documents/BLISS/test_june/mapping/')

t <- read.delim('R1_5prime/C1_seqlengths.txt', header=F)
pdf('R1_C1_distr.pdf', 5, 4)
hist(t$V1, xlim=c(0,80), col='grey', xlab='Read length (barcodes left)', main='C1_S1 R1 lengths')
dev.off()

t <- read.delim('R1_5prime/NB_seqlengths.txt', header=F)
pdf('R1_NB_distr.pdf', 5, 4)
hist(t$V1, xlim=c(0,80), col='grey', xlab='Read length (barcodes left)', main='NB_S2 R1 lengths')
dev.off()

t <- read.delim('R2_5prime/C1_seqlengths.txt', header=F)
pdf('R2_C1_distr.pdf', 5, 4)
hist(t$V1, xlim=c(0,80), col='grey', xlab='Read length (barcodes left)', main='C1_S1 R2 lengths')
dev.off()

t <- read.delim('R2_5prime/NB_seqlengths.txt', header=F)
pdf('R2_NB_distr.pdf', 5, 4)
hist(t$V1, xlim=c(0,80), col='grey', xlab='Read length (barcodes left)', main='NB_S2 R2 lengths')
dev.off()