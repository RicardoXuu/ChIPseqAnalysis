rm(list = ls())


library("ChIPseeker")
library("org.Hs.eg.db")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("clusterProfiler")

#为了方便调用，将TxDb.Hsapiens.UCSC.hg19.knownGene命名为txdb
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

#下面导入之前masc2做call peak 得到的.narrowPeak文件
rmdup_peaks <- readPeakFile("CP_rmdup_peaks.narrowPeak")
head(rmdup_peaks)
keepChr = !grepl('_',seqlevels(rmdup_peaks))
table(keepChr)
seqlevels(rmdup_peaks,pruning.mode="coarse") <- seqlevels(rmdup_peaks)[keepChr]


covplot(rmdup_peaks,weightCol=5)

#从数据库引入promoter并将peaks数据比对到promoter数据上
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
head(promoter)
tagMatrix <- getTagMatrix(rmdup_peaks, windows=promoter)
head(tagMatrix)

#下面两个画图函数，因为tagMatrix文件过大，非常耗费内存，不要轻易跑
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            conf=0.95,resample = 1000,
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

#利用annotatePeak函数进行peaks注释
peakAnno <- annotatePeak(rmdup_peaks, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

 #这个数据集不能用head函数查看，原因不了解,下面尝试格式转换并保存成.csv文件
library(dplyr)      
as.GRanges(peakAnno) %>% head(5)
save_peakAnno <- as.data.frame(peakAnno)
head(save_peakAnno)
write.table(save_peakAnno,file = "save_peakAnno.csv")

#对得到的注释信息做一些可视化
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
library("ggupset")
upsetplot(peakAnno)
plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")


