rm(list = ls())
# 值得注意的是有些包检查会有版本的问题，包括 ChIPseeker  

install.packages("devtools")
library(devtools) 
  
BiocManager::install("airway")
#这里还需要安装DEseq2,limma,edgeR,之前已安装过，所以没有列出

BiocManager::install(c('ChIPpeakAnno','ChIPseeker'))
library(ChIPpeakAnno) 
library(ChIPseeker) 

BiocManager::install('TxDb.Hsapiens.UCSC.hg19.knownGene',
                        ask=F,suppressUpdates=T)
BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene',
                        ask=F,suppressUpdates=T)
#在这里我们只做人的ChIPseq分析，所以只下载了人的数据库
library(TxDb.Hsapiens.UCSC.hg19.knownGene) 
library(TxDb.Hsapiens.UCSC.hg38.knownGene) 

