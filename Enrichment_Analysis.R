### load libraries
library(DESeq2)
library(Rsubread)
library(ggplot2)
library(annotate)
#library("pheatmap")
library("AnnotationDbi")
library('org.Hs.eg.db')
library('org.Mm.eg.db')



### path
setwd('your/directory')

#Run FeatureCount on bam files
fc <- featureCounts(files=c("Input_1.bam",'Input_2.bam','Input_3.bam', "Target_1.bam",'Target_2.bam','Target_3.bam'), isPairedEnd=TRUE, nthreads=10, annot.inbuilt="your reference genome", minMQS=20, countChimericFragments=FALSE)

#Read in the meta data table
sampleTable <- read.csv("meta.csv",row.names=1,header=T, stringsAsFactors=F)
sampleTable

#An example meta data table looks like this
#               NAME  BATCH      GROUP
#Input_1.bam  Input_1    B1 Whole_cell
#Input_2.bam  Input_2    B2 Whole_cell
#Input_3.bam  Input_3    B3 Whole_cell
#Target_1.bam Target_1   B1        ROI
#Target_2.bam Target_2   B2        R0I
#Target_3.bam Target_3   B3        ROI

#Filter genes
countdata <- fc$counts[ rowSums(fc$counts) >= 2, ]
head(countdata)

coldata <- sampleTable
head(coldata)

#Generate DESeq2 object
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ BATCH + GROUP)

#Variance stabilization
vsd <- vst(ddsMat)

#PCA plots
plotPCA(vsd, "GROUP")
plotPCA(vsd, "BATCH")

#Run DESeq2
dds <- DESeq(ddsMat)

#Output normalized gene count data
dds <- estimateSizeFactors(dds) 
cnt_table <- counts(dds, normalized=TRUE)
write.csv(cnt_table, "Normalized_gene_counts.csv")

#Calculate enrichment of genes in ROI over the whole cell
res <- results(dds, contrast=c("GROUP","ROI","Whole_cell"))

summary(res)
mcols(res, use.names=TRUE)

#Obtain the gene symbols
res$symbol <- getSYMBOL(row.names(res), data='org.Hs.eg') #or 'org.Mm.eg'

#Output data
write.csv(res, "Fold_enrichment_data.csv")

#Save data
saveRDS(res, file = "res.rds")
save.image(file = "DEseq.RData")

