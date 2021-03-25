library("DESeq2")
library("ggplot2")
library("EnsDb.Hsapiens.v75")
library("AnnotationDbi")
library("ComplexHeatmap")
library("circlize")

RawCount <- readRDS("data/RNA-seq_RawCounts.Rds")

conditions <- factor(c(rep("Stroma",10),rep("Normal",10),rep("Tumor",10)))
case <- factor(c(rep(c(rep("squamous",4),rep("adeno",6)),3)))
sp.table <- data.frame(sampleName=colnames(RawCount), condition=conditions, case=case)

RNA.data <- DESeqDataSetFromMatrix(countData=RawCount, colData=sp.table, design=~condition)

#filtering out 0 or low count genes
#removing genes with less than  30 counts in 30 samples or less than 11 counts in all samples
keep <- rowSums(counts(RNA.data) <=10) < 30
keep.index <- which(keep == TRUE)
ddsRNA <- RNA.data[keep.index,]

#Running DESeq2 
ddsRNA <- estimateSizeFactors(ddsRNA)

rna.all <- DESeq(ddsRNA)
#Getting DEGs results between tumor versus normal
res.rna.all <- results(rna.all, contrast = c("condition","Tumor","Normal"))

#Adding metadata: gene symbol, entrez ID, and gene biotype
res.rna.all$symbol <- mapIds(EnsDb.Hsapiens.v75,keys = rownames(res.rna.all),column = "SYMBOL", keytype = "GENEID", multiVals = "first")
res.rna.all$entrezid <- mapIds(EnsDb.Hsapiens.v75,keys = rownames(res.rna.all),column = "ENTREZID", keytype = "GENEID", multiVals = "first")
res.rna.all$genetype <- mapIds(EnsDb.Hsapiens.v75,keys = rownames(res.rna.all),column = "GENEBIOTYPE", keytype = "GENEID", multiVals = "first")
#replacing missing gene symbols with rownames
na.index <- which(is.na(res.rna.all$symbol) == T)
res.rna.all$symbol[na.index] = rownames(res.rna.all)[na.index]
#Getting statistically significant DEGs
SigDEG.all <- subset(res.rna.all, abs(res.rna.all$log2FoldChange) >= 1 & res.rna.all$log2FoldChange != "Inf" & 
                       res.rna.all$log2FoldChange != "-Inf"  & res.rna.all$padj <= 0.01)
norm.count.RNA <- counts(rna.all, normalized = T)
sig.nor.count.RNA <- subset(norm.count.RNA, row.names(norm.count.RNA) %in% rownames(SigDEG.all))
#4,118 DEGs were found to be significantly differential expressed between tumor versus normal tissues.
