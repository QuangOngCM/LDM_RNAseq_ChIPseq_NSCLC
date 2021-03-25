library("DESeq2")
library("ggplot2")

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

#Running DESeq2 log transformation
ddsRNA <- estimateSizeFactors(ddsRNA)
rld <- rlogTransformation(ddsRNA, blind = FALSE)

#Making plot
#svg("Figure1c.svg", width = 8, height = 7)

pcaData <- plotPCA(rld, intgroup=c("condition", "case"),  returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=case)) +
  scale_color_manual(values = c("#3288bd", "#1a9850", "#d53e4f")) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 20, face = "bold"),
        legend.title=element_text(size=13), 
        legend.text=element_text(size=15)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

#dev.off() 