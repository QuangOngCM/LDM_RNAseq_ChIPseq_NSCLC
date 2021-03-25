library("DESeq2")
library("ggplot2")
library("stringr")
library("EnsDb.Hsapiens.v86")
library("AnnotationDbi")
library("circlize")
library("cowplot")

TCGA.RawCount <- readRDS("data/TCGA101_RNA-seq_RawCounts.Rds")

#Adding data to DEseq2 format
conditions <- factor(c(rep("Normal",52),rep("Tumor",52),rep("Normal",49),rep("Tumor",49)))
case <- factor(c(rep("adeno",104),rep("squamous",98)))
sp.table <- data.frame(sampleName=colnames(TCGA.RawCount), condition=conditions, case=case)

TCGA.RNA.data <- DESeqDataSetFromMatrix(countData=TCGA.RawCount, colData=sp.table, design=~condition)

ddsTcgaRNA <- estimateSizeFactors(TCGA.RNA.data)
#Running Deseq2
rna.TCGA <- DESeq(ddsTcgaRNA)

res.rna.TCGA <- results(rna.TCGA, contrast = c("condition","Tumor","Normal"))
res.rna.TCGA$new_id <- str_replace(rownames(res.rna.TCGA),pattern = ".[0-9]+$",replacement = "")

res.rna.TCGA$symbol <- mapIds(EnsDb.Hsapiens.v86,keys = res.rna.TCGA$new_id,column = "SYMBOL", keytype = "GENEID", multiVals = "first")
res.rna.TCGA$entrezid <- mapIds(EnsDb.Hsapiens.v86,keys = res.rna.TCGA$new_id,column = "ENTREZID", keytype = "GENEID", multiVals = "first")
res.rna.TCGA$genetype <- mapIds(EnsDb.Hsapiens.v86,keys = res.rna.TCGA$new_id,column = "GENEBIOTYPE", keytype = "GENEID", multiVals = "first")

#Making plot
##P1
gene.data <- subset(res.rna.TCGA, res.rna.TCGA$symbol == "COL1A1") #COL1A1,COL1A2,IGLC3,IGHG3

gene_ID <- rownames(gene.data)
gene_symbol <- as.character(gene.data$symbol)
L2FC = format(gene.data$log2FoldChange, scientific = T)
FDR = format(gene.data$padj, scientific = T)

gene_counts <- counts(rna.TCGA[gene_ID,], normalized = TRUE)
m <- data.frame(counts = as.numeric(gene_counts), group = as.factor(rna.TCGA$condition), sub.group = as.factor(rna.TCGA$case))
p1 <- ggplot(m, aes(x = group, y = counts)) +
  geom_boxplot(alpha = 0.80) +
  geom_jitter(aes(fill = sub.group), size = 3, alpha=0.5, shape = 21, position = position_jitterdodge()) +
  theme(text = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 15),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "")+
  labs(x = "Experimental Group", y = "Normalized Counts ", title =paste0(gene_symbol))+
  scale_x_discrete(labels=c("Normal","Tumor"))
##P2
gene.data <- subset(res.rna.TCGA, res.rna.TCGA$symbol == "COL1A2") #COL1A1,COL1A2,IGLC3,IGHG3

gene_ID <- rownames(gene.data)
gene_symbol <- as.character(gene.data$symbol)
L2FC = format(gene.data$log2FoldChange, scientific = T)
FDR = format(gene.data$padj, scientific = T)

gene_counts <- counts(rna.TCGA[gene_ID,], normalized = TRUE)
m <- data.frame(counts = as.numeric(gene_counts), group = as.factor(rna.TCGA$condition), sub.group = as.factor(rna.TCGA$case))
p2 <- ggplot(m, aes(x = group, y = counts)) +
  geom_boxplot(alpha = 0.80) +
  geom_jitter(aes(fill = sub.group), size = 3, alpha=0.5, shape = 21, position = position_jitterdodge()) +
  theme(text = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 15),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "")+
  labs(x = "Experimental Group", y = "Normalized Counts ", title =paste0(gene_symbol))+
  scale_x_discrete(labels=c("Normal","Tumor"))
##P3
gene.data <- subset(res.rna.TCGA, res.rna.TCGA$symbol == "IGLC3") #COL1A1,COL1A2,IGLC3,IGHG3

gene_ID <- rownames(gene.data)
gene_symbol <- as.character(gene.data$symbol)
L2FC = format(gene.data$log2FoldChange, scientific = T)
FDR = format(gene.data$padj, scientific = T)

gene_counts <- counts(rna.TCGA[gene_ID,], normalized = TRUE)
m <- data.frame(counts = as.numeric(gene_counts), group = as.factor(rna.TCGA$condition), sub.group = as.factor(rna.TCGA$case))
p3 <- ggplot(m, aes(x = group, y = counts)) +
  geom_boxplot(alpha = 0.80) +
  geom_jitter(aes(fill = sub.group), size = 3, alpha=0.5, shape = 21, position = position_jitterdodge()) +
  theme(text = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 15),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "")+
  labs(x = "Experimental Group", y = "Normalized Counts ", title =paste0(gene_symbol))+
  scale_x_discrete(labels=c("Normal","Tumor"))
##P4
gene.data <- subset(res.rna.TCGA, res.rna.TCGA$symbol == "IGHG3") #COL1A1,COL1A2,IGLC3,IGHG3

gene_ID <- rownames(gene.data)
gene_symbol <- as.character(gene.data$symbol)
L2FC = format(gene.data$log2FoldChange, scientific = T)
FDR = format(gene.data$padj, scientific = T)

gene_counts <- counts(rna.TCGA[gene_ID,], normalized = TRUE)
m <- data.frame(counts = as.numeric(gene_counts), group = as.factor(rna.TCGA$condition), sub.group = as.factor(rna.TCGA$case))
p4 <- ggplot(m, aes(x = group, y = counts)) +
  geom_boxplot(alpha = 0.80) +
  geom_jitter(aes(fill = sub.group), size = 3, alpha=0.5, shape = 21, position = position_jitterdodge()) +
  theme(text = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 15),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "")+
  labs(x = "Experimental Group", y = "Normalized Counts ", title =paste0(gene_symbol))+
  scale_x_discrete(labels=c("Normal","Tumor"))
##Final plot
#svg("Figure2f.svg", width = 8, height = 5)
plot_grid(p1, p2, p3, p4, nrow = 2, ncol=2)
#dev.off()