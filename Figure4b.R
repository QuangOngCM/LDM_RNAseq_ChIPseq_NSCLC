library("DESeq2")
library("ggplot2")
library("EnsDb.Hsapiens.v75")
library("AnnotationDbi")
library("circlize")
library("cowplot")

source("RNA-seq_DEGsAnalysis.R")

#Making plot
gene.data <- res.rna.all[res.rna.all$symbol == "EZH2",] 
gene_ID <- rownames(gene.data)
gene_symbol <- as.character(gene.data$symbol)
L2FC = round(gene.data$log2FoldChange, digits = 2)
FDR = formatC(gene.data$padj, format = "e", digits = 2)

gene_counts <- counts(rna.all[gene_ID,], normalized = TRUE)
m <- data.frame(counts = as.numeric(gene_counts), group = as.factor(rna.all$condition), sub.group = as.factor(rna.all$case))
#svg("Figure4bEZH2.svg", width = 6, height = 5)
ggplot(m, aes(x = group, y = counts)) +
  geom_boxplot(alpha = 0.80) +
  geom_jitter(aes(fill = sub.group), size = 3, alpha=0.5, shape = 21, position = position_jitterdodge()) +
  labs(fill = "Subtype")+
  theme(text = element_text(size = 13),
        axis.text = element_text(size = 15),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = c(0.2, 0.8))+
  labs(x = "Samples from This Study", y = "Normalized Counts ", title =paste0("Normalized Expression of ",gene_symbol,
                                                                         "\n L2FC = ",L2FC," (Tumor versus Normal)",
                                                                         "\n FDR = ",FDR))+
  scale_x_discrete(labels=c("Normal","Stroma","Tumor"))
#dev.off()
