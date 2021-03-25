library("DESeq2")
library("ggplot2")
library("EnsDb.Hsapiens.v75")
library("AnnotationDbi")
library("circlize")
library("cowplot")

source("RNA-seq_DEGsAnalysis.R")
#Making plot
##P1
gene.data <- res.rna.all[res.rna.all$symbol == "COL1A1",] #COL1A1,COL1A2,IGLC3,IGHG3
gene_ID <- rownames(gene.data)
gene_symbol <- as.character(gene.data$symbol)

gene_counts <- counts(rna.all[gene_ID,], normalized = TRUE)
m <- data.frame(counts = as.numeric(gene_counts), group = as.factor(rna.all$condition), sub.group = as.factor(rna.all$case))
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
  scale_x_discrete(labels=c("Normal","Stroma","Tumor"))
##P2
gene.data <- res.rna.all[res.rna.all$symbol == "COL1A2",] #COL1A1,COL1A2,IGLC3,IGHG3
gene_ID <- rownames(gene.data)
gene_symbol <- as.character(gene.data$symbol)

gene_counts <- counts(rna.all[gene_ID,], normalized = TRUE)
m <- data.frame(counts = as.numeric(gene_counts), group = as.factor(rna.all$condition), sub.group = as.factor(rna.all$case))
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
  scale_x_discrete(labels=c("Normal","Stroma","Tumor"))
##P3
gene.data <- res.rna.all[res.rna.all$symbol == "IGLC3",] #COL1A1,COL1A2,IGLC3,IGHG3
gene_ID <- rownames(gene.data)
gene_symbol <- as.character(gene.data$symbol)

gene_counts <- counts(rna.all[gene_ID,], normalized = TRUE)
m <- data.frame(counts = as.numeric(gene_counts), group = as.factor(rna.all$condition), sub.group = as.factor(rna.all$case))
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
  scale_x_discrete(labels=c("Normal","Stroma","Tumor"))
##P4
gene.data <- res.rna.all[res.rna.all$symbol == "IGHG3",] #COL1A1,COL1A2,IGLC3,IGHG3
gene_ID <- rownames(gene.data)
gene_symbol <- as.character(gene.data$symbol)

gene_counts <- counts(rna.all[gene_ID,], normalized = TRUE)
m <- data.frame(counts = as.numeric(gene_counts), group = as.factor(rna.all$condition), sub.group = as.factor(rna.all$case))
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
  scale_x_discrete(labels=c("Normal","Stroma","Tumor"))
#Final plot
#svg("Figure2e.svg", width = 8, height = 5)
plot_grid(p1, p2, p3, p4, nrow = 2, ncol=2)
#dev.off()