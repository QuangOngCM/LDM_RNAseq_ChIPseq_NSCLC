library("DESeq2")
library("ggplot2")
library("EnsDb.Hsapiens.v75")
library("AnnotationDbi")
library("ggrepel")

source("RNA-seq_DEGsAnalysis.R")

#Making subset data by NSCLC subtypes: lung adenocarcinoma (LUAD) and lung squamous cell carcinoma (LUSC)

rna.LUAD <- ddsRNA[,colData(ddsRNA)$case %in% c("adeno") & colData(ddsRNA)$condition %in% c("Normal","Tumor")]
rna.LUSC <- ddsRNA[,colData(ddsRNA)$case %in% c("squamous") & colData(ddsRNA)$condition %in% c("Normal","Tumor")]
rna.LUAD$condition <- droplevels(rna.LUAD$condition)
rna.LUSC$condition <- droplevels(rna.LUSC$condition)

#Filtering out low count genes
##LUAD
keep <- rowSums(counts(rna.LUAD) <=30) < 12
keep.index <- which(keep == TRUE)

rna.LUAD.filtered <- rna.LUAD[keep.index,]
rna.LUAD.filtered <- DESeq(rna.LUAD.filtered)
count.LUAD <- counts(rna.LUAD.filtered, normalized = F)
norm.count.LUAD <- counts(rna.LUAD.filtered, normalized = T)
##LUSC
keep <- rowSums(counts(rna.LUSC) <=30) < 8
keep.index <- which(keep == TRUE)

rna.LUSC.filtered <- rna.LUSC[keep.index,]
rna.LUSC.filtered <- DESeq(rna.LUSC.filtered)
count.LUSC <- counts(rna.LUSC.filtered, normalized = F)
norm.count.LUSC <- counts(rna.LUSC.filtered, normalized = T)

##Running DESeq2
res.rna.LUAD <- results(rna.LUAD.filtered, contrast = c("condition","Tumor","Normal"))
res.rna.LUSC <- results(rna.LUSC.filtered, contrast = c("condition","Tumor","Normal"))

res.rna.LUAD$symbol <- mapIds(EnsDb.Hsapiens.v75,keys = rownames(res.rna.LUAD),column = "SYMBOL", keytype = "GENEID", multiVals = "first")
res.rna.LUAD$entrezid <- mapIds(EnsDb.Hsapiens.v75,keys = rownames(res.rna.LUAD),column = "ENTREZID", keytype = "GENEID", multiVals = "first")
res.rna.LUAD$genetype <- mapIds(EnsDb.Hsapiens.v75,keys = rownames(res.rna.LUAD),column = "GENEBIOTYPE", keytype = "GENEID", multiVals = "first")

res.rna.LUSC$symbol <- mapIds(EnsDb.Hsapiens.v75,keys = rownames(res.rna.LUSC),column = "SYMBOL", keytype = "GENEID", multiVals = "first")
res.rna.LUSC$entrezid <- mapIds(EnsDb.Hsapiens.v75,keys = rownames(res.rna.LUSC),column = "ENTREZID", keytype = "GENEID", multiVals = "first")
res.rna.LUSC$genetype <- mapIds(EnsDb.Hsapiens.v75,keys = rownames(res.rna.LUSC),column = "GENEBIOTYPE", keytype = "GENEID", multiVals = "first")

na.index <- which(is.na(res.rna.LUAD$symbol) == TRUE)
res.rna.LUAD$symbol[na.index] = rownames(res.rna.LUAD)[na.index]

na.index <- which(is.na(res.rna.LUSC$symbol) == TRUE)
res.rna.LUSC$symbol[na.index] = rownames(res.rna.LUSC)[na.index]

SigDEG.LUAD <- subset(res.rna.LUAD, abs(res.rna.LUAD$log2FoldChange) >= 1 & res.rna.LUAD$log2FoldChange != "Inf" & 
                        res.rna.LUAD$log2FoldChange != "-Inf"  & res.rna.LUAD$padj <= 0.01)
SigDEG.LUSC <- subset(res.rna.LUSC, abs(res.rna.LUSC$log2FoldChange) >= 1 & res.rna.LUSC$log2FoldChange != "Inf" & 
                        res.rna.LUSC$log2FoldChange != "-Inf"  & res.rna.LUSC$padj <= 0.01)

#Selecting significant up/down regulated genes in LUAD and LUSC
Sig.Up.Res.LUAD <- SigDEG.LUAD[SigDEG.LUAD$log2FoldChange > 0,]
up.genes.LUAD <- data.frame(rownames(Sig.Up.Res.LUAD), Sig.Up.Res.LUAD$symbol)
Sig.Down.Res.LUAD <- SigDEG.LUAD[SigDEG.LUAD$log2FoldChange < 0,]
down.genes.LUAD <- data.frame(rownames(Sig.Down.Res.LUAD), Sig.Down.Res.LUAD$symbol)

Sig.Up.Res.LUSC <- SigDEG.LUSC[SigDEG.LUSC$log2FoldChange > 0,]
up.genes.LUSC <- data.frame(rownames(Sig.Up.Res.LUSC), Sig.Up.Res.LUSC$symbol)
Sig.Down.Res.LUSC <- SigDEG.LUSC[SigDEG.LUSC$log2FoldChange < 0,]
down.genes.LUSC <- data.frame(rownames(Sig.Down.Res.LUSC), Sig.Down.Res.LUSC$symbol)

#making figure 2b plot
SigDEG.LUAD.df <- data.frame(res.rna.LUAD)
SigDEG.LUAD.df$thershold <- ifelse(SigDEG.LUAD.df$log2FoldChange >= 1 & SigDEG.LUAD.df$padj <= 0.05, "gained", 
                                   ifelse(SigDEG.LUAD.df$log2FoldChange <= -1 & SigDEG.LUAD.df$padj <= 0.05, "loss","unchange"))
cols <- c("gained" = "red", "loss" = "blue", "unchange" = "darkgrey")

Sig.Up.Res.LUAD <- subset(SigDEG.LUAD.df, SigDEG.LUAD.df$log2FoldChange > 0 & SigDEG.LUAD.df$padj != "NA")
Sig.Up.Res.LUAD_ordered <- Sig.Up.Res.LUAD[order(Sig.Up.Res.LUAD$log2FoldChange, decreasing = T), ] 
Sig.Up.Res.LUAD_ordered <- Sig.Up.Res.LUAD_ordered[Sig.Up.Res.LUAD_ordered$genetype != c("pseudogene","IG_V_gene"),]


Sig.Up.Res.LUAD_ordered$genelabels <- ""
Sig.Up.Res.LUAD_ordered$genelabels[1:15] <- Sig.Up.Res.LUAD_ordered$symbol[1:15]

Sig.Down.Res.LUAD <- subset(SigDEG.LUAD.df, SigDEG.LUAD.df$log2FoldChange < 0 & SigDEG.LUAD.df$padj != "NA")
Sig.Down.Res.LUAD_ordered <- Sig.Down.Res.LUAD[order(Sig.Down.Res.LUAD$log2FoldChange, decreasing = F), ] 
Sig.Down.Res.LUAD_ordered <- Sig.Down.Res.LUAD_ordered[Sig.Down.Res.LUAD_ordered$genetype != c("pseudogene","IG_V_gene"),]

Sig.Down.Res.LUAD_ordered$genelabels <- ""
Sig.Down.Res.LUAD_ordered$genelabels[1:15] <- Sig.Down.Res.LUAD_ordered$symbol[1:15]

res.rna.all.df.LUAD <- as.data.frame(rbind(Sig.Up.Res.LUAD_ordered,Sig.Down.Res.LUAD_ordered))

plot1 <- ggplot(res.rna.all.df.LUAD, aes(x=log2FoldChange, y=-log10(padj), label = genelabels)) +
  geom_text_repel() +
  ggtitle(label = "RNA-seq DEGs", subtitle = "LUAD vs Normal") +
  geom_point(aes(color=factor(thershold)), alpha=1/2, size=0.8) + xlim(-7.5, 7.5) +
  scale_colour_manual(values = cols) +
  xlab("log2 fold change") + ylab("-log10 adjusted P")+
  theme_bw(base_size = 14)+
  theme(legend.position = "none") +
  geom_hline(yintercept = 1.30103, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = 1, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = -1, colour="#990000", linetype="dashed") 

#svg("Figure2b.svg", width = 9, height = 6)
plot1
#dev.off() 

#making figure 2c plot

SigDEG.LUSC.df <- data.frame(res.rna.LUSC)
SigDEG.LUSC.df$thershold <- ifelse(SigDEG.LUSC.df$log2FoldChange >= 1 & SigDEG.LUSC.df$padj <= 0.05, "gained", 
                                   ifelse(SigDEG.LUSC.df$log2FoldChange <= -1 & SigDEG.LUSC.df$padj <= 0.05, "loss","unchange"))
cols <- c("gained" = "red", "loss" = "blue", "unchange" = "darkgrey")

Sig.Up.Res.LUSC <- subset(SigDEG.LUSC.df, SigDEG.LUSC.df$log2FoldChange > 0 & SigDEG.LUSC.df$padj != "NA")
Sig.Up.Res.LUSC_ordered <- Sig.Up.Res.LUSC[order(Sig.Up.Res.LUSC$log2FoldChange, decreasing = T), ] 
Sig.Up.Res.LUSC_ordered <- Sig.Up.Res.LUSC_ordered[Sig.Up.Res.LUSC_ordered$genetype != c("pseudogene","IG_V_gene"),]


Sig.Up.Res.LUSC_ordered$genelabels <- ""
Sig.Up.Res.LUSC_ordered$genelabels[1:15] <- Sig.Up.Res.LUSC_ordered$symbol[1:15]

Sig.Down.Res.LUSC <- subset(SigDEG.LUSC.df, SigDEG.LUSC.df$log2FoldChange < 0 & SigDEG.LUSC.df$padj != "NA")
Sig.Down.Res.LUSC_ordered <- Sig.Down.Res.LUSC[order(Sig.Down.Res.LUSC$log2FoldChange, decreasing = F), ] 
Sig.Down.Res.LUSC_ordered <- Sig.Down.Res.LUSC_ordered[Sig.Down.Res.LUSC_ordered$genetype != c("pseudogene","IG_V_gene"),]

Sig.Down.Res.LUSC_ordered$genelabels <- ""
Sig.Down.Res.LUSC_ordered$genelabels[1:15] <- Sig.Down.Res.LUSC_ordered$symbol[1:15]

res.rna.all.df.LUSC <- as.data.frame(rbind(Sig.Up.Res.LUSC_ordered,Sig.Down.Res.LUSC_ordered))

plot2 <- ggplot(res.rna.all.df.LUSC, aes(x=log2FoldChange, y=-log10(padj), label = genelabels)) +
  geom_text_repel() +
  ggtitle(label = "RNA-seq DEGs", subtitle = "LUSC vs Normal") +
  geom_point(aes(color=factor(thershold)), alpha=1/2, size=0.8) +
  scale_colour_manual(values = cols) +
  xlab("log2 fold change") + ylab("-log10 adjusted P")+
  theme_bw(base_size = 14)+
  theme(legend.position = "none") +
  geom_hline(yintercept = 1.30103, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = 1, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = -1, colour="#990000", linetype="dashed") 

#svg("Figure2c.svg", width = 9, height = 6)
plot2
#dev.off() 