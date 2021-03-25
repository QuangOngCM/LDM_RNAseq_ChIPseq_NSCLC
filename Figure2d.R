library("DESeq2")
library("ggplot2")
library("EnsDb.Hsapiens.v75")
library("AnnotationDbi")
library("ggrepel")
library("clusterProfiler")

source("Figure2b-c.R")

#Selecting significant up/down regulated genes in LUAD and LUSC
Sig.Up.Res.LUAD <- SigDEG.LUAD[SigDEG.LUAD$log2FoldChange > 0,]
Sig.Down.Res.LUAD <- SigDEG.LUAD[SigDEG.LUAD$log2FoldChange < 0,]

Sig.Up.Res.LUSC <- SigDEG.LUSC[SigDEG.LUSC$log2FoldChange > 0,]
Sig.Down.Res.LUSC <- SigDEG.LUSC[SigDEG.LUSC$log2FoldChange < 0,]
##

Up.LUAD <- as.character(Sig.Up.Res.LUAD$entrezid)
Down.LUAD <- as.character(Sig.Down.Res.LUAD$entrezid)
Up.LUSC <- as.character(Sig.Up.Res.LUSC$entrezid)
Down.LUSC <- as.character(Sig.Down.Res.LUSC$entrezid)

mylist <- list(Up.LUAD,Down.LUAD,Up.LUSC,Down.LUSC)
names(mylist) <- c("Up.LUAD", "Down.LUAD","Up.LUSC","Down.LUSC")

ck <- compareCluster(geneCluster = mylist, fun = "enrichKEGG", organism="hsa")

#svg("Figure2d.svg", width = 9.2, height = 7)
dotplot(ck, showCategory=10, includeAll = F) + ggtitle("KEGG pathway enrichment analysis")
#dev.off() 