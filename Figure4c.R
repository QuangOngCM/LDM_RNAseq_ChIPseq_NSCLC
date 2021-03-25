library("ChIPseeker")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("clusterProfiler")
library("org.Hs.eg.db")


PRC2.regions <- "data/Tumor-altered-overlapEZH2-SUZ12_TFBS.bed"

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
PRC2.regions.anno <- annotatePeak(PRC2.regions, tssRegion=c(-1000, 1000),TxDb=txdb, annoDb="EnsDb.Hsapiens.v75")
PRC2.regions.anno.df <- data.frame(PRC2.regions.anno)
PRC2.regions.anno.df <- PRC2.regions.anno.df[PRC2.regions.anno.df$annotation == "Promoter",]
g <- unique(PRC2.regions.anno.df$geneId)
ego <- enrichGO(gene          = g,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                #qvalueCutoff  = 0.05,
                readable      = TRUE)
#svg("Figure4c.svg", width = 8, height = 9)
barplot(ego,showCategory = 10) + ggtitle("Gene Ontology (BP) Enrichment Analysis")
#dev.off() 