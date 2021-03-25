library("ChIPseeker")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("EnsDb.Hsapiens.v75")
library("clusterProfiler")
library("org.Hs.eg.db")

source("RNA-seq_DEGsAnalysis.R")
source("Figure2b-c.R")

geneList <- res.rna.all$log2FoldChange
names(geneList) <- res.rna.all$entrezid
geneList <- sort(geneList, decreasing = TRUE)

altered.noTSS.geneList <- read.table(file = "data/GeneHancer-TargetGenes.txt")
altered.noTSS.gene <- unique(as.character(na.omit(mapIds(EnsDb.Hsapiens.v75,keys = as.character(altered.noTSS.geneList$V1),column = "ENTREZID", 
                                                         keytype = "SYMBOL", multiVals = "first"))))

LUAD.rna.DEGs <- na.omit(SigDEG.LUAD$entrezid)
LUSC.rna.DEGs <- na.omit(SigDEG.LUSC$entrezid)

a1 <- intersect(altered.noTSS.gene,LUAD.rna.DEGs)
a2 <- intersect(altered.noTSS.gene,LUSC.rna.DEGs)
a3 <- unique(union(a1,a2)) #126 genes were differentially expressed in LUAD and LUSC

ego <- enrichGO(gene          = a3,
                OrgDb         = org.Hs.eg.db,
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                readable      = TRUE)
#svg("SupplementaryFigureS7.svg", width = 8, height = 8)
cnetplot(ego,  showCategory = 5, categorySize="pvalue", foldChange=geneList)
#dev.off()
