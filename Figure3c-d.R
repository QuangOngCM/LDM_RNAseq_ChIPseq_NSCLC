###Running RNA-seq analysis###
library("DESeq2")
library("ggplot2")
library("EnsDb.Hsapiens.v75")
library("AnnotationDbi")
library("ComplexHeatmap")
library("circlize")
library("ChIPseeker")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("enrichplot")
library("clusterProfiler")
library("org.Hs.eg.db")

source("RNA-seq_DEGsAnalysis.R")
source("Figure2b-c.R")


###Getting ChIP-seq analysis results###
altered_TSS <- "data/Tumor-altered-overlapTSS.bed"
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

peakAnno <- annotatePeak(altered_TSS, tssRegion=c(-1000, 1000),TxDb=txdb, annoDb="EnsDb.Hsapiens.v75")
altered.TSS.genes <- na.exclude(as.data.frame(peakAnno)$geneId)
dup.genes <- altered.TSS.genes[duplicated(altered.TSS.genes) == TRUE]
altered_TSS.df <- as.data.frame(peakAnno)
altered.TSS.genes <- unique(na.exclude(as.data.frame(peakAnno)$geneId))

###Making heatmap###
Sig.LUAD.genes <- na.omit(SigDEG.LUAD$entrezid)
Sig.LUSC.genes <- na.omit(SigDEG.LUSC$entrezid)

altered.TSS.genes.DEGs <- union(intersect(Sig.LUAD.genes,altered.TSS.genes), intersect(Sig.LUSC.genes,altered.TSS.genes))
altered.TSS.genes.DEGs.data <- subset(res.rna.all, res.rna.all$entrezid %in% altered.TSS.genes.DEGs)
altered.TSS.genes.DEGs.count <- subset(norm.count.RNA, rownames(norm.count.RNA) %in% rownames(altered.TSS.genes.DEGs.data))
altered.TSS.genes.DEGs.count <- altered.TSS.genes.DEGs.count[,c(11:30)]
colnames(altered.TSS.genes.DEGs.count) <- colnames(norm.count.RNA[,c(11:30)])
altered.TSS.genes.DEGs.ezid <- as.data.frame(na.omit(mapIds(EnsDb.Hsapiens.v75,keys = altered.TSS.genes.DEGs,column = "GENEID", keytype = "ENTREZID", multiVals = "first")))
altered.TSS.df <- as.data.frame(peakAnno)
altered.TSS.df$status <- ""
altered.TSS.df$status[1:188] <- "gained.TSS" 
altered.TSS.df$status[189:444] <- "loss.TSS" 

#Checking gene names
a1 <- unique(subset(altered.TSS.df, altered.TSS.df$geneId %in% altered.TSS.genes.DEGs))[,c(15,20)]
a1$entrez <- mapIds(EnsDb.Hsapiens.v75,keys = a1$geneId,column = "GENEID", keytype = "ENTREZID", multiVals = "first")
a2 <- unique(subset(a1, a1$entrez %in% rownames(altered.TSS.genes.DEGs.count))) 
subset(altered.TSS.genes.DEGs.count, !(rownames(altered.TSS.genes.DEGs.count) %in% a2$entrez)) 
b <- data.frame("729238","loss.TSS","ENSG00000185303")
colnames(b) <- colnames(a2)
a3 <- rbind(a2,b)
a3 <- a3[,c(3,2)]

#There are 311 (NSCLC subtype specific) DEGs with their TSSs overlaped with tumor-altered H3K4me3 regions as shown in 'altered.TSS.genes.DEGs.data'

altered.TSS.genes.DEGs.type <- data.frame(type = a3[,2])
rownames(altered.TSS.genes.DEGs.type) <- a3$entrez
altered.TSS.genes.DEGs.type <- altered.TSS.genes.DEGs.type[match(rownames(altered.TSS.genes.DEGs.count), rownames(altered.TSS.genes.DEGs.type)), ]

#Making Figure 3c
peak_type_col <- c("#91cf60","#f7fcb9")
names(peak_type_col) <- c("gained.TSS","loss.TSS")

typeRNA = c(rep("Normal", 10), rep("Tumor",10))
CelltypeRNA = c(rep(c(rep("squamous", 4),rep("adeno", 6)),2))

haRNA = HeatmapAnnotation(df = data.frame(Sample = typeRNA, Case = CelltypeRNA), 
                          col = list(Sample = c("Normal" =  "#3288bd", "Tumor" = "#d53e4f"),
                                     Case = c("squamous" = "#fc8d59", "adeno" = "#b35806")),
                          annotation_name_side = "left",
                          annotation_legend_param = list(
                            Sample = list(nrow = 2),
                            Case = list(nrow = 2)))

scale.chip <- as.data.frame(t(apply(as.matrix(altered.TSS.genes.DEGs.count),1, scale)),stringsAsFactors = FALSE)
colnames(scale.chip) <- colnames(altered.TSS.genes.DEGs.count)
hm1 <- Heatmap(scale.chip, name = "z-score",col = colorRamp2(c(-0.7,1.2,2.7,4.2,5.7), c("black", "red", "orange","yellow","green")),
               top_annotation = haRNA, cluster_columns = T, column_names_gp = gpar(fontsize = 8), clustering_distance_rows = "euclidean",
               show_row_names = FALSE, show_column_names = TRUE, row_km = 2, heatmap_legend_param = list(direction = "horizontal"))
hm2 <- Heatmap(altered.TSS.genes.DEGs.type, name = "Type", show_row_names = FALSE,
               show_column_names = FALSE, col = peak_type_col,heatmap_legend_param = list(nrow = 2),
               width = unit(5, "mm"))
#svg("Figure3c.svg", width = 7, height = 7.5)
draw(hm1 + hm2, merge_legend = T, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
#dev.off() 
#Making Figure 3d
geneList <- altered.TSS.genes.DEGs.data$log2FoldChange
names(geneList) <- altered.TSS.genes.DEGs.data$entrezid
kk <- enrichKEGG(gene         = altered.TSS.genes.DEGs,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
edox <- setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
#svg("Figure3d.svg", width = 9, height = 7.5)
cnetplot(edox, foldChange=geneList)
#dev.off() 