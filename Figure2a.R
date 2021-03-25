library("DESeq2")
library("ggplot2")
library("EnsDb.Hsapiens.v75")
library("AnnotationDbi")
library("ComplexHeatmap")
library("circlize")

source("RNA-seq_DEGsAnalysis.R")

#Making heatmap
typeRNA = c(rep("Stroma", 10),rep("Normal", 10), rep("Tumor",10))
case = c(rep(c(rep("squamous", 4),rep("adeno", 6)),3))

haRNA = HeatmapAnnotation(df = data.frame(Sample = typeRNA, Case = case), 
                          col = list(Sample = c("Stroma" = "#1a9850", "Normal" =  "#3288bd", "Tumor" = "#d53e4f"),
                                     Case = c("squamous" = "#fc8d59", "adeno" = "#b35806")),
                          annotation_name_side = "left",
                          annotation_legend_param = list(
                            Sample = list(nrow = 1),
                            Case = list(nrow = 1)))

scale.chip <- as.data.frame(t(apply(as.matrix(sig.nor.count.RNA),1, scale)),stringsAsFactors = FALSE)
colnames(scale.chip) <- colnames(sig.nor.count.RNA)

hm1 <- Heatmap(scale.chip, name = "z-Score", col = colorRamp2(c(-0.7,1.2,2.7,4.2,5.7), c("black", "red", "orange","yellow","green")),
               top_annotation = haRNA, cluster_columns = T, column_names_gp = gpar(fontsize = 8),
               clustering_distance_columns = "euclidean",clustering_distance_rows = "pearson",
               show_row_names = FALSE, show_column_names = TRUE, row_km = 1, heatmap_legend_param = list(direction = "horizontal"))

#svg(filename="Figure2a.svg", 
#    width=6, 
#    height=9.5, 
#    pointsize=12)
draw(hm1, merge_legend = T, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
#dev.off()