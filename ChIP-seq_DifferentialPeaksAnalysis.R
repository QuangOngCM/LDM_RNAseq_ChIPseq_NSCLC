library("EDASeq")
library("DESeq")
ChIP.RawCount <- readRDS("data/ChIP-seq_RawCounts.Rds")

peak.length <- read.table("data/Total-unified-H3K4me3_ChIP-seq_length_in_Kb.txt")

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

Chip.tpms <- apply(ChIP.RawCount, 2, function(x) tpm(x, peak.length$V2))
Chip.dataSet <-newSeqExpressionSet(counts = as.matrix(round(Chip.tpms)),
                                   phenoData = data.frame(conditions=c(rep("Normal",6),rep("Tumor",6),rep("Encode",6)),
                                                          row.names = colnames(Chip.tpms)))
typeChip = factor(c(rep("squamous",3),rep("adeno", 3),rep("squamous",3),rep("adeno", 3),rep("lung",6)))

design = data.frame(row.names = colnames (Chip.tpms), 
                    condition = c(rep("Normal",6),rep("Tumor",6),rep("Encode",6)),
                    cellType = typeChip)

cds = newCountDataSet(counts(Chip.dataSet),design)
sizeFactors(cds) <- rep(1,18)

Chip.Nor.count <- estimateDispersions(cds)

res.Chip <- nbinomTest(Chip.Nor.count, "Encode", "Tumor")

SigResChip <- subset(res.Chip, abs(res.Chip$log2FoldChange) >= 1 & res.Chip$log2FoldChange != "Inf" & res.Chip$log2FoldChange != "-Inf" 
                     & res.Chip$padj <= 0.05)
#We identified 645 promoters exhibiting differential H3K4me3 modifications between NSCLCs and ENCODE normal lung tissues (tumor-associated promoters). 