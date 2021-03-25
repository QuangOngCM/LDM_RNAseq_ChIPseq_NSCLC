library("gghighlight")
library("ggplot2")
source("RNA-seq_DEGsAnalysis.R")

SigDEG.all.genes <- as.character(na.omit(SigDEG.all$entrezid))
geneList <- SigDEG.all$log2FoldChange
names(geneList) <- SigDEG.all$entrezid

##Generate random 1,000 gene list from 4,118 DEGs. Each gene list contains 127 genes
random.list <- list()
for(i in 1:1000) {
  set.seed(i)
  random.list[[paste0("genelist", i)]] <- sample(SigDEG.all.genes,127)
}
#Run KEGG pathway analysis for each gene list
ego.list <- list()
for(i in 1:1000) {
  ego.list[[paste0("ego", i)]] <- as.data.frame(enrichGO(gene          = random.list[[i]],
                                                OrgDb         = org.Hs.eg.db,
                                                ont           = "MF",
                                                pAdjustMethod = "BH",
                                                pvalueCutoff  = 0.05,
                                                #qvalueCutoff  = 0.05,
                                                readable      = TRUE))$Description
}
#Making Figure S9a
lengths(ego.list)
m <-data.frame(table(lengths(ego.list)))
#svg("SupplementaryFigure9a.svg", width = 5, height = 4)
ggplot(data=m,  aes(x=Var1, y=Freq, group=1, label = Freq)) +
  geom_text(aes(label=Freq),hjust=0, vjust=-1)+
  geom_line(linetype = "dashed")+
  geom_point() +
  ylim(0,1000)+
  labs(y = "Number of gene list", x = "Number of enriched pathway(s)") +
  ggtitle("1,000 Random Gene Lists with n = 126")
#dev.off() 

#Making Figure S9b
table.kegg <- data.frame(table(as.character(unlist(ego.list))))
table.kegg$Freq1000 <- table.kegg$Freq/10
my_data <- as_tibble(table.kegg)
table.kegg.order <- my_data[,c(1,3)]
K0 <- data.frame("No enrichment", 81.5)
colnames(K0) <- c("Var1","Freq1000")
table.kegg.order <-rbind(table.kegg.order,K0)
table.kegg.order <- table.kegg.order %>% arrange(desc(Freq1000)) #my_data %>% arrange(desc(xxx))
table.test1 <- table.kegg.order

#svg("SupplementaryFigure9b1.svg", width = 11, height = 9)
table.test1 %>%
  ggplot(aes(x=reorder(Var1, Freq1000), y=Freq1000)) +
  geom_col() + 
  gghighlight(Freq1000 < 5,
              unhighlighted_params = list(size = 1, fill = "red")) +
  theme_minimal()+
  labs(y = "Percentage of occurrence (%)", x = "GO-MF pathways (Total = 221)") +
  ggtitle("1,000 Random Gene Lists with n = 310")+
  geom_hline(yintercept=5, linetype="dashed", color = "red") + 
  coord_flip()
#dev.off() 
#Making Figure S9b2 top10
table.test.fig9 <- table.test1[-1,]
table.test.fig9.top10 <- table.test.fig9[c(1:10),]
#svg("SupplementaryFigure9b2.svg", width = 9, height = 7)
table.test.fig9.top10 %>%
  ggplot(aes(x=reorder(Var1, Freq1000), y=Freq1000)) +
  geom_col() + 
  gghighlight(Freq1000 < 5,
              unhighlighted_params = list(size = 1, fill = "red")) +
  theme_minimal()+
  labs(y = "Percentage of occurrence (%)", x = "GO-MF pathways") +
  ggtitle("1,000 Random Gene Lists with n = 310")+
  theme(plot.title = element_blank(),
        axis.title.y = element_blank()) +
  geom_hline(yintercept=5, linetype="dashed", color = "red") + 
  coord_flip()
#dev.off() 