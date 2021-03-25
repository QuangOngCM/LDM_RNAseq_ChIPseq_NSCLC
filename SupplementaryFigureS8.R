library("gghighlight")

source("RNA-seq_DEGsAnalysis.R")

SigDEG.all.genes <- as.character(na.omit(SigDEG.all$entrezid))

##Generate random 1,000 gene list from 4,118 DEGs. Each gene list contains 310 genes
random.list <- list()
for(i in 1:1000) {
  set.seed(i)
  random.list[[paste0("genelist", i)]] <- sample(SigDEG.all.genes,310)
}
#Run KEGG pathway analysis for each gene list
kk.list <- list()
for(i in 1:1000) {
  kk.list[[paste0("kk", i)]] <- as.data.frame(enrichKEGG(gene = random.list[[i]],
                                                         organism     = 'hsa',
                                                         pvalueCutoff = 0.05))$Description
}
#Making Figure S8a
lengths(kk.list)
m <-data.frame(table(lengths(kk.list)))
#svg("SupplementaryFigureS8a.svg", width = 5, height = 4)
ggplot(data=m,  aes(x=Var1, y=Freq, group=1, label = Freq)) +
  geom_text(aes(label=Freq),hjust=0, vjust=-1)+
  geom_line(linetype = "dashed")+
  geom_point() +
  ylim(0,1000)+
  labs(y = "Number of gene list", x = "Number of enriched pathway(s)") +
  ggtitle("1,000 Random Gene Lists with n = 310")
#dev.off() 
##Making Figure S8b1
table.kegg <- data.frame(table(as.character(unlist(kk.list))))
table.kegg$Freq1000 <- table.kegg$Freq/10
my_data <- as_tibble(table.kegg)
table.kegg.order <- my_data[,c(1,3)]
K0 <- data.frame("No enrichment", 80)
colnames(K0) <- c("Var1","Freq1000")
table.kegg.order <-rbind(table.kegg.order,K0)
table.kegg.order <- table.kegg.order %>% arrange(desc(Freq1000)) #my_data %>% arrange(desc(xxx))
table.test1 <- table.kegg.order
#Checking KEGG pathways related to Pathways in cancer (hsa05200)
kegg.cancer <- read.csv("KEGG-pw-cancer.txt", header = F)
table.test.fig6c <- table.test1
table.test.fig6c$Kegg.cancer <- "no"
table.test.fig6c[table.test.fig6c$Var1 %in% kegg.cancer$V1,]$Kegg.cancer <- "yes"
#
#svg("SupplementaryFigureS8b1.svg", width = 7, height = 9)
table.test.fig6c %>%
  ggplot(aes(x=reorder(Var1, Freq1000), y=Freq1000, fill = Kegg.cancer)) +
  scale_fill_manual( values = c( "yes"="blue", "no"="gray" ), guide = FALSE )+
  geom_col() + 
  gghighlight(Freq1000 < 5,
              unhighlighted_params = list(size = 1, fill = "red")) +
  theme_minimal()+
  labs(y = "Percentage of occurrence (%)", x = "KEGG pathways (Total = 99)") +
  ggtitle("1,000 Random Gene Lists with n = 310")+
  geom_hline(yintercept=5, linetype="dashed", color = "red") + 
  coord_flip()
#dev.off() 
#Making Figure S8b2 top10
table.test.fig6c2 <- table.test.fig6c[-1,]
table.test.fig6c2.top10 <- table.test.fig6c2[c(1:10),]
#svg("SupplementaryFigureS8b2.svg", width = 9, height = 7)
table.test.fig6c2.top10 %>%
  ggplot(aes(x=reorder(Var1, Freq1000), y=Freq1000, fill = Kegg.cancer)) +
  scale_fill_manual( values = c( "yes"="blue", "no"="gray" ), guide = FALSE )+
  geom_col() + 
  gghighlight(Freq1000 < 5,
              unhighlighted_params = list(size = 1, fill = "red")) +
  theme_minimal()+
  labs(y = "Percentage of occurrence (%)", x = "KEGG pathways (Top 10)") +
  ggtitle("1,000 Random Gene Lists with n = 310")+
  theme(plot.title = element_blank(),
        axis.title.y = element_blank()) +
  geom_hline(yintercept=5, linetype="dashed", color = "red") + 
  coord_flip()
#dev.off() 