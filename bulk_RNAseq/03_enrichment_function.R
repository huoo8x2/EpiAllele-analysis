library(tidyverse)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

GO_enrichment <- function(gene, name, title, results.path, species){
  if (species == "human") {
    orgDB <- org.Hs.eg.db
  } else if(species == "mouse"){
    orgDB <- org.Mm.eg.db
  } else {
    stop("species should be human or mouse")
  }
  
  
  switched_entrez <- bitr(unique(gene), "SYMBOL", "ENTREZID", orgDB, drop = TRUE)
  GO_results <- enrichGO(gene=switched_entrez$ENTREZID, 
                         OrgDb=orgDB, 
                         pvalueCutoff=1, 
                         qvalueCutoff=1, 
                         ont="all", readable=T)
  GO_results@result <- GO_results[order(GO_results@result$p.adjust)]
  GO_results_df <- as.data.frame(GO_results)
  write.csv(GO_results_df, paste0(results.path, name, ".GO.results.csv"))
  
  colorSel="qvalue"
  showNum <- ifelse(nrow(GO_results)<5, nrow(GO_results), 5)
  
  bar <- barplot(GO_results, drop=TRUE, showCategory=showNum, 
                 color=colorSel,label_format=60, 
                 title = paste0("GO enrichment of ", title), split = "ONTOLOGY") +
    facet_grid(ONTOLOGY~., scale = "free")
  
  bub <- dotplot(GO_results, showCategory=showNum, orderBy="GeneRatio", 
                 color=colorSel, label_format=60,
                 title = paste0("GO enrichment of ", title), split = "ONTOLOGY") +
    facet_grid(ONTOLOGY~., scale = "free")
  
  
  GO_results_df_top5 <- GO_results_df %>% 
    group_by(ONTOLOGY) %>% 
    do(head(., n = 5))
  GO_results_df_top5$Description <- as.factor(GO_results_df_top5$Description)
  GO_results_df_top5$Description <- fct_inorder(GO_results_df_top5$Description)
  p <- ggplot(GO_results_df_top5, aes(Description, Count,fill=ONTOLOGY))+
    geom_bar(stat = "identity")+
    geom_text(aes(label=Count, y=Count+0.2),size=3)+
    coord_flip()+
    labs(x='',y='Gene count', title = paste0("GO enrichment of ", title))+
    scale_fill_manual(values = c('#852f88',
                                 '#eb990c',
                                 '#0f8096'))+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 10))
  
  
  pdf(paste0(results.path, name, ".GO.result.pdf"), width = 10, height = 7)
  print(bar)
  print(bub)
  print(p)
  dev.off()
}

KEGG_enrichment <- function(gene, name, title, results.path, species){
  if (species == "human") {
    orgDB <- org.Hs.eg.db
    organism <- "hsa"
  } else if(species == "mouse"){
    orgDB <- org.Mm.eg.db
    organism <- "mmu"
  } else {
    stop("species should be human or mouse")
  }
  
  switched_entrez <- bitr(unique(gene), "SYMBOL", "ENTREZID", orgDB,drop = TRUE)
  KEGG_results <- enrichKEGG(gene=switched_entrez$ENTREZID, 
                             organism = organism,
                             pvalueCutoff=1, 
                             qvalueCutoff=1)
  KEGG_results@result <- KEGG_results[order(KEGG_results@result$p.adjust)]
  KEGG_results_df <- as.data.frame(KEGG_results)
  write.csv(KEGG_results_df, paste0(results.path, name, ".KEGG.results.csv"))
  
  colorSel="qvalue"
  showNum <- ifelse(nrow(KEGG_results)<10, nrow(KEGG_results), 10)
  bar <- barplot(KEGG_results, drop=TRUE, showCategory=showNum, 
                 color=colorSel,label_format=60, 
                 title = paste0("KEGG enrichment of ", title))
  
  bub <- dotplot(KEGG_results, showCategory=showNum, orderBy="GeneRatio", 
                 color=colorSel, label_format=60,
                 title = paste0("KEGG enrichment of ", title))
  
  pdf(paste0(results.path, name, ".KEGG.result.pdf"), width = 10, height = 6)
  print(bar)
  print(bub)
  dev.off()
}
