#!/usr/bin/env Rscript
#
# =============================================================================
#
# Description:
#   This script performs differential gene expression (DGE) and subsequent
#   functional enrichment analysis for specific cell types (Fibroblasts,
#   Cardiomyocytes, and Endothelial cells). 
#
# Analysis Workflow:
#   For each cell type (Fibroblasts, Cardiomyocytes, and Endothelial cells):
#     1.  Perform DGE analysis for the following four comparisons:
#         - SQ vs. NT
#         - WT vs. R404Q
#         - SQ vs. R404Q
#         - R404Q vs. NT
#     2.  For each comparison, identify the sets of significantly upregulated
#         and downregulated genes.
#     3.  Perform Gene Ontology (GO) and KEGG pathway enrichment analyses
#         separately on the upregulated and downregulated gene sets.
#     4.  Visualize the results using volcano plots for DGE and dot plots
#         or bar plots for enrichment results.
#     5.  Generate heatmaps to visualize the expression patterns of the 
#         differentially expressed genes across the relevant samples.
#
# =============================================================================


# Load packages and colorSet ----
for (p in c("tidyverse", "Seurat", "stringr", "ggrepel", "EnhancedVolcano",
            "biomaRt", "clusterProfiler", "org.Mm.eg.db")) {
  if (!require(p, character.only = T)) {
    stop(paste0("Please install ", p))
  } else {
    suppressMessages(library(p, quietly = T, character.only = T))
  }
}
source("~/project/EpiAllele/script/snRNAseq/06_enrichment_funtions.R")

data.path <- "~/project/EpiAllele/result/snRNAseq/integration/"


# run DEG ----
result.path <- "~/project/EpiAllele/result/snRNAseq/DE/"
if (!dir.exists(result.path)) {
  dir.create(result.path)
}

data.harmony <- readRDS(paste0(data.path, "harmony.annotation.rds"))
extract_celltypes <- c("Fibroblasts", "Endothelial cells", "Cardiomyocytes")
Idents(data.harmony) <- "cellType"
#NT vs Sg； WT vs R404Q： R404Q vs sg； R404Q vs Nt
group_list <- list(c("SQ", "NT"),
                   c("WT", "R404Q"),
                   c("SQ", "R404Q"),
                   c("R404Q", "NT"))

for (celltype in extract_celltypes) {
  subset.data <- subset(data.harmony, cellType == celltype)
  Idents(subset.data) <- "dataset"
  DefaultAssay(subset.data) <- "RNA"
  DEGs_list <- list()

  for (group in group_list) {
    group1 <- group[1]
    group2 <- group[2]
    group_name <- paste0(group1, ".vs.", group2)
    #Perform DEG analysis
    DEGs <- FindMarkers(subset.data, ident.1 = group1, ident.2 = group2, 
                        logfc.threshold = 0, min.pct = 0.01)
    DEGs$gene <- rownames(DEGs)
    DEGs <- DEGs %>% 
      mutate(group = ifelse((p_val_adj < 0.05 & avg_log2FC > 0.5), "UP",
                            ifelse((p_val_adj < 0.05 & avg_log2FC < -0.5), "DOWN", "NOT-SIG"))) %>% 
      mutate(group = ifelse(gene == "Myh6", "Myh6", group)) %>% 
      arrange(group == "Myh6")
    top10_genes <- DEGs %>% 
      filter(group %in% c("DOWN", "UP")) %>% 
      arrange(p_val_adj, -abs(avg_log2FC)) %>% 
      do(head(., n = 10)) %>% 
      pull(gene)
    ### Create volcano plot
    keyvals <- ifelse(
      DEGs$group == "UP", "red3",
      ifelse(DEGs$group == "DOWN", "blue4", 
            ifelse(DEGs$group == "Myh6", "green4", "gray"))
    )
    names(keyvals) <- DEGs$group

    p <- EnhancedVolcano(DEGs, x = "avg_log2FC", y = "p_val_adj",
                        title = group_name,
                        lab = DEGs$gene, 
                        pCutoff = 0.05, 
                        FCcutoff = 0.5,
                        max.overlaps = 6, 
                        selectLab = top10_genes, 
                        colCustom = keyvals,
                        drawConnectors = TRUE,
                        widthConnectors = 0.5,
                        colConnectors = 'black')
    pdf(paste0(result.path, str_replace(celltype, " ", "."), ".", group_name, ".VolcanoPlot.pdf"), width = 6, height = 6)
    print(p)
    dev.off()
    DEGs$celltype <- celltype
    row.names(DEGs) <- NULL
    write.csv(DEGs, 
              paste0(result.path, str_replace(celltype, " ", "."), ".", group_name, ".DEresult.csv"))
    DEGs$comp <- group_name
    DEGs_list[[group_name]] <- DEGs
  }
  DEGs_allgroup <- do.call(rbind, DEGs_list)
  row.names(DEGs_allgroup) <- NULL
  write.csv(DEGs_allgroup, paste0(result.path, str_replace(celltype, " ", "."), ".DEG.csv"))
}



# enrichment ----
enrich.path <- "~/project/EpiAllele/result/snRNAseq/enrichment/"
if (!dir.exists(enrich.path)) {
  dir.create(enrich.path)
}
species <- "mouse"

for (celltype in extract_celltypes) {
  DEGs <- read.csv(paste0(result.path, str_replace(celltype, " ", "."), ".DEG.csv"))
  for (group in group_list) {
    group_name <- paste0(group[1], ".vs.", group[2])
    DEGs_group <- DEGs %>% filter(comp == group_name)
    up_genes <- DEGs_group %>% 
      filter(group == "UP") %>% 
      pull(gene)
    down_genes <- DEGs_group %>% 
      filter(group == "DOWN") %>% 
      pull(gene)
    
    if (length(up_genes) > 20) {
      GO_enrichment(gene=up_genes, name = paste0(celltype, ".", group_name, ".up"), 
                    title = paste0(group_name, " (up-regulated)"), 
                    results.path = enrich.path, species = species)
      KEGG_enrichment(up_genes, name = paste0(celltype, ".", group_name, ".up"), 
                      title = paste0(group_name, " (up-regulated)"), 
                      results.path = enrich.path, species = species)
    } else {
      print(paste0(group, " up-regulated genes are less than 20"))
    }
    
    if (length(down_genes) > 20) {
      GO_enrichment(down_genes, name = paste0(celltype, ".", group_name,  ".down"), 
                    title = paste0(group_name, " (down-regulated)"), 
                    results.path = enrich.path, species = species)
      KEGG_enrichment(down_genes, name = paste0(celltype, ".", group_name,  ".down"), 
                      title = paste0(group_name, " (down-regulated)"), 
                      results.path = enrich.path, species = species)
    } else {
      print(paste0(group, " down-regulated genes are less than 20"))
    }
  }
}


# expression heatmap ----
DE.path <- "~/project/EpiAllele/result/snRNAseq/DE/"
result.path <- "~/project/EpiAllele/result/snRNAseq/DE_heatmap/"
if (!dir.exists(result.path)) {
  dir.create(result.path)
}

group_list <- list(c("R404Q", "NT"),
                   c("SQ", "R404Q"), 
                   c("SQ", "NT"),
                   c("WT", "R404Q")
)

for (celltype in extract_celltypes) {
  subset.data <- subset(data.harmony, cellType == celltype)
  DEGs <- read.csv(paste0(DE.path, str_replace(celltype, " ", "."), ".DEG.csv"))
  geneset <- c()
  
  for (group in group_list) {
    group_name <- paste0(group[1], ".vs.", group[2])
    DEGs_group <- DEGs %>% filter(comp == group_name)
    up_gene <- DEGs_group %>% 
      filter(group == "UP") %>% 
      pull(gene)
    down_gene <- DEGs_group %>% 
      filter(group == "DOWN") %>% 
      pull(gene)
    gene <- c(up_gene, down_gene)
    subset.data.avg <- AverageExpression(subset.data, assays = "RNA", slot = "data",
                                         features = gene, group.by = "dataset")
    subset.data.avg <- as.matrix(subset.data.avg$RNA)
    subset.data.avg <- subset.data.avg[,c("NT", "R404Q", "SQ", "WT")]
    pdf(paste0(result.path, str_replace(celltype, " ", "."), ".", group_name, ".showDEG.heatmap.pdf"), 
        width=5, height=8)
    p1 <- pheatmap(subset.data.avg, 
                   scale = "row",
                   cluster_rows = T, cluster_cols = F,
                   show_rownames = F,
                   color = colorRampPalette(c("blue","white","red"))(100),
                   border_color = F,
                   heatmap_legend_param = list(
                     title = "Expression level"
                   ))
    print(p1)
    dev.off()
    
    geneset <- c(geneset, gene)
  }
  print(paste0(celltype, ": ", length(unique(geneset))))
  subset.data.avg <- AverageExpression(subset.data, assays = "RNA", slot = "data",
                                       features = geneset, group.by = "dataset")
  subset.data.avg <- as.matrix(subset.data.avg$RNA)
  
  pdf(paste0(result.path, str_replace(celltype, " ", "."), ".showDEG.heatmap.pdf"), 
      width=5, height=8)
  p1 <- pheatmap(subset.data.avg, 
                 scale = "row",
                 cluster_rows = T, cluster_cols = F,
                 show_rownames = F,
                 color = colorRampPalette(c("blue","white","red"))(100),
                 border_color = F,
                 heatmap_legend_param = list(
                   title = "Expression level"
                 ))
  print(p1)
  dev.off()
}






