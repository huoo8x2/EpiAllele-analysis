#!/usr/bin/env Rscript
#
# =============================================================================
#
# Description:
#   This script performs detailed sub-clustering analysis on specific cell
#   lineages (Fibroblasts, Cardiomyocytes, and Endothelial cells). For each
#   cell type, it subsets the data, re-integrates across samples to remove
#   batch effects, performs subtype annotation, and analyzes the resulting
#   subtype compositions between samples.
#
# Analysis Workflow:
#   For each major cell type (Fibroblasts, Cardiomyocytes, Endothelial cells):
#     1. Subset the cells from the main integrated Seurat object.
#     2. Re-run integration (e.g., Harmony) on the subsetted data to correct
#        for batch effects specific to that lineage.
#     3. Perform fine-grained clustering and detailed subtype annotation using
#        known markers.
#     4. Calculate and visualize the proportions of each subtype across the
#        different original samples.
#     5. Analyze and visualize the similarity of subtype compositions between
#        samples.
#
# =============================================================================


# Load packages and colorSet ----
for (p in c("tidyverse", "Seurat", "ggpubr", "harmony", "ComplexHeatmap", "ggalluvial",
            "ggrepel", "clustree", "homologene", "GeneOverlap", "reshape2")) {
  if (!require(p, character.only = T)) {
    stop(paste0("Please install ", p))
  } else {
    suppressMessages(library(p, quietly = T, character.only = T))
  }
}
source("~/project/EpiAllele/script/snRNAseq/snRNA_colorSet.R")



data.path <- "~/project/EpiAllele/result/snRNAseq/integration/"
result.path <- "~/project/EpiAllele/result/snRNAseq/subcluster/"
if (!dir.exists(result.path)) {
  dir.create(result.path)
}

# Subclustering analysis ----
## for Fibroblasts, Cardiomyocytes and Endothelial cells
data.harmony <- readRDS(paste0(data.path, "harmony.annotation.rds"))
extract_celltypes <- c("Fibroblasts", "Endothelial cells", "Cardiomyocytes")
Idents(data.harmony) <- "cellType"
## runHarmony ----
for (celltype in extract_celltypes) {
  data <- subset(data.harmony, cellType == celltype)
  subset.data <- CreateSeuratObject(counts = data@assays$RNA@counts, 
                                    project = celltype, meta.data = data@meta.data)
  subset.data <- subset.data %>% 
    NormalizeData() %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(vars.to.regress = c("percent.mt", "S.Score", "G2M.Score")) %>% 
    RunPCA(npcs = 50, verbose = FALSE) %>% 
    RunHarmony(group.by.vars = "dataset", theta=2)
  
  ##Select optimal PCs for dimensional reduction
  pdf(paste0(result.path, "harmony.", celltype, ".ElbowPlot.pdf"),
      width = 7,height = 5)
  ElbowPlot(subset.data, ndims = 50)
  dev.off()
  
  subset.data <- subset.data %>% 
    RunUMAP(reduction = "pca", dims = 1:20) %>% 
    FindNeighbors(reduction = "pca", dims = 1:20)
  
  ##Select optimal resolution according to clustree
  # subset.res <- subset.data
  # subset.res <- FindClusters(
  #   object = subset.res,
  #   resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
  # )
  # pdf(paste0(result.path, "harmony.", celltype, ".clustree.pdf"),
  #     width = 12,height = 10)
  # p <- clustree(subset.res@meta.data, prefix = "RNA_snn_res.")
  # print(p)
  # dev.off()
  
  ##clustering
  rm(subset.res)
  subset.data <- FindClusters(subset.data, resolution = 0.4)
  Idents(subset.data) <- "seurat_clusters"
  
  pdf(paste0(result.path,"harmony.", celltype, ".umap.pdf"),width = 8,height = 6)
  p1 <- DimPlot(subset.data,reduction = "umap", group.by = "seurat_clusters", 
                label = T, label.size = 3)
  print(p1)
  p2 <- DimPlot(subset.data,reduction = "umap", group.by = "dataset")+
    scale_color_d3("category20")
  print(p2)
  if (celltype == "Fibroblasts") {
    p <- FeaturePlot(subset.data, features = "Postn")
    print(p)
    p <- FeaturePlot(subset.data, features = "Meox1")
    print(p)
  }
  if (celltype == "Cardiomyocytes") {
    p <- FeaturePlot(subset.data, features = "AAV-C")
    print(p)
    p <- FeaturePlot(subset.data, features = "AAV-N")
    print(p)
  }
  dev.off()
  
  pdf(paste0(result.path, "harmony.", celltype, ".bySample.umap.pdf"),width = 14,height = 10)
  p1 <- DimPlot(subset.data, reduction = "umap", group.by = "seurat_clusters",
                split.by = "dataset", ncol = 2)
  print(p1)
  if (celltype == "Fibroblasts") {
    p <- FeaturePlot(subset.data, ncol = 2, features = "Postn", split.by = "dataset")+
      patchwork::plot_layout(ncol = 2, nrow = 2)
    print(p)
    p <- FeaturePlot(subset.data, ncol = 2, features = "Meox1", split.by = "dataset")+
      patchwork::plot_layout(ncol = 2, nrow = 2)
    print(p)
  }
  if (celltype == "Cardiomyocytes") {
    p <- FeaturePlot(subset.data, ncol = 2, features = "AAV-C", split.by = "dataset")+
      patchwork::plot_layout(ncol = 2, nrow = 2)
    print(p)
    p <- FeaturePlot(subset.data, ncol = 2, features = "AAV-N", split.by = "dataset")+
      patchwork::plot_layout(ncol = 2, nrow = 2)
    print(p)
  }
  dev.off()
  
  saveRDS(subset.data, paste0(result.path, "harmony.", celltype, ".rds"))
  
}
## CCA integration ----
for (celltype in extract_celltypes) {
  data <- subset(data.harmony, cellType == celltype)
  subset.data <- CreateSeuratObject(counts = data@assays$RNA@counts, 
                                    project = celltype, meta.data = data@meta.data)
  data.list <- SplitObject(subset.data, split.by = "dataset")
  
  ### normalize and identify variable features for each dataset independently
  data.list.CCA <- lapply(X = data.list, FUN = function(x) {
    x <- x %>% 
      NormalizeData() %>% 
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
  })
  
  ### Select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = data.list.CCA)
  
  ### Perform integration
  anchors <- FindIntegrationAnchors(object.list = data.list.CCA, anchor.features = features)
  data.integrated <- IntegrateData(anchorset = anchors)
  DefaultAssay(data.integrated) <- "integrated"
  
  ### Run the standard workflow for visualization and clustering
  data.integrated <- data.integrated %>% 
    ScaleData(vars.to.regress = c("percent.mt", "S.Score", "G2M.Score")) %>% 
    RunPCA(npcs = 50, verbose = FALSE)
  
  ### Select optimal PCs for dimensional reduction
  pdf(paste0(result.path,"CCA.integrated.", celltype, ".ElbowPlot.pdf"),
      width = 7,height = 5)
  ElbowPlot(data.integrated, ndims = 50)
  dev.off()
  
  # ##Select optimal resolution according to clustree
  # data.integrated <- data.integrated %>% 
  #   RunUMAP(reduction = "pca", dims = 1:20) %>% 
  #   FindNeighbors(reduction = "pca", dims = 1:20)
  # 
  # subset.res <- data.integrated
  # subset.res <- FindClusters(
  #   object = subset.res,
  #   resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8)
  # )
  # pdf(paste0(result.path, "CCA.integrated.", celltype, ".clustree.pdf"),
  #     width = 12,height = 10)
  # p <- clustree(subset.res@meta.data, prefix = "integrated_snn_res.")
  # print(p)
  # dev.off()
  
  ### Clustering
  data.integrated <- data.integrated %>%
    RunUMAP(reduction = "pca", dims = 1:20) %>%
    FindNeighbors(reduction = "pca", dims = 1:20) %>%
    FindClusters(resolution = 0.3)

  ### UMAP for visualization
  DefaultAssay(data.integrated) <- "RNA"

  pdf(paste0(result.path,"CCA.integrated.", celltype, ".umap.pdf"),width = 8,height = 6)
  p1 <- DimPlot(data.integrated,reduction = "umap", group.by = "seurat_clusters",
                label = T, label.size = 3)
  print(p1)
  p2 <- DimPlot(data.integrated,reduction = "umap", group.by = "dataset")+
    scale_color_d3("category20")
  print(p2)
  if (celltype == "Fibroblasts") {
    p <- FeaturePlot(data.integrated, reduction = "umap", features = "Postn")
    print(p)
    p <- FeaturePlot(data.integrated, reduction = "umap", features = "Meox1")
    print(p)
  }
  if (celltype == "Cardiomyocytes") {
    p <- FeaturePlot(data.integrated, reduction = "umap", features = "AAV-C")
    print(p)
    p <- FeaturePlot(data.integrated, reduction = "umap", features = "AAV-N")
    print(p)
  }
  dev.off()

  pdf(paste0(result.path,"CCA.integrated.", celltype, ".bySample.umap.pdf"),width = 14,height = 10)
  p1 <- DimPlot(data.integrated, reduction = "umap", group.by = "seurat_clusters",
                split.by = "dataset", ncol = 2)
  print(p1)
  if (celltype == "Fibroblasts") {
    p <- FeaturePlot(data.integrated, ncol = 2, features = "Postn", split.by = "dataset")+
      patchwork::plot_layout(ncol = 2, nrow = 2)
    print(p)
    p <- FeaturePlot(data.integrated, ncol = 2, features = "Meox1", split.by = "dataset")+
      patchwork::plot_layout(ncol = 2, nrow = 2)
    print(p)
  }
  if (celltype == "Cardiomyocytes") {
    p <- FeaturePlot(data.integrated, ncol = 2, features = "AAV-C", split.by = "dataset")+
      patchwork::plot_layout(ncol = 2, nrow = 2)
    print(p)
    p <- FeaturePlot(data.integrated, ncol = 2, features = "AAV-N", split.by = "dataset")+
      patchwork::plot_layout(ncol = 2, nrow = 2)
    print(p)
  }
  dev.off()

  saveRDS(data.integrated, paste0(result.path, "CCA.integrated.", celltype, ".rds"))
}


# FindAllMarkers ----
for (celltype in extract_celltypes) {
  data.integrated <- readRDS(paste0(result.path, "CCA.integrated.", celltype, ".rds"))
  
  DefaultAssay(data.integrated) <- "RNA"
  Idents(data.integrated) <- "seurat_clusters"
  markers <- FindAllMarkers(data.integrated,logfc.threshold = 0)
  markers$expression <- ifelse(markers$avg_log2FC>0, "up-regulated", "down-regulated")
  rownames(markers) <- NULL
  write.csv(markers, file = paste0(result.path, "CCA.", celltype, ".subcluster.markers.csv"))
}

# Subtype annotation ----
# ## Calculate the Jaccard similarity of top100 marker genes across clusters ----
# jaccard_matrix_list <- list()
# for (celltype in extract_celltypes) {
#   markers <- read.csv(paste0(result.path, "CCA.", celltype, ".subcluster.markers.csv"),
#                       header=T, row.names = 1)
#   top.markers <- markers %>% 
#     filter(expression=="up-regulated") %>% 
#     group_by(cluster) %>% 
#     do(head(., n = 100))
#   top.markers.list <- split(top.markers$gene, top.markers$cluster)
#   gom.obj <- newGOM(top.markers.list,top.markers.list)
#   jaccard_matrix <- getMatrix(gom.obj, name = "Jaccard")
#   jaccard_matrix_list[[celltype]] <- jaccard_matrix
#   pdf(paste0(result.path, celltype, "_markergene_jaccrad_similarity.pdf"),
#       width = 8, height = 6)
#   p <- pheatmap(jaccard_matrix,
#                 heatmap_legend_param = list(
#                   title = "Jaccard Similarity"))
#   draw(p)
#   dev.off()
# }

## Fibroblasts ----
celltype <- "Fibroblasts"
data.integrated <- readRDS(paste0(result.path, "CCA.integrated.", celltype, ".rds"))
data.integrated$subtype <- paste0("Fibroblasts", "_", as.numeric(data.integrated$seurat_clusters))
saveRDS(data.integrated, paste0(result.path, "CCA.integrated.", celltype, ".rds"))

## Endothelial cells ----
celltype <- "Endothelial cells"
data.integrated <- readRDS(paste0(result.path, "CCA.integrated.", celltype, ".rds"))
data.integrated$subtype <- paste0("Endothelial_cells", "_", as.numeric(data.integrated$seurat_clusters))
saveRDS(data.integrated, paste0(result.path, "CCA.integrated.", celltype, ".rds"))

## Cardiomyocytes ----
celltype <- "Cardiomyocytes"
data.integrated <- readRDS(paste0(result.path, "CCA.integrated.", celltype, ".rds"))
subtype <- c(1, 2, 3, 4, 5, 2, 6, 7, 8, 9, 10)
Idents(data.integrated) <- "seurat_clusters"
names(subtype) <- levels(data.integrated)
data.integrated <- RenameIdents(data.integrated, subtype)
data.integrated$subtype <- data.integrated@active.ident
data.integrated$subtype <- paste0("Cardiomyocytes", "_", as.numeric(data.integrated$subtype))
data.integrated$subtype <- factor(data.integrated$subtype, levels=paste0("Cardiomyocytes_", 1:10))
saveRDS(data.integrated, paste0(result.path, "CCA.integrated.", celltype, ".rds"))

## FindAllMarkers ----
for (celltype in extract_celltypes){
  data.integrated <- readRDS(paste0(result.path, "CCA.integrated.", celltype, ".rds"))
  DefaultAssay(data.integrated) <- "RNA"
  Idents(data.integrated) <- "subtype"
  markers <- FindAllMarkers(data.integrated,logfc.threshold = 0)
  markers$expression <- ifelse(markers$avg_log2FC>0, "up-regulated", "down-regulated")
  rownames(markers) <- NULL
  write.csv(markers, file = paste0(result.path, "CCA.", celltype, ".subtype.markers.csv"))
}

## UMAP for visualization ----
for (celltype in extract_celltypes) {
  data.integrated <- readRDS(paste0(result.path, "CCA.integrated.", celltype, ".rds"))
  
  DefaultAssay(data.integrated) <- "RNA"
  pdf(paste0(result.path,"CCA.integrated.", celltype, ".umap.pdf"),width = 8,height = 6)
  p1 <- DimPlot(data.integrated,reduction = "umap", group.by = "seurat_clusters",
                label = T, label.size = 3)
  print(p1)
  p2 <- DimPlot(data.integrated,reduction = "umap", group.by = "subtype")+
    scale_color_manual(values = colorSet[[celltype]])
  print(p2)
  p3 <- DimPlot(data.integrated,reduction = "umap", group.by = "dataset")+
    scale_color_manual(values = colorSet$dataset)
  print(p3)
  if (celltype == "Fibroblasts") {
    p <- FeaturePlot(data.integrated, reduction = "umap", features = "Postn")
    print(p)
    p <- FeaturePlot(data.integrated, reduction = "umap", features = "Meox1")
    print(p)
  }
  if (celltype == "Cardiomyocytes") {
    p <- FeaturePlot(data.integrated, reduction = "umap", features = "AAV-C")
    print(p)
    p <- FeaturePlot(data.integrated, reduction = "umap", features = "AAV-N")
    print(p)
  }
  dev.off()
  
  pdf(paste0(result.path,"CCA.integrated.", celltype, ".bySample.umap.pdf"),width = 14,height = 10)
  p1 <- DimPlot(data.integrated, reduction = "umap", group.by = "seurat_clusters",
                split.by = "dataset", ncol = 2, pt.size = 1.0)
  print(p1)
  p2 <- DimPlot(data.integrated, reduction = "umap", group.by = "subtype",
                split.by = "dataset", ncol = 2, pt.size = 1.0)+
    scale_color_manual(values = colorSet[[celltype]])
  print(p2)
  if (celltype == "Fibroblasts") {
    p <- FeaturePlot(data.integrated, ncol = 2, features = "Postn", split.by = "dataset")+
      patchwork::plot_layout(ncol = 2, nrow = 2)
    print(p)
    p <- FeaturePlot(data.integrated, ncol = 2, features = "Meox1", split.by = "dataset")+
      patchwork::plot_layout(ncol = 2, nrow = 2)
    print(p)
  }
  if (celltype == "Cardiomyocytes") {
    p <- FeaturePlot(data.integrated, ncol = 2, features = "AAV-C", split.by = "dataset")+
      patchwork::plot_layout(ncol = 2, nrow = 2)
    print(p)
    p <- FeaturePlot(data.integrated, ncol = 2, features = "AAV-N", split.by = "dataset")+
      patchwork::plot_layout(ncol = 2, nrow = 2)
    print(p)
  }
  dev.off()
}



# Dotplot for top5 markers expression for each subtypes in interested celltypes ----
extract_celltypes <- c("Fibroblasts", "Endothelial cells", "Cardiomyocytes")
for (celltype in extract_celltypes) {
  data.integrated <- readRDS(paste0(result.path, "CCA.integrated.", celltype, ".rds"))
  markers <- read.csv(paste0(result.path, "CCA.", celltype, ".subtype.markers.csv"),
                      header=T, row.names = 1)
  markers$cluster <- factor(markers$cluster, levels = levels(data.integrated$subtype))
  top.markers <- markers %>% 
    filter(expression=="up-regulated") %>% 
    group_by(cluster) %>% 
    do(head(., n = 5))
  data.integrated$subtype <- factor(data.integrated$subtype)
  Idents(data.integrated) <- "subtype"
  pdf(paste0(result.path, celltype, ".top5_markers_dotplot.pdf"),width = 20,height = 10)
  p <- DotPlot(data.integrated,assay = "RNA",features = unique(top.markers$gene),
               cols = c(low="blue",high = "red"))+
    theme(axis.text.x = element_text(angle = 45,size = 18,vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 24),
          axis.title = element_text(face = "bold",size=30),
          legend.text = element_text(size = 18),
          legend.title=element_text(face = "bold",size = 20))
  print(p)
  dev.off()
}

# Subtype proportion statistics ----
prop_path <- "~/project/EpiAllele/result/snRNAseq/proportion_stat/"
if (!dir.exists(prop_path)) {
  dir.create(prop_path)
}


for (celltype in extract_celltypes) {
  data.integrated <- readRDS(paste0(result.path, "CCA.integrated.", celltype, ".rds"))
  meta <- data.integrated@meta.data
  celltype_stat <- table(meta$dataset, meta$subtype) %>% melt()
  colnames(celltype_stat) <- c("Dataset","CellType","Number")
  
  ## Calculate proportion ----
  celltype_stat <- celltype_stat %>%
    group_by(Dataset) %>%
    mutate(Proportion = Number / sum(Number))
  celltype_stat$Dataset <- factor(celltype_stat$Dataset, levels = c("NT", "R404Q", "SQ", "WT"))
  
  ## Generate bar plot and alluvial plot ----
  pdf(paste0(prop_path, celltype, "_subtype_stack_barplot.pdf"),width = 9,height = 6)
  ### Bar plot ----
  p <- ggplot(celltype_stat, aes(x=Dataset,y=Proportion,fill=CellType))+
    geom_bar(stat = "identity", width=0.8,position="fill")+
    theme_classic()+
    labs(x="Dataset",y="Ratio",fill="Cell type",
         title = "Celltype Proportion Bar Plot")+
    theme(axis.title = element_text(size=16,face="bold"),
          legend.title  = element_text(size=16,face="bold"),
          legend.text = element_text(size=14),
          axis.text = element_text(size=15),
          plot.title = element_text(size=18,face="bold"))+
    scale_fill_manual(values=colorSet[[celltype]])
  print(p)
  
  ### Alluvial plot ----
  p <- ggplot(celltype_stat,aes(x=Dataset,y=Proportion,fill=CellType,
                                stratum = CellType, alluvium = CellType))+
    geom_bar(stat = "identity", width=0.35,position="fill")+
    geom_alluvium(alpha = 0.3, knot.pos = 0)+
    theme_classic()+
    labs(x="Dataset",y="Ratio",fill="Cell type",
         title = "Celltype Proportion Alluvial Plot")+
    theme(axis.title = element_text(size=16,face="bold"),
          legend.title  = element_text(size=16,face="bold"),
          legend.text = element_text(size=14),
          axis.text = element_text(size=15),
          plot.title = element_text(size=18,face="bold"))+
    scale_fill_manual(values=colorSet[[celltype]])
  print(p)
  dev.off()
  
  ## Correlation of celltype proportion between groups ----
  prop_matrix <- celltype_stat %>% 
    select(Dataset, CellType, Proportion) %>% 
    pivot_wider(names_from = Dataset, values_from = Proportion) %>%
    column_to_rownames("CellType")
  
  cor_matrix <- cor(prop_matrix, method = "spearman")
  
  pdf(paste0(prop_path, celltype, ".subtype.proportion.cor.heatmap.pdf"), width = 4, height = 3)
  p1 <- pheatmap(cor_matrix, show_rownames = T,
                 cluster_rows = F, cluster_cols = F,
                 color = colorRampPalette(c("white", "#D5676D"))(100),
                 breaks = seq(0.5, 1, length.out = 100),
                 display_numbers = T, number_color = "black",
                 heatmap_legend_param = list(
                   title = "correlation"
                 ))
  draw(p1)
  dev.off()
  
}








