#!/usr/bin/env Rscript
#
# =============================================================================
#
# Description:
#   This script integrates single-cell RNA-seq datasets from multiple samples.
#   It performs and compares two batch correction methods (CCA and Harmony),
#   proceeds with the Harmony-corrected data to perform cell type annotation, 
#   and analyzes the resulting cell type compositions across samples.
#   
# Analysis Workflow:
#   1.  Merge individual, pre-processed Seurat objects from multiple samples.
#   2.  Perform batch correction using both Seurat's CCA workflow and Harmony.
#   3.  Visualize and compare the integration results to select the optimal
#       method (Harmony).
#   4.  Perform clustering and broad cell type annotation on the integrated data.
#   5.  Calculate and plot the proportions of each cell type across all samples.
#   6.  Check and compare the expression of AAV-C and AAV-N in different samples.
#   7.  Analyze and visualize the similarity of cell type compositions between
#       the different samples.
#
# =============================================================================

# Load packages and colorSet ----
for (p in c("tidyverse", "Seurat", "ggpubr", "harmony", "scCustomize", "stringr",
            "ggrepel", "clustree", "homologene", "dittoSeq", "ggalluvial", "scales")) {
  if (!require(p, character.only = T)) {
    stop(paste0("Please install ", p))
  } else {
    suppressMessages(library(p, quietly = T, character.only = T))
  }
}


source("~/project/EpiAllele/script/snRNAseq/snRNA_colorSet.R")

# Data integration ----
data.path <- "~/project/EpiAllele/result/snRNAseq/doubletfinder/"
result.path <- "~/project/EpiAllele/result/snRNAseq/integration/"
if (!dir.exists(result.path)) {
  dir.create(result.path)
}

## Quality control ----
data.list <- list()
samples <- c("WT", "R404Q", "NT", "SQ")
for (sample in samples) {
  data <- readRDS(paste0(data.path, sample, ".afterDoubeltFinder.rds"))
  
  data$dataset <- sample
  data$percent.mt <- PercentageFeatureSet(data, pattern = "^mt-")
  
  data <- subset(data, subset = nFeature_RNA > 200 & 
                   nFeature_RNA < quantile(nFeature_RNA, 0.99) & 
                   nCount_RNA < 20000 &
                   percent.mt < 20 &
                   DoubletFinder == "Singlet")
  
  data.list[[sample]] <- data
}

##cell cycle gene sets
s.genes <- homologene(cc.genes$s.genes, inTax = 9606, outTax = 10090) %>% 
  pull(`10090`)
g2m.genes <- homologene(cc.genes$g2m.genes, inTax = 9606, outTax = 10090) %>% 
  pull(`10090`)

## Check the batch among four samples ----
merged.data <- merge(data.list[[1]], y = data.list[2:4], add.cell.ids = names(data.list))
data.unintegrated <- merged.data %>% 
  NormalizeData() %>% 
  CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(vars.to.regress = c("percent.mt", "S.Score", "G2M.Score")) %>% 
  RunPCA(npcs = 50, verbose = FALSE) %>% 
  RunUMAP(reduction = "pca", dims = 1:20) %>% 
  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
  FindClusters(resolution = 1.0)

pdf(paste0(result.path,"unintegrated.umap.pdf"),width = 8,height = 6)
p1 <- DimPlot(data.unintegrated,reduction = "umap", group.by = "seurat_clusters", 
              label = T, label.size = 3)
print(p1)
p2 <- DimPlot(data.unintegrated,reduction = "umap", group.by = "dataset")+
  scale_color_manual(values=colorSet$dataset)
print(p2)
p <- FeaturePlot(data.unintegrated, features = "AAV-C")
print(p)
p <- FeaturePlot(data.unintegrated, features = "AAV-N")
print(p)
dev.off()

pdf(paste0(result.path,"unintegrated.bySample.umap.pdf"),width = 14,height = 10)
p1 <- DimPlot(data.unintegrated, reduction = "umap", group.by = "seurat_clusters",
              split.by = "dataset", ncol = 2)
print(p1)
p <- FeaturePlot(data.unintegrated, ncol = 2, features = "AAV-C", split.by = "dataset")+
  patchwork::plot_layout(ncol = 2, nrow = 2)
print(p)
p <- FeaturePlot(data.unintegrated, ncol = 2, features = "AAV-N", split.by = "dataset")+
  patchwork::plot_layout(ncol = 2, nrow = 2)
print(p)
dev.off()

saveRDS(data.unintegrated, paste0(result.path, "unintegrated.rds"))

## CCA ----
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
  CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) %>% 
  ScaleData(vars.to.regress = c("percent.mt", "S.Score", "G2M.Score")) %>% 
  RunPCA(npcs = 50, verbose = FALSE)

### Select optimal PCs for dimensional reduction
pdf(paste0(result.path,"CCA.integrated.ElbowPlot.pdf"),
    width = 7,height = 5)
ElbowPlot(data.integrated, ndims = 50)
dev.off()

### Clustering
data.integrated <- data.integrated %>% 
  RunUMAP(reduction = "pca", dims = 1:20) %>% 
  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
  FindClusters(resolution = 1.0)

### UMAP for visualization
DefaultAssay(data.integrated) <- "RNA"

pdf(paste0(result.path,"CCA.integrated.umap.pdf"),width = 8,height = 6)
p1 <- DimPlot(data.integrated,reduction = "umap", group.by = "seurat_clusters", 
              label = T, label.size = 3)
print(p1)
p2 <- DimPlot(data.integrated,reduction = "umap", group.by = "dataset")+
  scale_color_manual(values=colorSet$dataset)
print(p2)
p <- FeaturePlot(data.integrated, features = "AAV-C")
print(p)
p <- FeaturePlot(data.integrated, features = "AAV-N")
print(p)
dev.off()

pdf(paste0(result.path,"CCA.integrated.bySample.umap.pdf"),width = 14,height = 10)
p1 <- DimPlot(data.integrated, reduction = "umap", group.by = "seurat_clusters",
              split.by = "dataset", ncol = 2)
print(p1)
p <- FeaturePlot(data.integrated, ncol = 2, features = "AAV-C", split.by = "dataset")+
  patchwork::plot_layout(ncol = 2, nrow = 2)
print(p)
p <- FeaturePlot(data.integrated, ncol = 2, features = "AAV-N", split.by = "dataset")+
  patchwork::plot_layout(ncol = 2, nrow = 2)
print(p)
dev.off()

### store the data.integrated object
saveRDS(data.integrated, paste0(result.path, "CCA.integrated.rds"))

## Harmony ----
data.harmony <- merged.data %>% 
  NormalizeData() %>% 
  CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(vars.to.regress = c("percent.mt", "S.Score", "G2M.Score")) %>% 
  RunPCA(npcs = 50, verbose = FALSE) %>% 
  RunHarmony(group.by.vars = "dataset", theta=2)

pdf(paste0(result.path,"harmony.integrated.ElbowPlot.pdf"),
    width = 7,height = 5)
ElbowPlot(data.harmony, ndims = 50)
dev.off()

data.harmony <- data.harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 1.0)
Idents(data.harmony) <- "seurat_clusters"


pdf(paste0(result.path,"harmony.umap.pdf"),width = 8,height = 6)
p1 <- DimPlot(data.harmony,reduction = "umap", group.by = "seurat_clusters", 
              label = T, label.size = 3)
print(p1)
p2 <- DimPlot(data.harmony,reduction = "umap", group.by = "dataset")+
  scale_color_manual(values=colorSet$dataset)
print(p2)
p <- FeaturePlot(data.harmony, features = "AAV-C")
print(p)
p <- FeaturePlot(data.harmony, features = "AAV-N")
print(p)
dev.off()

pdf(paste0(result.path,"harmony.bySample.umap.pdf"),width = 14,height = 10)
p1 <- DimPlot(data.harmony, reduction = "umap", group.by = "seurat_clusters",
              split.by = "dataset", ncol = 2)
print(p1)
p <- FeaturePlot(data.harmony, ncol = 2, features = "AAV-C", split.by = "dataset")+
  patchwork::plot_layout(ncol = 2, nrow = 2)
print(p)
p <- FeaturePlot(data.harmony, ncol = 2, features = "AAV-N", split.by = "dataset")+
  patchwork::plot_layout(ncol = 2, nrow = 2)
print(p)
dev.off()
saveRDS(data.harmony, paste0(result.path, "harmony.rds"))

# FindAllMarkers ----
DefaultAssay(data.harmony) <- "RNA"
markers <- FindAllMarkers(data.harmony,logfc.threshold = 0)
markers$expression <- ifelse(markers$avg_log2FC>0, "up-regulated", "down-regulated")
rownames(markers) <- NULL
write.csv(markers, file = paste0(result.path, "harmony.markers.csv"))

DefaultAssay(data.integrated) <- "RNA"
markers <- FindAllMarkers(data.integrated,logfc.threshold = 0)
markers$expression <- ifelse(markers$avg_log2FC>0, "up-regulated", "down-regulated")
rownames(markers) <- NULL
write.csv(markers, file = paste0(result.path, "CCA.markers.csv"))



# Celltype annotation ----
## Check canonical marker genes expression ----
features1 <- list("Cardiomyocytes" = c("Ryr2", "Myh6", "Actc1", "Tpm", "Atp2a2", "Nppa", "Tnnc1", "Acta1", "Trdn"),
                 "Fibroblasts" = c("Col1a2", "Col3a1", "Vim", "Fstl1", "Gsn", "Fbln2", "Sparc", "Mmp2"),
                 "Endothelial cells" = c("Fabp4", "Tie1", "Egfl7", "Pecam1", "Flt1", "Epas", "Emcn", "Ednrb"),
                 "Macrophages" = c("Cd74", "Cd163", "Lyz1", "Csfr1", "Emr1", "Cd68", "Lgals3", "Itgam", "Trem2"))
features2 <- list("Schwann cells" = c("Plp1", "Kcna1", "Scn7a", "Gfra3", "Gpr37l1"),
                  "Cardiomyocytes" = c("Ttn", "Mhrt", "Myh6", "Tnnc1", "Pde4d"),
                  "Fibroblasts" = c("Pdgfra", "Col1a1", "Crispld2", "Ogn", "Islr"),
                  "Pericytes" = c("Pdgfrb", "Vtn", "Colec11", "Steap4", "Kcnj8", "Abcc9"),
                  "Lymphatic ECs" = c("Lyve1", "Mmrn1", "Ccl21a", "Cldn5", "Flt4"),
                  "Endothelial cells" = c("Pecam1", "Ly6c1", "Lims2", "Rgcc", "C1qtnf9", "Cyyr1"),
                  "Endocardial cells" = c("Npr3", "Cytl1", "Vwf", "Ptgs1", "Cgnl1", "Plvap"),
                  "Epicardial cells" = c("Msln", "Upk3b", "Gpc3", "Fmod"),
                  "Macrophages" = c("Fcgr1", "Csf1r", "Adgre1", "Pld4", "Ms4a6c", "Mgl2"),
                  "DCs" = c("Cd209a"), 
                  "Granulocytes" = c("S100a9", "Lmnb1", "Slpi", "Clec4d", "Retnlg", "Il1r2"),
                  "T cells" = c("Cd3e", "Hcst", "Gimap3", "Il7r", "Skap1"),
                  "NK cells" = c("Ncr1"),
                  "B cells" = c("Ms4a1", "Cd79a", "Ly6d", "H2-DMb2", "Cd79b"),
                  "Smooth muscle cells" = c("Acta2", "Tagln", "Myh11", "Olfr558", "Lmod1", "Nrip2"))

features3 <- list("Cardiomyocytes" = c("Trdn", "Nkx2-5", "Actn2", "Lmod2", "Myh6", "Lrrc10", "Ckmt2", "Trm63"),
                 "Endothelial cells" = c("Pcdh17", "Ldb2", "Mall", "She", "Ccdc85a", "Fam181b", "Gpr4"),
                 "Fibroblasts" = c("Lama2", "Twist1", "Ebf2", "Slc1a3", "Cdh11", "Tenm3", "Gng8", "Kcnk2", "Kcne4"),
                 "Macrophages" = c("Hpgds", "Slamf7", "Tbxas1", "Arhgap9", "Ly9", "Nxpe5", "Il7r", "Gpr65", "Kcnk13"))

features4 <- list("T cells" = c("Cd3e", "Cd3d", "Cd8a", "Trbc1", "Trgc1", "Foxp3"),
                  "NK cells" = c("Ncr1", "Nkg7", "Gzma", "Gzmb", "Klrb1c", "Klra4", "Prf1"),
                  "B cells" = c("Ms4a1", "Cd79a", "Cd79b", "Ighd", "Igha", "Igkc", "Jchain", "Ly6d", "Ly6g"),
                  "Schwann cells" = c("Plp1", "Kcna1", "Scn7a", "Gfra3", "Gpr37l1", "Mpz"),
                  "Neural-like cells" = c("Sox10", "S100b", "Fabp7", "Rbfox1", "Efemp1", "Pcdh15"),
                  "Fibroblasts" = c("Pdgfra", "Col1a1", "Col1a2", "Dkk3", "Tcf21"),
                  "Adipocyte" = c("Plin1", "Fabp4", "Acsl1", "Cidec", "Pparg", "Lpl"),
                  "Pericytes" = c("Pdgfrb", "Vtn", "Colec11", "Steap4", "Kcnj8", "Abcc9", "Rgs5", "Notch3", "Cspg4"),
                  "Smooth muscle cells" = c("Acta2", "Tagln", "Myh11", "Lmod1", "Nrip2", "Cnn1", "Smtn"),
                  "Endothelial cells" = c("Pecam1", "Ly6c1", "Lims2", "Rgcc", "C1qtnf9", "Cyyr1"),
                  "Endocardial cells" = c("Npr3", "Gja5", "Tbx3", "Cytl1", "Vwf", "Ptgs1", "Cgnl1", "Plvap"),
                  "Epicardial cells" = c("Msln", "Upk3b", "Wt1", "Tbx18", "Gpc3", "Fmod"),
                  "Granulocytes" = c("S100a9", "S100a8", "Lmnb1", "Mpo"),
                  "Progenitor/stem" = c("Cd34", "Kit", "Gata2", "Runx1"))

data.harmony <- readRDS(paste0(result.path, "harmony.rds"))

markers <- read.csv(paste0(result.path, "harmony.markers.csv"), header = T, 
                    check.names = F, row.names = 1)


pdf(paste0(result.path, "harmony.celltype.marker1.pdf"), width = 14, height = 9)
p <- DotPlot(data.harmony, features=features1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
print(p)
dev.off()

pdf(paste0(result.path, "harmony.celltype.marker2.pdf"), width = 30, height = 9)
p <- DotPlot(data.harmony, features=features2)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
print(p)
dev.off()

pdf(paste0(result.path, "harmony.celltype.marker3.pdf"), width = 16, height = 9)
p <- DotPlot(data.harmony, features=features3)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
print(p)
dev.off()

pdf(paste0(result.path, "harmony.celltype.marker4.pdf"), width = 30, height = 9)
p <- DotPlot(data.harmony, features=features4)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
print(p)
dev.off()

## Celltype annotation ----
annotation <- read.csv(paste0(result.path, "harmony_celltype_annotation.csv"), header=T)
cellTypes <- annotation$celltype
Idents(data.harmony) <- "seurat_clusters"
names(cellTypes) <- levels(data.harmony)
data.harmony <- RenameIdents(data.harmony, cellTypes)
data.harmony$cellType <- data.harmony@active.ident
data.harmony$cellType <- factor(data.harmony$cellType, 
                                levels = c("Cardiomyocytes", "Endothelial cells", "Endocardial cells",
                                           "Lymphatic endothelial cells", "Proliferative endothelial cells",
                                           "Fibroblasts", "Smooth muscle cells", "Adipocytes",
                                           "Macrophages", "Lymphocytes", "Pericytes", 
                                           "Schwann cells", "Neural-like cells"))

Idents(data.harmony) <- "cellType"
pdf(paste0(result.path,"harmony.celltype.umap.pdf"),width = 8,height = 6)
p <- DimPlot(data.harmony,reduction = "umap", group.by = "seurat_clusters", 
              label = T, label.size = 3)
print(p)
p <- DimPlot(data.harmony,reduction = "umap", group.by = "cellType")+
  scale_color_manual(values=colorSet$cellType)
print(p)
p <- DimPlot(data.harmony,reduction = "umap", group.by = "dataset")+
  scale_color_manual(values=colorSet$dataset)
print(p)
p <- FeaturePlot(data.harmony, features = "AAV-C",
                cols = c("gray90", "red"), pt.size=1.0)
print(p)
p <- FeaturePlot(data.harmony, features = "AAV-N",
                cols = c("gray90", "red"), pt.size=1.0)
print(p)
dev.off()

pdf(paste0(result.path,"harmony.celltype.bySample.umap.pdf"),width = 14,height = 10)
p1 <- DimPlot(data.harmony, reduction = "umap", group.by = "seurat_clusters",
              split.by = "dataset", ncol = 2, pt.size = 0.6)
print(p1)
p <- DimPlot(data.harmony,reduction = "umap", group.by = "cellType", 
             split.by = "dataset", ncol = 2, pt.size = 0.6)+
  scale_color_manual(values=colorSet$cellType)
print(p)
p <- FeaturePlot_scCustom(data.harmony, num_columns = 2, features = "AAV-C", split.by = "dataset",
                          colors_use  = c("gray90", "red"), pt.size=1.0)
print(p)
p <- FeaturePlot_scCustom(data.harmony, num_columns = 2, features = "AAV-N", split.by = "dataset",
                          colors_use  = c("gray90", "red"), pt.size=1.0)
print(p)
dev.off()


## top5 makers for each seurat_cluster ----
markers_top5 <- markers %>% 
  filter(expression == "up-regulated") %>% 
  group_by(cluster) %>% 
  do(head(., n = 5))

Idents(data.harmony) <- "seurat_clusters"
pdf(paste0(result.path, "harmony.cluster.top5marker.Dotplot.pdf"), width = 35, height = 18)
p <- DotPlot(data.harmony,assay = "RNA",features = unique(markers_top5$gene),
             cols = c(low="blue",high = "red"))+
  theme(axis.text.x = element_text(angle = 45,size = 18,vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 24),
        axis.title = element_text(face = "bold",size=30),
        legend.text = element_text(size = 18),
        legend.title=element_text(face = "bold",size = 20))
print(p)
dev.off()

## top5 makers for each celltype ----
DefaultAssay(data.harmony) <- "RNA"
celltype_markers <- FindAllMarkers(data.harmony, group.by="cellType", logfc.threshold = 0)
celltype_markers$expression <- ifelse(celltype_markers$avg_log2FC>0, 
                                      "up-regulated", "down-regulated")
celltype_markers_top5 <- celltype_markers %>% 
  filter(expression == "up-regulated") %>% 
  group_by(cluster) %>% 
  do(head(., n = 5))

Idents(data.harmony) <- "cellType"
pdf(paste0(result.path, "harmony.celltype.top5marker.Dotplot.pdf"), width = 35, height = 18)
p <- DotPlot(data.harmony,assay = "RNA",features = unique(celltype_markers_top5$gene),
             cols = c(low="blue",high = "red"))+
  theme(axis.text.x = element_text(angle = 45,size = 18,vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 24),
        axis.title = element_text(face = "bold",size=30),
        legend.text = element_text(size = 18),
        legend.title=element_text(face = "bold",size = 20))
print(p)
dev.off()

saveRDS(data.harmony, paste0(result.path, "harmony.annotation.rds"))


# Check AAV expression ----
pdf(paste0(result.path, "harmony.AAV.check.Dotplot.pdf"), width = 10, height = 4)
p <- dittoDotPlot(data.harmony, vars = c("AAV-C", "AAV-N"), group.by = "cellType",
                  split.by = "dataset", split.ncol = 4)
print(p)
dev.off()

data.harmony$AAV_C <- data.harmony@assays$RNA@counts[c("AAV-C"),]
data.harmony$AAV_N <- data.harmony@assays$RNA@counts[c("AAV-N"),]
data.harmony$AAV_C <- ifelse(data.harmony@assays$RNA@counts[c("AAV-C"),]>0, 
                             "AAV_C express", "AAV_C notexpress")
data.harmony$AAV_N <- ifelse(data.harmony@assays$RNA@counts[c("AAV-N"),]>0, 
                             "AAV_N express", "AAV_N notexpress")
AAV_stat <- table(data.harmony$cellType, data.harmony$dataset, 
                  data.harmony$AAV_C, data.harmony$AAV_N) %>% 
  as.data.frame()
colnames(AAV_stat) <- c("cellType", "dataset", "AAV_C", "AAV_N", "Freq")
AAV_stat <- AAV_stat %>% 
  group_by(dataset, cellType) %>% 
  mutate(fraction=Freq/sum(Freq), percentage=percent(Freq/sum(Freq)), 
         AAV=paste0(AAV_C, ", ", AAV_N))
cardio_stat <- subset(AAV_stat, cellType=="Cardiomyocytes")

pdf(paste0(result.path, "AAV_stat.pdf"), width = 12, height = 8)
p <- ggplot(AAV_stat, aes(dataset, fraction, fill=AAV))+
  geom_bar(stat = "identity", position = "stack")+
  theme_bw()+
  facet_wrap(~cellType)
print(p)
dev.off()

pdf(paste0(result.path, "cardio_stat.pdf"), width = 8, height = 6)
p <- ggplot(cardio_stat, aes(dataset, fraction, fill=AAV))+
  geom_bar(stat = "identity", position = "stack")+
  theme_bw()+
  facet_wrap(~cellType)
print(p)
dev.off()

write.csv(AAV_stat, paste0(result.path, "AAV_stat.csv"))

# Celltype proportion statistics ----
metadata <- data.harmony@meta.data
result.path <- "~/result/snRNAseq/proportion_stat/"
if (!dir.exists(result.path)) {
  dir.create(result.path)
}

celltype_stat <- table(metadata$dataset, metadata$cellType) %>% melt()
colnames(celltype_stat) <- c("Dataset", "CellType", "Number")

## Calculate celltype proportion ----
celltype_stat <- celltype_stat %>%
  group_by(Dataset) %>%
  mutate(Proportion = Number / sum(Number))
celltype_stat$Dataset <- factor(celltype_stat$Dataset, levels = c("NT", "R404Q", "SQ", "WT"))

## Generate bar plot and alluvial plot ----
pdf(paste0(result.path, "stack_barplot.pdf"),width = 9,height = 6)
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
  scale_fill_manual(values=colorSet$cellType)
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
  scale_fill_manual(values=colorSet$cellType)
print(p)
dev.off()

## Correlation of celltype proportion between groups ----
prop_matrix <- celltype_stat %>% 
  select(Dataset, CellType, Proportion) %>% 
  pivot_wider(names_from = Dataset, values_from = Proportion) %>%
  column_to_rownames("CellType")

cor_matrix <- cor(prop_matrix, method = "spearman")

pdf(paste0(result.path, "proportion.cor.heatmap.pdf"), width = 4, height = 3)
p1 <- pheatmap(cor_matrix, show_rownames = T,
               cluster_rows = F, cluster_cols = F,
               color = colorRampPalette(c("white", "#D5676D"))(100),
               breaks = seq(0.5, 1, length.out = 100),
               display_numbers = T, number_color = "white",
               heatmap_legend_param = list(
                 title = "correlation"
               ))
draw(p1)
dev.off()





