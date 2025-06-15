#!/usr/bin/env Rscript
# quality check ----
runQC <- function(sample, data.path, result.path) {
  seu <- CreateSeuratObject(Read10X(paste0(data.path, sample, "/outs/filtered_feature_bc_matrix/") , 
                                    gene.column = 2), min.cells = 1)
  seu$percent.mt <- PercentageFeatureSet(seu, pattern = "^mt-")
  pdf(paste0(result.path, sample, ".qc.pdf"), width = 8, height = 4)
  p <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(p)
  dev.off()
  
  ## remove low quality cells ----
  seu <- subset(seu, subset = nFeature_RNA > 200 & 
                  nFeature_RNA < quantile(nFeature_RNA, 0.99) & 
                  nCount_RNA < 20000 &
                  percent.mt < 20)
  
  s.genes <- homologene(cc.genes$s.genes, inTax = 9606, outTax = 10090) %>% 
    pull(`10090`)
  g2m.genes <- homologene(cc.genes$g2m.genes, inTax = 9606, outTax = 10090) %>% 
    pull(`10090`)
  seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  seu <- seu %>% 
    SCTransform(vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), verbose = FALSE) %>% 
    RunPCA(verbose = FALSE) %>% 
    RunUMAP(dims = 1:30, verbose = FALSE) %>% 
    FindNeighbors(dims = 1:30, verbose = FALSE) %>% 
    FindClusters(verbose = FALSE)
  
  ## check ambient RNA: IG and HB genes ----
  DefaultAssay(seu) <- "RNA"
  IG <- rownames(seu)[grepl(pattern = paste(c("Jchain","Igkc.*","Igha.*","Iglc.*","Igkv.*","Ighv.*","Iglv.*"), collapse = "|"),rownames(seu))]
  HB <- rownames(seu)[grepl(pattern = paste(c("Hbb-.*","Hba-.*"), collapse = "|"),rownames(seu))]
  if (length(IG) > 0) {
    IG_score <- AddModuleScore(seu, features = IG, name = "IG")
  }
  if (length(HB) > 0) {
    HB_score <- AddModuleScore(seu, features = HB, name = "HB")
  }
  
  pdf(paste0(result.path, sample, ".umap.pdf"), width = 8, height = 6)
  p <- DimPlot(seu, label = TRUE)
  print(p)
  p <- FeaturePlot(seu, features = "percent.mt")
  print(p)
  p <- FeaturePlot(IG_score, features = "IG1")
  print(p)
  p <- FeaturePlot(HB_score, features = "HB1")
  print(p)
  dev.off()
  
  saveRDS(seu, paste0(result.path, sample, ".raw.rds"))
}

# Soupx to remove ambient RNA ----
RunSoupX <- function(sample, data.path, result.path){
  raw_matrix <- paste0(data.path, sample, "/outs/raw_feature_bc_matrix/")
  filter_matrix <- paste0(data.path, sample, "/outs/filtered_feature_bc_matrix/")
  tod <- Read10X(raw_matrix, gene.column = 2)
  toc <- Read10X(filter_matrix, gene.column = 2)
  tod <- tod[rownames(toc),]
  data <- CreateSeuratObject(toc)
  
  ## Get the genelists for checking ----
  IG <- rownames(toc)[grepl(pattern = paste(c("Jchain","Igkc.*","Igha.*","Iglc.*","Igkv.*","Ighv.*","Iglv.*"), collapse = "|"),rownames(toc))]
  HB <- rownames(toc)[grepl(pattern = paste(c("Hbb-.*","Hba-.*"), collapse = "|"),rownames(toc))]
  
  ## Run Seurat to get cluster information ----
  data <- data %>% 
    NormalizeData(normalization.method = "LogNormalize") %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA() %>% 
    FindNeighbors(dims = 1:30) %>% 
    FindClusters() %>% 
    RunUMAP(dims = 1:30)
  

  ## Adding extra meta data to the SoupChannel object ----
  sc = SoupChannel(tod, toc)
  matx <- data@meta.data
  sc = setClusters(sc,setNames(matx$seurat_clusters, rownames(matx)))
  
  ## Background noise correction ----
  umap = data@reductions$umap@cell.embeddings
  sc <- setDR(sc, umap)
  sc = setContaminationFraction(sc, 0.2)
  out = adjustCounts(sc)
  seurat_obj <- CreateSeuratObject(out)
  
  ## Plot change map to check the correction ----
  p1 <- plotChangeMap(sc, out, geneSet = IG)+
    labs(title = "Change in Immunoglobin genes expression")
  p2 <- plotChangeMap(sc, out, geneSet = HB)+
    labs(title = "Change in Hemoglobin genes expression")
  
  pdf(paste0(result.path, sample, ".changeplot.pdf"), width = 8, height = 6)
  print(p1)
  print(p2)
  dev.off()
  
  saveRDS(seurat_obj, paste0(result.path, sample, ".afterSoupX.rds"))
}


# DoubletFinder to detect and remove doublets ----
RunDoubletFinder <- function(sample, data.path, result.path){
  ## Seurat standard workflow ----
  data <- readRDS(paste0(data.path, sample, ".afterSoupX.rds"))
  
  data <- data %>% 
    NormalizeData(normalization.method = "LogNormalize") %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA() %>% 
    FindNeighbors(dims = 1:30) %>% 
    FindClusters() %>% 
    RunUMAP(dims = 1:30)
  
  ## pK Identification ----
  sweep.res.list <- paramSweep(data, PCs = 1:30, sct = F)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
  pdf(paste0(result.path, sample, ".Find_pk.pdf"),width = 10,height = 8)
  bcmvn <- find.pK(sweep.stats)
  dev.off()
  pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])
  
  ## Homotypic Doublet Proportion Estimate ----
  annotations <- data@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp_poi <- round(0.075*nrow(data@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----
  seurat_obj <- doubletFinder(data, PCs = 1:10, pN = 0.25, pK = pK_bcmvn, 
                              nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  seurat_obj$DoubletFinder <- seurat_obj[[paste0('DF.classifications_0.25_',
                                                 pK_bcmvn,'_',nExp_poi.adj)]]
  
  pdf(paste0(result.path, sample, ".doubletFinder.pdf"), width = 8, height = 6)
  p <- DimPlot(seurat_obj, group.by = "DoubletFinder")
  print(p)
  dev.off()
  
  meta <- seurat_obj@meta.data %>% 
    select(orig.ident, nCount_RNA, nFeature_RNA, DoubletFinder)
  seurat_obj <- CreateSeuratObject(counts = seurat_obj@assays$RNA@counts,
                                   meta.data = meta, project = sample)
  
  saveRDS(seurat_obj, paste0(result.path, sample, ".afterDoubeltFinder.rds"))
}

# remove doublets and perform standard cluster steps ----
RunCluster <- function(sample, data.path, result.path, 
                       checkAAV = FALSE, Findmarker = FALSE){
  ## Seurat standard workflow ----
  data <- readRDS(paste0(data.path, sample, ".afterDoubeltFinder.rds"))
  data$percent.mt <- PercentageFeatureSet(data, pattern = "^mt-")
  s.genes <- homologene(cc.genes$s.genes, inTax = 9606, outTax = 10090) %>% 
    pull(`10090`)
  g2m.genes <- homologene(cc.genes$g2m.genes, inTax = 9606, outTax = 10090) %>% 
    pull(`10090`)
  
  data <- subset(data, subset = nFeature_RNA > 200 & 
                   nFeature_RNA < quantile(nFeature_RNA, 0.99) & 
                   nCount_RNA < 20000 &
                   percent.mt < 20 &
                   DoubletFinder == "Singlet")
  
  data <- data %>% 
    NormalizeData(normalization.method = "LogNormalize") %>% 
    CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) %>% 
    FindVariableFeatures() %>% 
    ScaleData(vars.to.regress = c("percent.mt", "S.Score", "G2M.Score")) %>% 
    RunPCA() %>% 
    FindNeighbors(dims = 1:30) %>% 
    FindClusters() %>% 
    RunUMAP(dims = 1:30)
  
  pdf(paste0(result.path, sample, ".cluster.umap.pdf"), width = 8, height = 6)
  p <- DimPlot(data, label = TRUE)
  print(p)
  if (checkAAV) {
    p <- FeaturePlot(data, features = "AAV-C")
    print(p)
    p <- FeaturePlot(data, features = "AAV-N")
    print(p)
  }
  dev.off()
  
  ## find markers ----
  if (Findmarker) {
    DefaultAssay(data) <- "RNA"
    markers <- FindAllMarkers(data,logfc.threshold = 0)
    markers$expression <- ifelse(markers$avg_log2FC>0,"up-regulated","down-regulated")
    rownames(markers) <- NULL
    write.csv(markers,file = paste0(result.path, sample, ".markers.csv"))
  }
  
  saveRDS(data, paste0(result.path, sample, ".normalized.rds"))
}


