#!/usr/bin/env Rscript

## {{{ Install required libraries if required
for (p in c("tidyverse", "stringr", "ggrepel", "ggsci", "ggVennDiagram", "ggvenn")) {
  if (!require(p, character.only = T)) {
    stop(paste0("Please install ", p))
  } else {
    suppressMessages(library(p, quietly = T, character.only = T))
  }
}
## }}}

# intersect: D2/D5/snRNA cardiomyocytes DEG ----
cardio_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/snRNAseq/DE/Cardiomyocytes.SQ.vs.NT.DEresult.csv",
                      header=T, row.names=1)
D2_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/20250211_mouse_RNAseq/DE/D2sg.vs.D2NT.DEresult.csv",
                    header=T)
D5_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/20250211_mouse_RNAseq/DE/D5sg.vs.D5NT.DEresult.csv",
                    header=T)

up.list <- list("D2.bulk.up" = D2_DEG %>% filter(group=="UP") %>% pull(gene),
                "D5.bulk.up" = D5_DEG %>% filter(group=="UP") %>% pull(gene),
                "snRNA.cardiomyocytes.up" = cardio_DEG %>% filter(group=="UP") %>% pull(gene))
down.list <- list("D2.bulk.down" = D2_DEG %>% filter(group=="DOWN") %>% pull(gene),
                "D5.bulk.down" = D5_DEG %>% filter(group=="DOWN") %>% pull(gene),
                "snRNA.cardiomyocytes.down" = cardio_DEG %>% filter(group=="DOWN") %>% pull(gene))
pdf("/data02/hukaijie/EpiAllele/result/combine/sg.vs.nt.rnaseq.vennplot.pdf",
      width = 5, height = 4)
p1 <- ggvenn(
  up.list, fill_alpha = 0.5, stroke_color = "white",
  fill_color = c("#8FBBD9", "#F4A582", "#92C47D"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
)
print(p1)
p2 <- ggvenn(
  down.list, fill_alpha = 0.5, stroke_color = "white",
  fill_color = c("#8FBBD9", "#F4A582", "#92C47D"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
)
print(p2)
dev.off()

rna_up_intersect <- intersect(up.list$`D2.bulk.up`, up.list$`D5.bulk.up`)
rna_down_intersect <- intersect(down.list$`D2.bulk.down`, down.list$`D5.bulk.down`)
write.table(data.frame(Gene = rna_up_intersect), 
            file = "/data02/hukaijie/EpiAllele/result/combine/d2.vs.d5.rnaseq.up_intersect.txt",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

write.table(data.frame(Gene = rna_down_intersect), 
            file = "/data02/hukaijie/EpiAllele/result/combine/d2.vs.d5.rnaseq.down_intersect.txt",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


# intersect: D5 RNAseq + chipseq + atacseq ----
D5_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/20250211_mouse_RNAseq/DE/D5sg.vs.D5NT.DEresult.csv",
                    header=T)
GATA4_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/20250304_mouse_chipseq/DE/GATA4.sg.vs.NT.peak.DEresult.csv",
                      header=T)
Nkx2_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/20250304_mouse_chipseq/DE/Nkx2-5.sg.vs.NT.peak.DEresult.csv",
                    header=T)
D5_ATAC <- read.csv("/data02/hukaijie/EpiAllele/result/20250211_result_cuttag/combine_atac/DE/ATAC.sg.vs.NT.peak.DEresult.csv",
                    header=T)
up.list <- list("GATA4.up" = GATA4_DEG %>% filter(group=="UP") %>% pull(SYMBOL) %>% unique(.),
                "Nkx2.5.up" = Nkx2_DEG %>% filter(group=="UP") %>% pull(SYMBOL) %>% unique(.),
                "D5.ATAC.up" = D5_ATAC %>% filter(group== "UP") %>% pull(SYMBOL) %>% unique(.),
                "D5.RNAseq.up" = D5_DEG %>% filter(group=="UP") %>% pull(gene))
down.list <- list("GATA4.down" = GATA4_DEG %>% filter(group=="DOWN") %>% pull(SYMBOL) %>% unique(.),
                "Nkx2.5.down" = Nkx2_DEG %>% filter(group=="DOWN") %>% pull(SYMBOL) %>% unique(.),
                "D5.ATAC.down" = D5_ATAC %>% filter(group== "DOWN") %>% pull(SYMBOL) %>% unique(.),
                "D5.RNAseq.down" = D5_DEG %>% filter(group=="DOWN") %>% pull(gene))
total.list <- list("GATA4" = GATA4_DEG %>% filter(group != "NOT-SIG") %>% pull(SYMBOL) %>% unique(.),
                "Nkx2.5" = Nkx2_DEG %>% filter(group != "NOT-SIG") %>% pull(SYMBOL) %>% unique(.),
                "D5.ATAC" = D5_ATAC %>% filter(group != "NOT-SIG") %>% pull(SYMBOL) %>% unique(.),
                "D5.RNAseq" = D5_DEG %>% filter(group != "NOT-SIG") %>% pull(gene))
pdf("/data02/hukaijie/EpiAllele/result/combine/sg.vs.nt.tf.vs.rna.vennplot.pdf",
      width = 5, height = 4)
p1 <- ggvenn(
  up.list,  fill_alpha = 0.5, stroke_color = "white",
  fill_color = c("#EB6A6C", "#6FAEDB", "#73C270", "#B97FC9"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
)
print(p1)
p2 <- ggvenn(
  down.list,  fill_alpha = 0.5, stroke_color = "white",
  fill_color = c("#EB6A6C", "#6FAEDB", "#73C270", "#B97FC9"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
)
print(p2)
p3 <- ggvenn(
  total.list,  fill_alpha = 0.5, stroke_color = "white",
  fill_color =  c("#EB6A6C", "#6FAEDB", "#73C270", "#B97FC9"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
)
print(p3)
dev.off()

## D2: D2 RNA-seq+ D2 WGBS + D2 histone cuttag ----
D2_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/20250211_mouse_RNAseq/DE/D2sg.vs.D2NT.DEresult.csv",
                    header=T)
H3K9me3_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/202407_mouse_d2_cuttag_myh6/combine_histone/DE/H3K9me3.sg.vs.NT.peak.DEresult.csv",
                      header=T)
H3K4me3_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/202407_mouse_d2_cuttag_myh6/combine_histone/DE/H3K4me3.sg.vs.NT.peak.DEresult.csv",
                    header=T)
H3K27me3_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/202407_mouse_d2_cuttag_myh6/combine_histone/DE/H3K27me3.sg.vs.NT.peak.DEresult.csv",
                    header=T)
H3K27ac_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/202407_mouse_d2_cuttag_myh6/combine_histone/DE/H3K27ac.sg.vs.NT.peak.DEresult.csv",
                    header=T)
H3K9ac_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/202407_mouse_d2_cuttag_myh6/combine_histone/DE/H3K9ac.sg.vs.NT.peak.DEresult.csv",
                    header=T)
WGBS_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/20250310_wgbs/results/anno_dmrs.csv", 
                    header=T, row.names = NULL, check.names = TRUE)
total.list <- list("RNA" = D2_DEG %>% filter(group != "NOT-SIG") %>% pull(gene) %>% unique(.),
                  "DNA methylation" =  WGBS_DEG %>% pull(SYMBOL) %>% unique(.),
                  "H3K27ac" = H3K27ac_DEG %>% filter(group != "NOT-SIG") %>% pull(SYMBOL) %>% unique(.),
                  "H3K27me3" = H3K27me3_DEG %>% filter(group != "NOT-SIG") %>% pull(SYMBOL) %>% unique(.),
                  "H3K9ac" = H3K9ac_DEG %>% filter(group != "NOT-SIG") %>% pull(SYMBOL) %>% unique(.),
                  "H3K9me3" = H3K9me3_DEG %>% filter(group != "NOT-SIG") %>% pull(SYMBOL) %>% unique(.),
                  "H3K4me3" = H3K4me3_DEG %>% filter(group != "NOT-SIG") %>% pull(SYMBOL) %>% unique(.)) 

pdf("/data02/hukaijie/EpiAllele/result/combine/sg.vs.nt.d2.combine.vennplot.pdf",
      width = 8, height = 8)
p <- ggVennDiagram(
  total.list,
  label_alpha = 0,
  label = "count", 
  set_color = c("#028BB8", "#A55A41", "#6D9F4B", "#692B78", "#9B2951", "#1D4B6F", "#7E65A2")
) +
  scale_fill_gradient(low = "white", high = "#993333") + 
  theme(legend.position = "none")

print(p)
dev.off()

# day5 ----
D5_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/20250211_mouse_RNAseq/DE/D5sg.vs.D5NT.DEresult.csv",
                    header=T)
H3K9me3_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/20250211_result_cuttag/combine_histone/DE/H3K9me3.sg.vs.NT.peak.DEresult.csv",
                      header=T)
H3K4me3_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/20250211_result_cuttag/combine_histone/DE/H3K4me3.sg.vs.NT.peak.DEresult.csv",
                    header=T)
H3K27me3_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/20250211_result_cuttag/combine_histone/DE/H3K27me3.sg.vs.NT.peak.DEresult.csv",
                    header=T)
H3K27ac_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/20250211_result_cuttag/combine_histone/DE/H3K27ac.sg.vs.NT.peak.DEresult.csv",
                    header=T)
H3K9ac_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/20250211_result_cuttag/combine_histone/DE/H3K9ac.sg.vs.NT.peak.DEresult.csv",
                    header=T)
WGBS_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/20250212_wgbs/results/anno_dmrs.csv", 
                    header=T, row.names = NULL, check.names = TRUE)
total.list <- list("RNA" = D5_DEG %>% filter(group != "NOT-SIG") %>% pull(gene) %>% unique(.),
                  "DNA methylation" =  WGBS_DEG %>% pull(SYMBOL) %>% unique(.),
                  "H3K27ac" = H3K27ac_DEG %>% filter(group != "NOT-SIG") %>% pull(SYMBOL) %>% unique(.),
                  "H3K27me3" = H3K27me3_DEG %>% filter(group != "NOT-SIG") %>% pull(SYMBOL) %>% unique(.),
                  "H3K9ac" = H3K9ac_DEG %>% filter(group != "NOT-SIG") %>% pull(SYMBOL) %>% unique(.),
                  "H3K9me3" = H3K9me3_DEG %>% filter(group != "NOT-SIG") %>% pull(SYMBOL) %>% unique(.),
                  "H3K4me3" = H3K4me3_DEG %>% filter(group != "NOT-SIG") %>% pull(SYMBOL) %>% unique(.)) 

pdf("/data02/hukaijie/EpiAllele/result/combine/sg.vs.nt.d5.combine.vennplot.pdf",
      width = 8, height = 8)
p <- ggVennDiagram(
  total.list,
  label_alpha = 0,
  label = "count", 
  set_color = c("#028BB8", "#A55A41", "#6D9F4B", "#692B78", "#9B2951", "#1D4B6F", "#7E65A2")
) +
  scale_fill_gradient(low = "white", high = "#993333") + 
  theme(legend.position = "none")

print(p)
dev.off()

# invivo ----
invivo_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/snRNAseq/DE/Cardiomyocytes.SQ.vs.NT.DEresult.csv",
                    header=T)
H3K9me3_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/202405_mouse_invivo_cuttag/combine_histone/DE/H3K9me3.sg.vs.NT.peak.DEresult.csv",
                      header=T)
H3K4me3_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/202405_mouse_invivo_cuttag/combine_histone/DE/H3K4me3.sg.vs.NT.peak.DEresult.csv",
                    header=T)
H3K27me3_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/202405_mouse_invivo_cuttag/combine_histone/DE/H3K27me3.sg.vs.NT.peak.DEresult.csv",
                    header=T)
H3K27ac_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/202405_mouse_invivo_cuttag/combine_histone/DE/H3K27ac.sg.vs.NT.peak.DEresult.csv",
                    header=T)
H3K9ac_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/202405_mouse_invivo_cuttag/combine_histone/DE/H3K9ac.sg.vs.NT.peak.DEresult.csv",
                    header=T)
WGBS_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/20240521_wgbs/results/sg.vs.NT.anno_dmrs.csv", 
                    header=T, row.names = NULL, check.names = TRUE)
total.list <- list("RNA" = invivo_DEG %>% filter(group != "NOT-SIG") %>% pull(gene) %>% unique(.),
                  "DNA methylation" =  WGBS_DEG %>% pull(SYMBOL) %>% unique(.),
                  "H3K27ac" = H3K27ac_DEG %>% filter(group != "NOT-SIG") %>% pull(SYMBOL) %>% unique(.),
                  "H3K27me3" = H3K27me3_DEG %>% filter(group != "NOT-SIG") %>% pull(SYMBOL) %>% unique(.),
                  "H3K9ac" = H3K9ac_DEG %>% filter(group != "NOT-SIG") %>% pull(SYMBOL) %>% unique(.),
                  "H3K9me3" = H3K9me3_DEG %>% filter(group != "NOT-SIG") %>% pull(SYMBOL) %>% unique(.),
                  "H3K4me3" = H3K4me3_DEG %>% filter(group != "NOT-SIG") %>% pull(SYMBOL) %>% unique(.)) 

pdf("/data02/hukaijie/EpiAllele/result/combine/sg.vs.nt.invivo.combine.vennplot.pdf",
      width = 8, height = 8)
p <- ggVennDiagram(
  total.list,
  label_alpha = 0,
  label = "count", 
  set_color = c("#028BB8", "#A55A41", "#6D9F4B", "#692B78", "#9B2951", "#1D4B6F", "#7E65A2")
) +
  scale_fill_gradient(low = "white", high = "#993333") + 
  theme(legend.position = "none")

print(p)
dev.off()

# snRNAseq ----
## DE for WT VS NT ----
result.path <- "/data02/hukaijie/EpiAllele/result/snRNAseq/DE/"
data <- readRDS("/data02/hukaijie/EpiAllele/result/snRNAseq/Rdata/harmony.annotation.rds")
celltype="Cardiomyocytes"
subset.data <- subset(data, cellType == celltype)
Idents(subset.data) <- "dataset"
DefaultAssay(subset.data) <- "RNA"
group1 <- "WT"
group2 <- "NT"
group_name <- paste0(group1, ".vs.", group2)
#Perform DEG analysis
DEGs <- FindMarkers(subset.data, ident.1 = group1, ident.2 = group2, 
                    logfc.threshold = 0)
DEGs$gene <- rownames(DEGs)
DEGs <- DEGs %>% 
  mutate(group = ifelse((p_val_adj < 0.05 & avg_log2FC > 0.5), "UP",
                        ifelse((p_val_adj < 0.05 & avg_log2FC < -0.5), "DOWN", "NOT-SIG"))) %>% 
  mutate(group_new = ifelse(gene == "Myh6", "Myh6", group)) %>% 
  arrange(group_new == "Myh6")
DEGs$celltype <- celltype
row.names(DEGs) <- NULL
write.csv(DEGs, 
          paste0(result.path, str_replace(celltype, " ", "."), ".", group_name, ".DEresult.csv"))


## Vennplot ----
DE_results <- read.csv("/data02/hukaijie/EpiAllele/result/snRNAseq/DE/Cardiomyocytes.DEG.csv", header=T, row.names=1)
SG_NT <- read.csv("/data02/hukaijie/EpiAllele/result/snRNAseq/DE/Cardiomyocytes.SQ.vs.NT.DEresult.csv", header=T, row.names=1)
SG_R404Q <- read.csv("/data02/hukaijie/EpiAllele/result/snRNAseq/DE/Cardiomyocytes.SQ.vs.R404Q.DEresult.csv", header=T, row.names=1)
WT_R404Q <- read.csv("/data02/hukaijie/EpiAllele/result/snRNAseq/DE/Cardiomyocytes.WT.vs.R404Q.DEresult.csv", header=T, row.names=1)
WT_NT <- read.csv("/data02/hukaijie/EpiAllele/result/snRNAseq/DE/Cardiomyocytes.WT.vs.NT.DEresult.csv", header=T, row.names=1)
#SG VS R404Q/ WT VS R404Q
up.list <- list("SG.vs.R404Q.up" = SG_R404Q %>% filter(group=="UP") %>% pull(gene),
                "WT.vs.R404Q.up" = WT_R404Q %>% filter(group=="UP") %>% pull(gene))
down.list <- list("SG.vs.R404Q.down" = SG_R404Q %>% filter(group=="DOWN") %>% pull(gene),
                "WT.vs.R404Q.down" = WT_R404Q %>% filter(group=="DOWN") %>% pull(gene))
total.list <- list("SG.vs.R404Q.diff" = SG_R404Q %>% filter(group != "NOT-SIG") %>% pull(gene),
                "WT.vs.R404Q.diff" = WT_R404Q %>% filter(group != "NOT-SIG") %>% pull(gene))

pdf("/data02/hukaijie/EpiAllele/result/combine/snRNA.comp1.vennplot.pdf",
      width = 5, height = 4)
p1 <- ggvenn(
  up.list, fill_alpha = 0.5, stroke_color = "white",
  fill_color = c("#E78AC3", "#A6D854"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
)
print(p1)
p2 <- ggvenn(
  down.list, fill_alpha = 0.5, stroke_color = "white",
  fill_color = c("#E78AC3", "#A6D854"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
)
print(p2)
p3 <- ggvenn(
  total.list, fill_alpha = 0.5, stroke_color = "white",
  fill_color = c("#E78AC3", "#A6D854"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
)
print(p3)
dev.off()

#SG VS NT/ WT VS NT
up.list <- list("SG.vs.NT.up" = SG_NT %>% filter(group=="UP") %>% pull(gene),
                "WT.vs.NT.up" = WT_NT %>% filter(group=="UP") %>% pull(gene))
down.list <- list("SG.vs.NT.down" = SG_NT %>% filter(group=="DOWN") %>% pull(gene),
                "WT.vs.NT.down" = WT_NT %>% filter(group=="DOWN") %>% pull(gene))
total.list <- list("SG.vs.NT.diff" = SG_NT %>% filter(group != "NOT-SIG") %>% pull(gene),
                "WT.vs.NT.diff" = WT_NT %>% filter(group != "NOT-SIG") %>% pull(gene))

pdf("/data02/hukaijie/EpiAllele/result/combine/snRNA.comp2.vennplot.pdf",
      width = 5, height = 4)
p1 <- ggvenn(
  up.list, fill_alpha = 0.5, stroke_color = "white",
  fill_color = c("#AB82FF", "#66C2A5"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
)
print(p1)
p2 <- ggvenn(
  down.list, fill_alpha = 0.5, stroke_color = "white",
  fill_color = c("#AB82FF", "#66C2A5"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
)
print(p2)
p3 <- ggvenn(
  total.list, fill_alpha = 0.5, stroke_color = "white",
  fill_color = c("#AB82FF", "#66C2A5"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
)
print(p3)
dev.off()

#SG VS NT/SG VS R404Q/ WT VS R404Q
up.list <- list("SG.vs.NT.up" = SG_NT %>% filter(group=="UP") %>% pull(gene),
                "SG.vs.R404Q.up" = SG_R404Q %>% filter(group=="UP") %>% pull(gene),
                "WT.vs.R404Q.up" = WT_R404Q %>% filter(group=="UP") %>% pull(gene))
down.list <- list("SG.vs.NT.down" = SG_NT %>% filter(group=="DOWN") %>% pull(gene),
                "SG.vs.R404Q.down" = SG_R404Q %>% filter(group=="DOWN") %>% pull(gene),
                "WT.vs.R404Q.down" = WT_R404Q %>% filter(group=="DOWN") %>% pull(gene))
total.list <- list("SG.vs.NT.diff" = SG_NT %>% filter(group != "NOT-SIG") %>% pull(gene),
                "SG.vs.R404Q.diff" = SG_R404Q %>% filter(group != "NOT-SIG") %>% pull(gene),
                "WT.vs.R404Q.diff" = WT_R404Q %>% filter(group != "NOT-SIG") %>% pull(gene))

pdf("/data02/hukaijie/EpiAllele/result/combine/snRNA.comp3.vennplot.pdf",
      width = 5, height = 4)
p1 <- ggvenn(
  up.list, fill_alpha = 0.5, stroke_color = "white",
  fill_color = c("#AB82FF", "#E78AC3", "#A6D854"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
)
print(p1)
p2 <- ggvenn(
  down.list, fill_alpha = 0.5, stroke_color = "white",
  fill_color = c("#AB82FF", "#E78AC3", "#A6D854"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
)
print(p2)
p3 <- ggvenn(
  total.list, fill_alpha = 0.5, stroke_color = "white",
  fill_color = c("#AB82FF", "#E78AC3", "#A6D854"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
)
print(p3)
dev.off()

## comparison with Nkx2.5 target gene/ GATA4 target gene ----
nkx2_5_target <- read.table("/data02/hukaijie/EpiAllele/data/tf_target_gene/Nkx2-5_targets.mouse.tsv")
colnames(nkx2_5_target) <- c("TF", "Target", "Type", "Reference")
gata4_target <- read.table("/data02/hukaijie/EpiAllele/data/tf_target_gene/Gata4_targets.mouse.tsv")
colnames(gata4_target) <- c("TF", "Target", "Type", "Reference")
### SG.vs.NT/ Nkx2.5 target gene/ GATA4 target gene ----
up.list <- list("Nkx2-5 targets" = nkx2_5_target %>% pull(Target),
                "GATA4 targets" = gata4_target %>% pull(Target),
                "SG.vs.NT.up" = SG_NT %>% filter(group =="UP") %>% pull(gene))
down.list <- list("Nkx2-5 targets" = nkx2_5_target %>% pull(Target),
                "GATA4 targets" = gata4_target %>% pull(Target),
                "SG.vs.NT.down" = SG_NT %>% filter(group =="DOWN") %>% pull(gene))
total.list1 <- list("Nkx2-5 targets" = nkx2_5_target %>% pull(Target),
                "GATA4 targets" = gata4_target %>% pull(Target),
                "SG.vs.NT.diff" = SG_NT %>% filter(group != "NOT-SIG") %>% pull(gene))

#### Vennplot ----
pdf("/data02/hukaijie/EpiAllele/result/combine/snRNA.SG.vs.NT.vs.tf_target.vennplot.pdf",
      width = 5, height = 4)
p1 <- ggvenn(
  up.list, fill_alpha = 0.5, stroke_color = "white",
  fill_color = c("#FDB462", "#B3DE69", "#BC80BD"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
)
print(p1)
p2 <- ggvenn(
  down.list, fill_alpha = 0.5, stroke_color = "white",
  fill_color = c("#FDB462", "#B3DE69", "#BC80BD"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
)
print(p2)
p3 <- ggvenn(
  total.list1, fill_alpha = 0.5, stroke_color = "white",
  fill_color = c("#FDB462", "#B3DE69", "#BC80BD"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
)
print(p3)
dev.off()

#### save intersection results ----
intersection_list <- list()
list_collection <- list(
  "UP-regulated" = up.list,
  "DOWN-regulated" = down.list,
  "DIFF-regulated" = total.list1
)
for (collection_name in names(list_collection)) {
   collection <- list_collection[[collection_name]]
   pairwise_intersections <- combn(names(collection), 2, simplify = FALSE) %>%
    map(~{
      a <- .x[1]
      b <- .x[2]
      inter <- intersect(collection[[a]], collection[[b]])
      tibble(comp=paste0(a, ".vs.", b), intersection_size = length(inter), genes = paste(inter, collapse = ","))
    }) %>%
    bind_rows()
  all_sets_intersection <- reduce(collection, intersect)
  intersection_list[[collection_name]] <- rbind(pairwise_intersections, 
                      c("all sets intersections", length(all_sets_intersection), paste(all_sets_intersection, collapse = ",")))
  intersection_list[[collection_name]]$type <- collection_name
}

all_intersection <- do.call(rbind, intersection_list)

write.csv(all_intersection, "/data02/hukaijie/EpiAllele/result/combine/snRNA.SG.vs.NT.vs.tf_target.intersection.csv")

### SG.vs.R404Q/ Nkx2.5 target gene/ GATA4 target gene ----
up.list <- list("Nkx2-5 targets" = nkx2_5_target %>% pull(Target),
                "GATA4 targets" = gata4_target %>% pull(Target),
                "SG.vs.R404Q.up" = SG_R404Q %>% filter(group =="UP") %>% pull(gene))
down.list <- list("Nkx2-5 targets" = nkx2_5_target %>% pull(Target),
                "GATA4 targets" = gata4_target %>% pull(Target),
                "SG.vs.R404Q.down" = SG_R404Q %>% filter(group =="DOWN") %>% pull(gene))
total.list2 <- list("Nkx2-5 targets" = nkx2_5_target %>% pull(Target),
                "GATA4 targets" = gata4_target %>% pull(Target),
                "SG.vs.R404Q.diff" = SG_R404Q %>% filter(group != "NOT-SIG") %>% pull(gene))

#### Vennplot ----
pdf("/data02/hukaijie/EpiAllele/result/combine/snRNA.SG.vs.R404Q.vs.tf_target.vennplot.pdf",
      width = 5, height = 4)
p1 <- ggvenn(
  up.list, fill_alpha = 0.5, stroke_color = "white",
  fill_color = c("#FDB462", "#B3DE69", "#BC80BD"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
)
print(p1)
p2 <- ggvenn(
  down.list, fill_alpha = 0.5, stroke_color = "white",
  fill_color = c("#FDB462", "#B3DE69", "#BC80BD"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
)
print(p2)
p3 <- ggvenn(
  total.list2, fill_alpha = 0.5, stroke_color = "white",
  fill_color = c("#FDB462", "#B3DE69", "#BC80BD"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
)
print(p3)
dev.off()


#### save intersection results ----
intersection_list <- list()
list_collection <- list(
  "UP-regulated" = up.list,
  "DOWN-regulated" = down.list,
  "DIFF-regulated" = total.list2
)
for (collection_name in names(list_collection)) {
   collection <- list_collection[[collection_name]]
   pairwise_intersections <- combn(names(collection), 2, simplify = FALSE) %>%
    map(~{
      a <- .x[1]
      b <- .x[2]
      inter <- intersect(collection[[a]], collection[[b]])
      tibble(comp=paste0(a, ".vs.", b), intersection_size = length(inter), genes = paste(inter, collapse = ","))
    }) %>%
    bind_rows()
  all_sets_intersection <- reduce(collection, intersect)
  intersection_list[[collection_name]] <- rbind(pairwise_intersections, 
                      c("all sets intersections", length(all_sets_intersection), paste(all_sets_intersection, collapse = ",")))
  intersection_list[[collection_name]]$type <- collection_name
}

all_intersection <- do.call(rbind, intersection_list)

write.csv(all_intersection, "/data02/hukaijie/EpiAllele/result/combine/snRNA.SG.vs.R404Q.vs.tf_target.intersection.csv")


### SG.vs.R404Q (specific)/ Nkx2.5 target gene/ GATA4 target gene ----
SG_R404Q_specific_up <- setdiff(SG_R404Q %>% filter(group == "UP") %>% pull(gene),
                              WT_R404Q %>% filter(group == "UP") %>% pull(gene))
SG_R404Q_specific_down <- setdiff(SG_R404Q %>% filter(group == "DOWN") %>% pull(gene),
                              WT_R404Q %>% filter(group == "DOWN") %>% pull(gene))
SG_R404Q_specific_diff <- setdiff(SG_R404Q %>% filter(group != "NOT-SIG") %>% pull(gene),
                              WT_R404Q %>% filter(group != "NOT-SIG") %>% pull(gene))

up.list <- list("Nkx2-5 targets" = nkx2_5_target %>% pull(Target),
                "GATA4 targets" = gata4_target %>% pull(Target),
                "SG.vs.R404Q.up\n(exclude intersection with WT.vs.R404Q.up)" = SG_R404Q_specific_up)
down.list <- list("Nkx2-5 targets" = nkx2_5_target %>% pull(Target),
                "GATA4 targets" = gata4_target %>% pull(Target),
                "SG.vs.R404Q.down\n(exclude intersection with WT.vs.R404Q.down)" = SG_R404Q_specific_down)
total.list <- list("Nkx2-5 targets" = nkx2_5_target %>% pull(Target),
                "GATA4 targets" = gata4_target %>% pull(Target),
                "SG.vs.R404Q.diff\n(exclude intersection with WT.vs.R404Q.diff)" = SG_R404Q_specific_diff)

#### Vennplot ----
pdf("/data02/hukaijie/EpiAllele/result/combine/snRNA.specific.SG.vs.R404Q.vs.tf_target.vennplot.pdf",
      width = 5, height = 4)
p1 <- ggvenn(
  up.list, fill_alpha = 0.5, stroke_color = "white",
  fill_color = c("#FDB462", "#B3DE69", "#BC80BD"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
)
print(p1)
p2 <- ggvenn(
  down.list, fill_alpha = 0.5, stroke_color = "white",
  fill_color = c("#FDB462", "#B3DE69", "#BC80BD"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
)
print(p2)
p3 <- ggvenn(
  total.list, fill_alpha = 0.5, stroke_color = "white",
  fill_color = c("#FDB462", "#B3DE69", "#BC80BD"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
)
print(p3)
dev.off()

#### save intersection results ----
intersection_list <- list()
list_collection <- list(
  "UP-regulated" = up.list,
  "DOWN-regulated" = down.list,
  "DIFF-regulated" = total.list
)
for (collection_name in names(list_collection)) {
   collection <- list_collection[[collection_name]]
   pairwise_intersections <- combn(names(collection), 2, simplify = FALSE) %>%
    map(~{
      a <- .x[1]
      b <- .x[2]
      inter <- intersect(collection[[a]], collection[[b]])
      tibble(comp=paste0(a, ".vs.", b), intersection_size = length(inter), genes = paste(inter, collapse = ","))
    }) %>%
    bind_rows()
  all_sets_intersection <- reduce(collection, intersect)
  intersection_list[[collection_name]] <- rbind(pairwise_intersections, 
                      c("all sets intersections", length(all_sets_intersection), paste(all_sets_intersection, collapse = ",")))
  intersection_list[[collection_name]]$type <- collection_name
}

all_intersection <- do.call(rbind, intersection_list)

write.csv(all_intersection, "/data02/hukaijie/EpiAllele/result/combine/snRNA.specific.SG.vs.R404Q.vs.tf_target.intersection.csv")
