library(tidyverse)
library(EnhancedVolcano)

source("/data02/hukaijie/EpiAllele/final_script/wgbs/07_DMR_DE_functions.R")
# day2 DMRs expression ----
result.path <- "/data02/hukaijie/EpiAllele/result/20250310_wgbs/results/"
D2_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/20250211_mouse_RNAseq/DE/D2sg.vs.D2NT.DEresult.csv",
                    header=T)
D2_DMR <- read.csv("/data02/hukaijie/EpiAllele/result/20250310_wgbs/results/anno_dmrs.csv", header=T, row.names=1)
D2_DEG_in_DMR <- D2_DEG %>%
  filter(gene %in% D2_DMR$SYMBOL)
rownames(D2_DEG_in_DMR) <- D2_DEG_in_DMR$gene

generate_volcano_plots(D2_DEG_in_DMR, result.path, filename="D2sg.vs.D2NT", gene_of_interest = "Myh6")


# day5 DMRs expression ----
result.path <- "/data02/hukaijie/EpiAllele/result/20250212_wgbs/results/"
D5_DEG <- read.csv("/data02/hukaijie/EpiAllele/result/20250211_mouse_RNAseq/DE/D5sg.vs.D5NT.DEresult.csv",
                    header=T)
D5_DMR <- read.csv("/data02/hukaijie/EpiAllele/result/20250212_wgbs/results/anno_dmrs.csv", header=T, row.names=1)
D5_DEG_in_DMR <- D5_DEG %>%
  filter(gene %in% D5_DMR$SYMBOL)
rownames(D5_DEG_in_DMR) <- D5_DEG_in_DMR$gene

generate_volcano_plots(D5_DEG_in_DMR, result.path, filename="D5sg.vs.D5NT", gene_of_interest = "Myh6")

# invivo DMRs expression ----
result.path <- "/data02/hukaijie/EpiAllele/result/20240521_wgbs/results/"
group_list <- list(c("sg.vs.NT", "SQ.vs.NT"))
for (group in group_list){
  DMR_group <- group[1]
  DE_group <- group[2]
  invivo_DEG <- read.csv(paste0("/data02/hukaijie/EpiAllele/result/snRNAseq/DE/Cardiomyocytes.", DE_group, ".DEresult.csv"),
                      header=T, row.names=1)
  invivo_DMR <- read.csv(paste0("/data02/hukaijie/EpiAllele/result/20240521_wgbs/results/", DMR_group, ".anno_dmrs.csv"), 
                        header=T, row.names=1)
  invivo_DEG_in_DMR <- invivo_DEG %>%
    filter(gene %in% invivo_DMR$SYMBOL)
  rownames(invivo_DEG_in_DMR) <- invivo_DEG_in_DMR$gene

  generate_volcano_plots_sn(invivo_DEG_in_DMR, result.path, filename=DMR_group, gene_of_interest = "Myh6")
}
