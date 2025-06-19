library(tidyverse)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.refGene)   
library(TxDb.Mmusculus.UCSC.mm39.refGene)
library(DESeq2)
library(EnhancedVolcano)
source("/data02/hukaijie/EpiAllele/final_script/off_target/DE_functions.R")

# offtarget reult ----
mouse_off <- read.table("/data02/hukaijie/EpiAllele/result/off-target/offinder/mouse_offtarget.txt", header=F)
colnames(mouse_off) <- c("Id","Bulge_Type","crRNA","DNA","Chromosome","Location","Direction","Mismatches","Bulge_Size")
human_off <- read.table("/data02/hukaijie/EpiAllele/result/off-target/offinder/human_offtarget.txt", header=F)
colnames(human_off) <- c("Id","Bulge_Type","crRNA","DNA","Chromosome","Location","Direction","Mismatches","Bulge_Size")

## annotation top 10 offtarget site ----
### mouse ----
txdb <- TxDb.Mmusculus.UCSC.mm39.refGene
anndb ="org.Mm.eg.db"
mouse_off_gr <- GRanges(seqnames = mouse_off$Chromosome,
                  ranges = IRanges(start = mouse_off$Location, end = mouse_off$Location),
                  strand = mouse_off$Direction, Mismatches=mouse_off$Mismatches, 
                  Bulge_Type = mouse_off$Bulge_Type, Bulge_Size=mouse_off$Bulge_Size
                  )
mouse_peakAnno <- annotatePeak(mouse_off_gr, TxDb = txdb, annoDb = anndb)
pdf("/data02/hukaijie/EpiAllele/result/off-target/offinder/mouse_offtarget.peakAnnoPie.pdf")
plotAnnoPie(mouse_peakAnno)
dev.off()

mouse_peakAnno_df <- as.data.frame(mouse_peakAnno) %>% 
  filter(SYMBOL != "Myh6") %>% 
  arrange(Mismatches, Bulge_Size)
mouse_top10 <- head(unique(mouse_peakAnno_df$SYMBOL),10)
write.csv(mouse_peakAnno_df, "/data02/hukaijie/EpiAllele/result/off-target/offinder/mouse_offtarget.annotation.csv")

# mouse_peakAnno_target_region_df <- mouse_peakAnno_df %>% 
#   filter(str_detect(annotation, "Exon"))

# length(unique(mouse_peakAnno_target_region_df$SYMBOL))

### human ----
txdb <- TxDb.Hsapiens.UCSC.hg38.refGene
anndb ="org.Hs.eg.db"
human_off_gr <- GRanges(seqnames = human_off$Chromosome,
                  ranges = IRanges(start = human_off$Location, end = human_off$Location),
                  strand = human_off$Direction, Mismatches=human_off$Mismatches, 
                  Bulge_Type = human_off$Bulge_Type, Bulge_Size=human_off$Bulge_Size
                  )
human_peakAnno <- annotatePeak(human_off_gr, TxDb = txdb, annoDb = anndb)
pdf("/data02/hukaijie/EpiAllele/result/off-target/offinder/human_offtarget.peakAnnoPie.pdf")
plotAnnoPie(human_peakAnno)
dev.off()

human_peakAnno_df <- as.data.frame(human_peakAnno) %>% 
  filter(SYMBOL != "MYH7") %>% 
  arrange(Mismatches, Bulge_Size)
human_top10 <- head(unique(human_peakAnno_df$SYMBOL),10)
write.csv(human_peakAnno_df, "/data02/hukaijie/EpiAllele/result/off-target/offinder/human_offtarget.annotation.csv")

# human_peakAnno_target_region_df <- human_peakAnno_df %>%
#   filter(str_detect(annotation, "Exon"))
# length(unique(human_peakAnno_target_region_df$SYMBOL))

# volcano plot ----
# in vivo RNA-seq/snRNA-seq cardiomyocytes/2025 in vitro D2 and D5
out_path <- "/data02/hukaijie/EpiAllele/result/off-target/DE_volcano"
## in vivo RNA-seq ----
DE_path <- "/data02/hukaijie/EpiAllele/result/20240416_mouse_RNAseq/DE"
result_path <- paste0(out_path, "/", "20240416_mouse_RNAseq")
dir.create(result_path, showWarnings=F)

### reverse some DE result ----
groups <- c("NT.vs.SG", "R404Q.vs.SG")
generate_reverse_results(
  DE_path, 
  groups
)
generate_volcano_plots(DE_path, genes_of_interest = c("Myh6", "Myh6_C57", "Myh6_DBA"))

### generate offtarget gene volcano plot ----
generate_top10_volcano_from_DEresults(DE_path, result_path, mouse_peakAnno_df)
#generate_offtarget_volcano_from_DEresults(DE_path, result_path, mouse_peakAnno_target_region_df)

## 2025 invitro D2 & D5 ----
DE_path <- "/data02/hukaijie/EpiAllele/result/20250211_mouse_RNAseq/DE"
result_path <- paste0(out_path, "/", "20250211_mouse_RNAseq")
dir.create(result_path, showWarnings=F)

### reverse some DE result ----
groups <- c("D5NT.vs.D5sg", "D2NT.vs.D2sg")
generate_reverse_results(
  DE_path, 
  groups
)
generate_volcano_plots(DE_path,genes_of_interest = c("Myh6", "Myh6_C57", "Myh6_DBA"))

### generate offtarget gene volcano plot ----
generate_top10_volcano_from_DEresults(DE_path, result_path, mouse_peakAnno_df)
#generate_offtarget_volcano_from_DEresults(DE_path, result_path, mouse_peakAnno_df)
#generate_offtarget_volcano_from_DEresults(DE_path, result_path, mouse_peakAnno_target_region_df)

## 2025 ipsc ----
DE_path <- "/data02/hukaijie/EpiAllele/result/20250212_ipsc_RNAseq/DE"
result_path <- paste0(out_path, "/", "20250212_ipsc_RNAseq")
dir.create(result_path, showWarnings=F)

### reverse some DE result ----
groups <- c("HCMNT.vs.HCMsg")
generate_reverse_results(
  DE_path, 
  groups
)
generate_volcano_plots(DE_path,genes_of_interest = c("MYH7", "MYH7_WT", "MYH7_MUTANT"))

### generate offtarget gene volcano plot ----
generate_top10_volcano_from_DEresults(DE_path, result_path, human_peakAnno_df)
#generate_offtarget_volcano_from_DEresults(DE_path, result_path, human_peakAnno_target_region_df)


## snRNA-seq cardiomyocytes ----
DE_path <- "/data02/hukaijie/EpiAllele/result/snRNAseq/DE"
result_path <- paste0(out_path, "/", "snRNAseq")
dir.create(result_path, showWarnings=F)

#### generate offtarget gene volcano plot ----
generate_top10_volcano_from_DEresults_sn(DE_path, result_path, mouse_peakAnno_df)
#generate_offtarget_volcano_from_DEresults_sn(DE_path, result_path, mouse_peakAnno_target_region_df)
