#!/usr/bin/env Rscript

## {{{ Install required libraries if required
for (p in c( 'stringr', 'ggplot2', 'GenomicRanges', 'tidyverse',"getopt",
             'chromVAR', 'DESeq2', 'ggpubr', 'patchwork', 'EnhancedVolcano',
             'data.table','GenomicAlignments','ggrepel', 'ChIPseeker', 
             'TxDb.Mmusculus.UCSC.mm39.refGene',
             'TxDb.Hsapiens.UCSC.hg38.knownGene')) {
  if (!require(p, character.only = T)) {
    stop(paste0("Please install ", p))
  } else {
    suppressMessages(library(p, quietly = T, character.only = T))
  }
}
## }}}
source("/data02/hukaijie/EpiAllele/final_script/chromatin_profiling/07_DE_function.R")
spec <- matrix(c(
  'result_dir', 'r', 2, "character",
  'species', 's', 2, "character",
  'gene', 'g', 2, "character"
), byrow=TRUE, ncol=4)

opt <- getopt(spec =spec)
out_dir <- opt$result_dir
species <- opt$species
gene <- opt$gene
projPath <- paste0(out_dir, "/macs_single")
bamDir <- paste0(out_dir, "/mapping")
resDir <- paste0(out_dir, "/DE")
dir.create(resDir, showWarnings = F)

peak_files <- list.files(projPath, pattern = "_peaks.xls$", full.names = F)
histList <- unique(sub("(.*)_.*_rep.*", "\\1", peak_files))

for (histone in histList) {
  peak_files <- list.files(projPath, pattern = paste0(histone, ".*_peaks.xls$"), full.names = F)
  sampleList <- unique(sub("(.*)_peaks.xls", "\\1", peak_files))
  ## create a master peak list ----
  mPeak = GRanges()
  for (sample in sampleList) {
    peakRes <- read.table(paste0(projPath, "/", sample, "_peaks.xls"), header = T, check.names = F)
    mPeak = GRanges(seqnames = peakRes$chr, IRanges(start = peakRes$start, end = peakRes$end), strand = "*") %>% append(mPeak, .)
  }
  masterPeak = GenomicRanges::reduce(mPeak)
  saveRDS(masterPeak, paste0(resDir, "/", histone, ".masterPeak.rds"))

  ## peak annotation for masterPeak ----
  if (species == "mouse") {
    txdb <- TxDb.Mmusculus.UCSC.mm39.refGene
    peakAnno <- annotatePeak(masterPeak, tssRegion=c(-3000, 3000),
                          TxDb=txdb, annoDb="org.Mm.eg.db")
  } else if (species == "human") {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    peakAnno <- annotatePeak(masterPeak, tssRegion=c(-3000, 3000),
                          TxDb=txdb, annoDb="org.Hs.eg.db")
  } else {
  stop(paste0(species, " not available"))
  }

  ### peakAnnoPie ----
  pdf(paste0(resDir, "/", histone, ".peakAnnoPie.pdf"))
  plotAnnoPie(peakAnno)
  dev.off()
  ### 
  peakAnno_df <- as.data.frame(peakAnno)
  saveRDS(peakAnno, paste0(resDir, "/", histone, ".peakAnno.rds"))

  ## Get the fragment counts for each peak in the master peak list ----
  countMat = matrix(NA, length(masterPeak), length(sampleList))
  i = 1
  for (sample in sampleList) {
    bamFile = paste0(bamDir, "/", sample, "/", sample, "_bowtie2.mapped.sorted.bam")
    fragment_counts <- getCounts(bamFile, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
    countMat[, i] = counts(fragment_counts)[,1]
    i = i + 1
  }
  colnames(countMat) <- sampleList
  masterPeak_df <- as.data.frame(masterPeak)
  rownames(countMat) <- paste0(masterPeak_df$seqnames,':',masterPeak_df$start,'-',masterPeak_df$end)
  dim(countMat)
  
  
  ## DE analysis ----
  selectR = which(rowSums(countMat) > 5) ## remove low count genes
  dataS = countMat[selectR,]
  dim(dataS)
  sampleTable <- data.frame(condition = factor(sub(".*_(.*)_rep.*", "\\1", sampleList)))
  sampleTable$condition <- relevel(sampleTable$condition, ref = "sg")
  rownames(sampleTable) <- colnames(countMat)
  
  dds = DESeqDataSetFromMatrix(countData = dataS,
                               colData = sampleTable,
                               design = ~ condition)
  DDS = DESeq(dds)
  normDDS = counts(DDS, normalized = TRUE) ## normalization with respect to the sequencing depth
  colnames(normDDS) = paste0(colnames(normDDS), "_norm")
  countMat = cbind(dataS, normDDS)
  head(countMat)
  write.table(countMat, paste0(resDir, "/", histone, ".peak.fragCount.txt"))
  
  ### Volcano plot & whole genome dotplot ----
  comparisons <- combn(unique(levels(sampleTable$condition)), 2, simplify = F)
  comparisons <- lapply(comparisons, function(x) c("condition", x))
  for (comp in comparisons) {
    perform_DE_analysis(DDS, contrast = comp, histone = histone, results.path = resDir, peakAnno_df=peakAnno_df, gene=gene)
  }
}
