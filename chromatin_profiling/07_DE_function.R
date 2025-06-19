library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)

perform_DE_analysis <- function(dds, contrast, histone, results.path, peakAnno_df, gene) {
  res <- results(dds, contrast = contrast, independentFiltering = FALSE, altHypothesis = "greaterAbs")
  summary(res)
  diff_res <- as.data.frame(res)
  diff_res <- na.omit(diff_res)
  diff_res$peaks <- rownames(diff_res)
  diff_res <- diff_res %>% 
    mutate(group = ifelse((padj<0.05 & log2FoldChange>1), "UP",
                          ifelse((padj<0.05 & log2FoldChange< -1), "DOWN", "NOT-SIG")))
  output_file <- paste0(results.path, "/", histone, ".", contrast[2], ".vs.", contrast[3], ".peak.DEresult.csv")
  write.table(diff_res, file = output_file, sep = ',', row.names = FALSE)
  

  ## add annotation
  diff_res$chr <- sub("(chr.*):.*-.*", "\\1", diff_res$peaks)
  peakAnno_df$peaks <- paste0(peakAnno_df$seqnames, ":", peakAnno_df$start, "-", peakAnno_df$end)
  diff_gene_res <- diff_res %>% 
    inner_join(peakAnno_df, by="peaks") %>% 
    na.omit() %>% 
    dplyr::select(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, peaks, chr, start, end, width, strand, annotation,
          geneChr, geneStart, geneEnd, geneLength, geneStrand, geneId, transcriptId, distanceToTSS,
          ENSEMBL, SYMBOL, group) %>% 
    mutate(start=as.numeric(start),
          end=as.numeric(end),
          mid=(start+end)/2,
          chr=factor(chr, levels=c(paste0("chr", 1:19), "chrX", "chrY", "chrM")))
  
  ## volcano plot ----
  rownames(diff_gene_res) <- diff_gene_res$peaks
  top_region <- diff_gene_res %>%
    filter(group != "NOT-SIG") %>% 
    arrange(pvalue, desc(abs(log2FoldChange))) %>%
    slice_head(n = 4) %>% 
    pull(peaks)
  
  diff_gene_res <- diff_gene_res %>% 
    mutate(group_new=ifelse(SYMBOL==gene, paste0(gene, "-related"), group)) %>% 
    arrange(group_new == paste0(gene, "-related"))
  # keyvals <- ifelse(
  #   diff_gene_res$group_new == "UP", "red3",
  #   ifelse(diff_gene_res$group_new == "DOWN", "blue4",
  #          ifelse(diff_gene_res$group_new == paste0(gene, "-related"), "green4", "gray")))
  # names(keyvals) <- diff_gene_res$group_new
  
  # target_gene_related_peaks <- diff_gene_res %>% 
  #   filter(group_new==paste0(gene, "-related")) %>% 
  #   pull(peaks)

  keyvals <- ifelse(
    diff_gene_res$group == "UP", "red3",
    ifelse(diff_gene_res$group == "DOWN", "blue4", "gray"))
  names(keyvals) <- diff_gene_res$group
  
  target_gene_related_peaks <- diff_gene_res %>% 
    filter(group_new==paste0(gene, "-related")) %>% 
    pull(peaks)
  if (length(target_gene_related_peaks) < 5) {
    top_target_region <- c(top_region, target_gene_related_peaks)
  } else{
    top_target_region <- c(top_region, head(target_gene_related_peaks,4))
  }
  ## show target gene related peaks
  p1 <- EnhancedVolcano(diff_gene_res,
                        x = 'log2FoldChange',
                        y = 'padj',
                        title = paste0(contrast[2], ".vs.", contrast[3]), 
                        lab = rownames(diff_gene_res),
                        selectLab = top_region,
                        colCustom = keyvals,
                        pCutoff = 0.05, FCcutoff = 1,
                        drawConnectors = TRUE,
                        widthConnectors = 0.5,
                        colConnectors = 'black',
                        labSize = 3)
  ## show target gene related peaks + Myh6 related peaks
  keyvals <- ifelse(
    diff_gene_res$group_new == "UP", "red3",
    ifelse(diff_gene_res$group_new == "DOWN", "blue4",
           ifelse(diff_gene_res$group_new == paste0(gene, "-related"), "green4", "gray")))
  names(keyvals) <- diff_gene_res$group_new
  p2 <- EnhancedVolcano(diff_gene_res,
                        x = 'log2FoldChange',
                        y = 'padj',
                        title = paste0(contrast[2], ".vs.", contrast[3]), 
                        lab = rownames(diff_gene_res),
                        selectLab = top_target_region,
                        colCustom = keyvals,
                        pCutoff = 0.05, FCcutoff = 1,
                        drawConnectors = TRUE,
                        widthConnectors = 0.5,
                        colConnectors = 'black',
                        labSize = 3)

  pdf(paste0(results.path, "/", histone, ".",  contrast[2], ".vs.", contrast[3], ".onlyTopRegion.DE.volcano.pdf"),
      width = 6, height = 6)
  print(p1)
  dev.off()

  pdf(paste0(results.path, "/", histone, ".",  contrast[2], ".vs.", contrast[3], ".TopRegionPlus", gene, ".DE.volcano.pdf"),
      width = 6, height = 6)
  print(p2)
  dev.off()

  png(paste0(results.path, "/", histone, ".",  contrast[2], ".vs.", contrast[3], ".onlyTopRegion.DE.volcano.png"),
      width = 6, height = 6, res = 300, units = 'in')
  print(p1)
  dev.off()

  png(paste0(results.path, "/", histone, ".",  contrast[2], ".vs.", contrast[3], ".TopRegionPlus", gene, ".DE.volcano.png"),
      width = 6, height = 6, res = 300, units = 'in')
  print(p2)
  dev.off()

  ## whole genome dot plot
  p1 <- ggplot(diff_gene_res, aes(mid, log2FoldChange))+
    geom_point(aes(col = group_new), size = 0.5) +
    theme_bw()+
    scale_color_manual(values=c("UP"="red3", "DOWN"="blue4", "NOT-SIG"="gray", 
                                setNames("green4", paste0(gene, "-related"))))+
    geom_hline(yintercept = 1, linetype = 'dashed')+
    geom_hline(yintercept = -1, linetype = 'dashed')+
    facet_wrap(~chr, nrow=1)+
    theme(axis.title.x = element_blank())+
    scale_x_continuous(breaks = NULL)

  png(paste0(results.path, "/", histone, ".",  contrast[2], ".vs.", contrast[3], ".DE.logFC.png"), 
  width = 12, height = 6, res = 300, units = 'in')
  print(p1)
  dev.off()

  pdf(paste0(results.path, "/", histone, ".",  contrast[2], ".vs.", contrast[3], ".DE.logFC.pdf"), 
    width = 12, height = 6)
  print(p1)
  dev.off()

  write.table(diff_gene_res, file = output_file, sep = ',', row.names = FALSE)
}
