library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)

perform_DE_analysis <- function(dds, contrast, target_geneset, results.path) {
  res <- results(dds, contrast = contrast)
  summary(res)
  diff_res <- as.data.frame(res)
  diff_res <- na.omit(diff_res)
  diff_res$gene <- rownames(diff_res)
  diff_res <- diff_res %>% 
    mutate(group = ifelse((padj<0.05 & log2FoldChange>1), "UP",
                          ifelse((padj<0.05 & log2FoldChange< -1), "DOWN", "NOT-SIG")))
  output_file <- paste0(results.path, contrast[2], ".vs.", contrast[3], ".DEresult.csv")
  write.table(diff_res, file = output_file, sep = ',', row.names = FALSE)
  
  ## volcano plot
  keyvals <- ifelse(
    diff_res$group == "UP", "red3",
    ifelse(diff_res$group == "DOWN", "blue4",
           "gray"))
  names(keyvals) <- diff_res$group
  p1 <- EnhancedVolcano(diff_res,
                        x = 'log2FoldChange',
                        y = 'padj',
                        title = paste0(contrast[2], ".vs.", contrast[3]), 
                        lab = rownames(diff_res),
                        colCustom = keyvals,
                        pCutoff = 0.05, FCcutoff = 1.0,
                        drawConnectors = TRUE,
                        widthConnectors = 0.5,
                        colConnectors = 'black')
  
  
  pdf(paste0(results.path, contrast[2], ".vs.", contrast[3], ".DE.volcano.pdf"),
      width = 6, height = 6)
  print(p1)
  dev.off()
  
  return(diff_res)
}