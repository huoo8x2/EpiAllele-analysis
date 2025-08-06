library(dplyr)
library(EnhancedVolcano)

generate_reverse_results <- function(DE_path, groups) {
  for (filename in groups) {
    
    diff_res <- read.csv(file.path(DE_path, paste0(filename, ".DEresult.csv")))
    rownames(diff_res) <- diff_res$gene

    diff_res_rev <- diff_res %>% 
      mutate(log2FoldChange = -log2FoldChange,
             group = ifelse(group == "UP", "NOT-UP",
                            ifelse(group == "DOWN", "UP", group))) %>%
      mutate(group = ifelse(group == "NOT-UP", "DOWN", group))

    comp1 <- sub("(.*).vs.(.*)", "\\1", filename)
    comp2 <- sub("(.*).vs.(.*)", "\\2", filename)
    new_filename <- paste0(comp2, ".vs.", comp1)

    write.table(diff_res_rev, 
                file = file.path(DE_path, paste0(new_filename, ".DEresult.csv")),
                sep = ',', row.names = FALSE)
    
    # delete previous DE results
    file.remove(file.path(DE_path, paste0(filename, ".DEresult.csv")))
  }
}

generate_volcano_plots <- function(DE_path, genes_of_interest = c("Myh6", "Myh6_C57", "Myh6_DBA")) {
  DEresult_files <- list.files(DE_path, pattern = ".DEresult.csv$", full.names = TRUE)
  for (file in DEresult_files) {
    DE_results <- read.csv(file, header = TRUE)
    rownames(DE_results) <- DE_results$gene
    filename <- sub(".*/(.*).DEresult.csv", "\\1", file)
    keyvals <- ifelse(
      DE_results$group == "UP", "red3",
      ifelse(DE_results$group == "DOWN", "blue4", "gray")
    )
    names(keyvals) <- DE_results$group
    top5_genes <- DE_results %>% 
      filter(group %in% c("DOWN", "UP")) %>% 
      arrange(padj, -abs(log2FoldChange)) %>% 
      do(head(., n = 5)) %>% 
      pull(gene)
    # volcano plot 1：all genes
    p1 <- EnhancedVolcano(DE_results,
                          x = 'log2FoldChange',
                          y = 'padj',
                          title = filename, 
                          selectLab = top5_genes,
                          lab = rownames(DE_results),
                          colCustom = keyvals,
                          pCutoff = 0.05, FCcutoff = 1.0,
                          drawConnectors = TRUE,
                          widthConnectors = 0.5,
                          colConnectors = 'black')
    # volcano plot 2：specific genes
    if (length(genes_of_interest) > 0){
      DE_results <- DE_results %>% 
        mutate(group_new=ifelse(gene %in% genes_of_interest, "Myh6-related", group)) %>% 
        arrange(gene %in% genes_of_interest)
      keyvals <- ifelse(
        DE_results$group_new == "UP", "red3",
        ifelse(DE_results$group_new == "DOWN", "blue4", 
        ifelse(DE_results$group_new == "Myh6-related", "green4","gray"))
      )
      names(keyvals) <- DE_results$group_new
      p2 <- EnhancedVolcano(DE_results,
                          x = 'log2FoldChange',
                          y = 'padj',
                          selectLab = c(top5_genes, genes_of_interest),
                          title = filename, 
                          lab = rownames(DE_results),
                          colCustom = keyvals,
                          pCutoff = 0.05, FCcutoff = 1.0,
                          drawConnectors = TRUE,
                          widthConnectors = 0.5,
                          colConnectors = 'black')
      pdf(str_replace(file, ".DEresult.csv", ".DE.volcano.pdf"),
        width = 6, height = 6)
      print(p1)
      print(p2)
      dev.off()
    } else{
      pdf(str_replace(file, ".DEresult.csv", ".DE.volcano.pdf"),
          width = 6, height = 6)
      print(p1)
      dev.off()
    }
  }
}


generate_top10_volcano_from_DEresults <- function(DE_path, result_path, peakAnno_df) {
  DEresult_files <- list.files(DE_path, pattern = ".DEresult.csv$", full.names = TRUE)
  
  for (file in DEresult_files) {
    DE_results <- read.csv(file, header = TRUE)
    rownames(DE_results) <- DE_results$gene

    mouse_top10 <- head(intersect(unique(peakAnno_df$SYMBOL), rownames(DE_results)), 10)

    top10_DE_results <- DE_results %>% 
      filter(gene %in% mouse_top10)

    filename <- sub(".*/(.*).DEresult.csv", "\\1", file)

    keyvals <- ifelse(
      top10_DE_results$group == "UP", "red3",
      ifelse(top10_DE_results$group == "DOWN", "blue4", "gray")
    )
    names(keyvals) <- top10_DE_results$group

    p1 <- EnhancedVolcano(top10_DE_results,
                          x = 'log2FoldChange',
                          y = 'padj',
                          selectLab = top10_DE_results$gene,
                          title = filename, 
                          max.overlaps=3,
                          lab = rownames(top10_DE_results),
                          colCustom = keyvals,
                          pCutoff = 0.05, FCcutoff = 1.0,
                          drawConnectors = TRUE,
                          widthConnectors = 0.5,
                          colConnectors = 'black')

    pdf(file.path(result_path, paste0(filename, ".top10.offtargetGene.DE.volcano.pdf")),
        width = 6, height = 6)
    print(p1)
    dev.off()

    write.csv(top10_DE_results, 
              file.path(result_path, paste0(filename, ".top10.offtargetGene.DEresult.csv")),
              row.names = FALSE)
  }
}

generate_top10_volcano_from_DEresults_sn <- function(DE_path, result_path, peakAnno_df) {
  DEresult_files <- list.files(DE_path, pattern = ".DEresult.csv$", full.names = TRUE)
  
  for (file in DEresult_files) {
    DE_results <- read.csv(file, header = TRUE)
    rownames(DE_results) <- DE_results$gene

    mouse_top10 <- head(intersect(unique(peakAnno_df$SYMBOL), rownames(DE_results)), 10)

    top10_DE_results <- DE_results %>% 
      filter(gene %in% mouse_top10)

    filename <- sub(".*/(.*).DEresult.csv", "\\1", file)

    keyvals <- ifelse(
      top10_DE_results$group == "UP", "red3",
      ifelse(top10_DE_results$group == "DOWN", "blue4", "gray")
    )
    names(keyvals) <- top10_DE_results$group

    p1 <- EnhancedVolcano(top10_DE_results,
                          x = "avg_log2FC", 
                          y = "p_val_adj",
                          selectLab = top10_DE_results$gene,
                          title = filename, 
                          max.overlaps=3,
                          lab = rownames(top10_DE_results),
                          colCustom = keyvals,
                          pCutoff = 0.05, FCcutoff = 0.5,
                          drawConnectors = TRUE,
                          widthConnectors = 0.5,
                          colConnectors = 'black')

    # 输出 PDF
    pdf(file.path(result_path, paste0(filename, ".top10.offtargetGene.DE.volcano.pdf")),
        width = 6, height = 6)
    print(p1)
    dev.off()

    write.csv(top10_DE_results, 
              file.path(result_path, paste0(filename, ".top10.offtargetGene.DEresult.csv")),
              row.names = FALSE)
  }
}

generate_offtarget_volcano_from_DEresults <- function(DE_path, result_path, peakAnno_df) {
  DEresult_files <- list.files(DE_path, pattern = ".DEresult.csv$", full.names = TRUE)
  
  for (file in DEresult_files) {
    # DE results
    DE_results <- read.csv(file, header = TRUE)
    rownames(DE_results) <- DE_results$gene

    # extract intersected offtarget genes
    offtarget_genes <- intersect(unique(peakAnno_df$SYMBOL), rownames(DE_results))
    offtarget_DE_results <- DE_results %>% 
      filter(gene %in% offtarget_genes)
    if (nrow(offtarget_DE_results) > 0) {
      filename <- sub(".*/(.*).DEresult.csv", "\\1", file)
      keyvals <- ifelse(
        offtarget_DE_results$group == "UP", "red3",
        ifelse(offtarget_DE_results$group == "DOWN", "blue4", "gray")
      )
      names(keyvals) <- offtarget_DE_results$group

      # volcano plot
      p1 <- EnhancedVolcano(offtarget_DE_results,
                            x = 'log2FoldChange',
                            y = 'padj',
                            selectLab = offtarget_DE_results$gene,
                            title = filename, 
                            max.overlaps=3,
                            lab = rownames(offtarget_DE_results),
                            colCustom = keyvals,
                            pCutoff = 0.05, FCcutoff = 1.0,
                            drawConnectors = TRUE,
                            widthConnectors = 0.5,
                            colConnectors = 'black')
      pdf(file.path(result_path, paste0(filename, ".offtargetGene.DE.volcano.pdf")),
          width = 6, height = 6)
      print(p1)
      dev.off()

      write.csv(offtarget_DE_results, 
                file.path(result_path, paste0(filename, ".offtargetGene.DEresult.csv")),
                row.names = FALSE)
    }
    
  }
}

generate_offtarget_volcano_from_DEresults_sn <- function(DE_path, result_path, peakAnno_df) {
  DEresult_files <- list.files(DE_path, pattern = ".DEresult.csv$", full.names = TRUE)
  
  for (file in DEresult_files) {
    # DE results
    DE_results <- read.csv(file, header = TRUE)
    rownames(DE_results) <- DE_results$gene

    # extract intersected offtarget genes
    offtarget_genes <- intersect(unique(peakAnno_df$SYMBOL), rownames(DE_results))
    offtarget_DE_results <- DE_results %>% 
      filter(gene %in% offtarget_genes)
    
    if (nrow(offtarget_DE_results) > 0) {
      filename <- sub(".*/(.*).DEresult.csv", "\\1", file)
      keyvals <- ifelse(
        offtarget_DE_results$group == "UP", "red3",
        ifelse(offtarget_DE_results$group == "DOWN", "blue4", "gray")
      )
      names(keyvals) <- offtarget_DE_results$group

      p1 <- EnhancedVolcano(offtarget_DE_results,
                            x = "avg_log2FC", 
                            y = "p_val_adj",
                            selectLab = offtarget_DE_results$gene,
                            title = filename, 
                            max.overlaps=3,
                            lab = rownames(offtarget_DE_results),
                            colCustom = keyvals,
                            pCutoff = 0.05, FCcutoff = 0.5,
                            drawConnectors = TRUE,
                            widthConnectors = 0.5,
                            colConnectors = 'black')

      pdf(file.path(result_path, paste0(filename, ".offtargetGene.DE.volcano.pdf")),
          width = 6, height = 6)
      print(p1)
      dev.off()

      write.csv(offtarget_DE_results, 
                file.path(result_path, paste0(filename, ".offtargetGene.DEresult.csv")),
                row.names = FALSE)
    }
  }
}
