generate_volcano_plots <- function(DE_results, result.path, filename, gene_of_interest = "Myh6") {
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
  # volcano plot 2：Only marked gene of interest
  if (sum(DE_results$gene == gene_of_interest) > 0) {
     DE_results <- DE_results %>% 
      mutate(group_new = ifelse(gene == gene_of_interest, gene_of_interest, group)) %>% 
      arrange(gene == gene_of_interest)
    keyvals <- ifelse(
      DE_results$group_new == "UP", "red3",
      ifelse(DE_results$group_new == "DOWN", "blue4", 
      ifelse(DE_results$group_new == gene_of_interest, "green4","gray"))
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
    pdf(paste0(result.path, filename, ".DMRs.expr.DE.volcano.pdf"),
        width = 6, height = 6)
    print(p1)
    print(p2)
    dev.off()
  } else{
    pdf(paste0(result.path, filename, ".DMRs.expr.DE.volcano.pdf"),
        width = 6, height = 6)
    print(p1)
    dev.off()
  }
}

generate_volcano_plots_sn <- function(DE_results, result.path, filename, gene_of_interest = "Myh6") {
  keyvals <- ifelse(
    DE_results$group == "UP", "red3",
    ifelse(DE_results$group == "DOWN", "blue4", "gray")
  )
  names(keyvals) <- DE_results$group
  top5_genes <- DE_results %>% 
    filter(group %in% c("DOWN", "UP")) %>% 
    arrange(p_val_adj, -abs(avg_log2FC)) %>% 
    do(head(., n = 5)) %>% 
    pull(gene)
  # volcano plot 1：all genes
  p1 <- EnhancedVolcano(DE_results,
                        x = 'avg_log2FC',
                        y = 'p_val_adj',
                        title = filename, 
                        selectLab = top5_genes,
                        lab = rownames(DE_results),
                        colCustom = keyvals,
                        pCutoff = 0.05, FCcutoff = 0.5,
                        drawConnectors = TRUE,
                        widthConnectors = 0.5,
                        colConnectors = 'black')
  # volcano plot 2：Only marked gene of interest
  if (sum(DE_results$gene == gene_of_interest) > 0) {
     DE_results <- DE_results %>% 
      mutate(group_new = ifelse(gene == gene_of_interest, gene_of_interest, group)) %>% 
      arrange(gene == gene_of_interest)
    keyvals <- ifelse(
      DE_results$group_new == "UP", "red3",
      ifelse(DE_results$group_new == "DOWN", "blue4", 
      ifelse(DE_results$group_new == gene_of_interest, "green4","gray"))
    )
    names(keyvals) <- DE_results$group_new
    p2 <- EnhancedVolcano(DE_results,
                        x = 'avg_log2FC',
                        y = 'p_val_adj',
                        selectLab = c(top5_genes, genes_of_interest),
                        title = filename, 
                        lab = rownames(DE_results),
                        colCustom = keyvals,
                        pCutoff = 0.05, FCcutoff = 0.5,
                        drawConnectors = TRUE,
                        widthConnectors = 0.5,
                        colConnectors = 'black')
    pdf(paste0(result.path, filename, ".DMRs.expr.DE.volcano.pdf"),
        width = 6, height = 6)
    print(p1)
    print(p2)
    dev.off()
  } else{
    pdf(paste0(result.path, filename, ".DMRs.expr.DE.volcano.pdf"),
        width = 6, height = 6)
    print(p1)
    dev.off()
  }
}