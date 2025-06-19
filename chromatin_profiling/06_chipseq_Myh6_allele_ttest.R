library(tidyverse)
library(ggplot2)
library(reshape2)
library(DESeq2)
library(getopt)
library(plyr)


spec <- matrix(c(
  'data_path', 'd', 2, "character",
  'result_path', 'r', 2, "character",
  'bin_size', 'b', 2, "integer"
), byrow=TRUE, ncol=4)

opt <- getopt(spec =spec)
data.path <- opt$data_path
result.path <- opt$result_path
dir.create(result.path, showWarnings=F)
bin_size <- opt$bin_size
if (bin_size %in% c(50, 100, 500)){
  bg_file <- list.files(data.path, pattern = paste0(".bin", bin_size, ".ratio.bedgraph"), full.names = F)
} else{
  stop(paste0("bin_size should be 50, 100 or 500"))
}


bg_file <- list.files(data.path, paste0(".bin", bin_size, ".ratio.bedgraph"), full.names = F)
samples <- unique(sub("(.*?)_.*.bedgraph", "\\1", bg_file))



for (sample in samples) {
  for (group in c("total", "c57", "dba")) {
    bedgraphs <- list.files(data.path, pattern = paste0(sample, ".*", group, ".bin", bin_size, ".ratio.bedgraph"), 
                            full.names = TRUE)
    matrix_list <- list()
    for (bedgraph in bedgraphs) {
      bin_value <- read.table(bedgraph, header = FALSE) %>% 
        mutate(bin = paste0(V1, ":", V2, "-", V3),
              V4 = log2(V4 + 1)) %>%
        select(bin, V4) %>%
        column_to_rownames("bin")
      print(dim(bin_value))
      
      samplename <- sub(paste0("(.*).", group, ".bin", bin_size, ".ratio.bedgraph"), "\\1", basename(bedgraph))
      colnames(bin_value) <- samplename
      matrix_list[[samplename]] <- bin_value
    }
    matrix <- do.call(cbind, matrix_list)
    matrix <- matrix[rowSums(matrix)>0, ]
    

    # sampleinfo
    sampleTable <- data.frame(condition = factor(sub(".*_(.*)_rep.*", "\\1", colnames(matrix))))
    sampleTable$condition <- relevel(sampleTable$condition, ref = "sg")
    rownames(sampleTable) <- colnames(matrix)
    comparisons <- combn(unique(levels(sampleTable$condition)), 2, simplify = F)

    for (comp in comparisons) {
      comp_sample <- sampleTable %>% 
        filter(condition %in% comp)
      comp_matrix <- matrix[, rownames(comp_sample)]
      expr_t <- as.data.frame(t(comp_matrix))
      expr_t$group <- comp_sample$condition
      res <- ddply(melt(expr_t), "variable", function(x) {
        if (length(unique(x$group)) < 2 || any(tapply(x$value, x$group, function(y) var(y) == 0))) {
          return(data.frame(statistic = NA, p.value = NA))
        } else {
          # w <- wilcox.test(value ~ group, data = x, exact = FALSE)
          w <- t.test(value ~ group, data = x)
          means <- tapply(x$value, x$group, mean)
      
          if (length(means) == 2) {
            log2FC <- log2(means[1] / means[2])
          } else {
            log2FC <- NA
          }
          return(with(w, data.frame(statistic, p.value, log2FC)))
        }
      })

      diff_res <- res %>% 
        na.omit() %>% 
        mutate(group = ifelse((p.value<0.05 & log2FC > 1), "UP",
                              ifelse((p.value<0.05 & log2FC < -1), "DOWN", "NOT-SIG"))) %>% 
        mutate(p_adj = p.adjust(p.value, method = "fdr")) %>% 
        mutate(logpval = ifelse(log2FC>0, -log10(p.value), log10(p.value))) %>% 
        mutate(chr=sub("(.*):.*-.*", "\\1", variable),
              start=as.numeric(sub(".*:(.*)-.*", "\\1", variable)),
              end=as.numeric(sub(".*:.*-(.*)", "\\1", variable))) %>% 
        mutate(mid=(start+end)/2)
      
      pdf(paste0(result.path, "/", sample, ".", group, ".", comp[1], ".vs.", comp[2], 
                ".bin", bin_size, ".dotplot.pdf"), width = 6, height = 4)
      p <- ggplot(diff_res, aes(mid, log2FC))+
        geom_point(aes(col=group))+
        theme_bw()+
        scale_color_manual(values=c("UP"="red3", "DOWN"="blue4", "NOT-SIG"="gray"))+
        geom_hline(yintercept = 0)+
        labs(x=NULL,y='log2FC')+
        geom_hline(yintercept = 1, linetype = 'dashed')+
        geom_hline(yintercept = -1, linetype = 'dashed')
      print(p)

      p <- ggplot(diff_res, aes(mid, logpval))+
        geom_point(aes(col=group))+
        theme_bw()+
        scale_color_manual(values=c("UP"="red3", "DOWN"="blue4", "NOT-SIG"="gray"))+
        geom_hline(yintercept = 0)+
        labs(x=NULL,y='-log10(P.value)')+
        geom_hline(yintercept = 1, linetype = 'dashed')+
        geom_hline(yintercept = -1, linetype = 'dashed')
      print(p)

      dev.off()
    }
}
}
