#!/usr/bin/env Rscript

## {{{ Install required libraries if required
for (p in c("tidyverse", "stringr", "ggrepel", "EnhancedVolcano", "DESeq2", "ggsci", "getopt",
            "biomaRt", "clusterProfiler", "org.Hs.eg.db", "ggpubr", "ComplexHeatmap")) {
  if (!require(p, character.only = T)) {
    print(paste0("Please install ", p))
  } else {
    suppressMessages(library(p, quietly = T, character.only = T))
  }
}
## }}}
getOption('timeout')
options(timeout=10000)

spec <- matrix(c(
  'result_path', 'r', 2, "character",
  'species', 's', 2, "character",
  'snp_path', 'p', 2, "character"
), byrow=TRUE, ncol=4)


opt <- getopt(spec =spec)
results.path <- opt$result_path
species <- opt$species
snp_path <- opt$snp_path

dir.create(results.path, showWarnings = F)
featureCount.path <- paste0(results.path, "/hisat2/featureCount")

source("/data02/hukaijie/EpiAllele/final_script/RcolorSet.R")
# Quality check ----
summary.files <- list.files(featureCount.path, pattern = ".summary$", full.names = T)
mapping.path <- paste0(results.path, "/mapping_quality/")
dir.create(mapping.path, showWarnings = F)

summary_list <- list()
for (file in summary.files) {
  sample <- sub(".*/(.*?-.*?)-.*.summary$", "\\1", file)
  summary <- read.table(file, header = T, check.names = F)
  colnames(summary) <- c("status", sample)
  summary_list[[sample]] <- summary
}
summary_info <- Reduce(function(x, y) merge(x, y, by = "status"), summary_list)
head(summary_info)

rownames(summary_info) <- summary_info$status
summary_info$status <- NULL
summary_info <- t(summary_info)
summary_info <- as.data.frame(summary_info)
summary_info$mapping_rate <- summary_info$Assigned/rowSums(summary_info)
summary_info$sample <- rownames(summary_info)
head(summary_info)

write.csv(summary_info, paste0(mapping.path, "mapping_summary.csv"), row.names = F)

## plot the number of reads per sample----
pdf(paste0(mapping.path, "qc_plot_reads_per_sample.pdf"), width=8, height=6)
p <- ggplot(summary_info, aes(x=sample, y=Assigned)) + 
  geom_col() + 
  theme_bw() + 
  ggtitle("Number of reads per sample") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(p)
dev.off()

## plot the mapping rate per sample ----
pdf(paste0(mapping.path, "qc_plot_mapping_rate_per_sample.pdf"), width=8, height=6)
p <- ggplot(summary_info, aes(x=sample, y=mapping_rate)) + 
  geom_col() + 
  theme_bw() + 
  ggtitle("Mapping rate per sample") + 
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(p)
dev.off()

# Generate expr matrix ----
expr.path <- paste0(results.path, "/expr/")
dir.create(expr.path, showWarnings = F)
files <- list.files(featureCount.path, pattern = ".txt$", full.names = T)

counts_list <- list()
for (file in files) {
  sample <- sub(".*/(.*?-.*?)-.*.txt$", "\\1", file)
  counts <- read.table(file, header = T, check.names = F)
  counts <- counts[, c(1, ncol(counts))]
  colnames(counts) <- c("geneid", sample)
  counts_list[[sample]] <- counts
}
expr <- Reduce(function(x, y) merge(x, y, by = "geneid"), counts_list)
head(expr)

# geneid -> gene symbol ----
## geneid <-> gene symbol index ----
if (species == "human") {
  gtf <- rtracklayer::import("/data02/hukaijie/EpiAllele/ref/human/GRCh38/Homo_sapiens.GRCh38.113.chr.gtf") %>% 
    as.data.frame()
} else if (species == "mouse") {
  gtf <- rtracklayer::import("/data02/hukaijie/EpiAllele/ref/mouse/GRCm39/Mus_musculus.GRCm39.113.chr.gtf") %>% 
    as.data.frame()
} else {
  stop(paste0(species, " not available"))
}
id_symbol_idx <- gtf %>% 
  filter(gene_biotype=="protein_coding", type=="gene") %>% 
  dplyr::select(gene_id, gene_name)
head(id_symbol_idx)
id_symbol_idx <- na.omit(id_symbol_idx)
id_symbol_idx <- unique(id_symbol_idx)

expr <- merge(expr, id_symbol_idx, by.x = "geneid", by.y="gene_id")
head(expr)
dim(expr)
length(unique(expr$gene_name))
length(unique(expr$geneid))
table(duplicated(expr$gene_name))

## remove duplicated genes (preserve max value) ----
expr <- aggregate(.~gene_name, max, data=expr)

### prepare genelength 
gene_length <- gtf %>% 
  filter(gene_biotype=="protein_coding", type=="gene", gene_id %in% expr$geneid) %>% 
  dplyr::select(gene_id, gene_name, width)

rownames(expr) <- expr$gene_name
expr$geneid <- NULL
expr$gene_name <- NULL
head(expr)
expr <- expr %>% mutate(across(where(is.character), as.numeric))
expr <- expr[rowSums(expr) > 0,]
write.table(expr, paste0(expr.path, "counts.txt"))


# calculate TPM ----
rownames(gene_length) <- gene_length$gene_name
gene_length <- gene_length[rownames(expr),]
countsTPM <- function(count, efflength){
  RPK <- count/(efflength/1000)
  PMSC_rpk <- sum(RPK)/1e6
  RPK/PMSC_rpk
}
tpm <- data.frame(apply(expr, 2, countsTPM, efflength=gene_length$width), check.names = F)
write.table(tpm, paste0(expr.path, "tpm.txt"))

# expression correlation ----
cor.result <- cor(scale(tpm), method = "spearman")

## sampel correaltion ----
my_colors <- colorRampPalette(c(
"#619DB8", "#AECDD7", "#E3EEEF",
"#FAE7D9", "#F0B79A", "#C85D4D"
))(100)
pdf(paste0(expr.path, "sample.cor.result.pdf"), width = 10, height = 8)
p <- pheatmap(cor.result, cluster_rows = T, cluster_cols = T, 
              display_numbers = T, number_format = "%.3f", color = my_colors)
print(p)
dev.off()

## mean expression correaltion ----
sampleTable <- data.frame(condition = factor(sub("(.*)-.*", "\\1", colnames(expr))))
rownames(sampleTable) <- colnames(expr)
sampleTable


###create sampleTable: sample-condition mapping
sample_list <- list()
sampleTable$sample <- rownames(sampleTable)
for (group in unique(sampleTable$condition)) {
  samples <- sampleTable %>%
    filter(condition == group) %>%
    pull(sample)
  sample_list[[group]] <- samples
}

tpm_mean <- sapply(sample_list, function(cols) rowMeans(tpm[, cols]))
cor.result <- cor(scale(tpm_mean), method = "spearman")
p <- pheatmap(cor.result, cluster_rows = F, cluster_cols = F,
              display_numbers = T, number_format = "%.3f")
pdf(paste0(expr.path, "meanexpr.cor.result.pdf"), width = 6, height = 4)
print(p)
dev.off()

# Run PCA ----
tpm <- tpm[which(rowSums(tpm) != 0),]
tpm4pca <- log10(tpm+1)
tpm4pca <- t(tpm4pca)
tpm4pca <- as.data.frame(tpm4pca)

tpm4pca$sample <- sampleTable$condition
fit <- prcomp(tpm4pca[,-ncol(tpm4pca)], center = T, scale. = T)
pca_pcs <- cbind(as.data.frame(fit$x),
                 Group=tpm4pca$sample,
                 sample=rownames(tpm4pca))
pca_var <- data.frame(pc=colnames(pca_pcs)[1:(ncol(pca_pcs)-2)], 
                      variance=round((fit$sdev^2)/sum(fit$sdev^2) *100, 2))
pca_var$pc <- factor(pca_var$pc, levels = pca_var$pc)

pdf(paste0(expr.path, "pcs_contributions_plot.pdf"), widt=10, height=6)
p <- ggplot(pca_var, aes(pc, variance, fill = pc)) +
  geom_bar(stat = 'identity')+
  scale_fill_d3(palette = "category20")+
  theme_bw() + labs(x = 'PCs', y = "Pcs Contributions(%)")
print(p)
dev.off()

pdf(paste0(expr.path, "pca.pdf"),width = 6,height = 4)
p <- ggplot(pca_pcs, aes(x=PC1, y=PC2, col=Group)) + 
  geom_point(size=2.5) + 
  theme_classic() + 
  scale_color_npg()+
  geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE, 
                  box.padding = unit(0.5, 'lines'))+
  #stat_ellipse(type = "t", linetype = 1)+
  ggtitle("PCA plot") + 
  xlab(sprintf("PC1 (%0.2f%%)", pca_var[1, 2])) + 
  ylab(sprintf("PC2 (%0.2f%%)", pca_var[2, 2]))
print(p)
p <- ggplot(pca_pcs, aes(x=PC1, y=PC3, col=Group)) + 
  geom_point(size=2.5) + 
  theme_classic() + 
  scale_color_npg()+
  geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE, 
                  box.padding = unit(0.5, 'lines'))+
  #stat_ellipse(type = "t", linetype = 1)+
  ggtitle("PCA plot") + 
  xlab(sprintf("PC1 (%0.2f%%)", pca_var[1, 2])) + 
  ylab(sprintf("PC3 (%0.2f%%)", pca_var[3, 2]))
print(p)
p <- ggplot(pca_pcs, aes(x=PC2, y=PC3, col=Group)) + 
  geom_point(size=2.5) + 
  theme_classic() + 
  scale_color_npg()+
  geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE, 
                  box.padding = unit(0.5, 'lines'))+
  #stat_ellipse(type = "t", linetype = 1)+
  ggtitle("PCA plot") + 
  xlab(sprintf("PC2 (%0.2f%%)", pca_var[2, 2])) + 
  ylab(sprintf("PC3 (%0.2f%%)", pca_var[3, 2]))
print(p)
p <- ggplot(pca_pcs, aes(x=PC2, y=PC4, col=Group)) + 
  geom_point(size=2.5) + 
  theme_classic() + 
  scale_color_npg()+
  geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE, 
                  box.padding = unit(0.5, 'lines'))+
  #stat_ellipse(type = "t", linetype = 1)+
  ggtitle("PCA plot") + 
  xlab(sprintf("PC2 (%0.2f%%)", pca_var[2, 2])) + 
  ylab(sprintf("PC4 (%0.2f%%)", pca_var[4, 2]))
print(p)
dev.off()



# DE analysis ----
de.path <- paste0(results.path, "/DE/")
dir.create(de.path, showWarnings = F)
allele_expr <- read.csv(snp_path, header = T, row.names = 1)
rownames(allele_expr) <- sub("(.*?-.*?)-.*", "\\1", rownames(allele_expr))


if (species == "mouse") {
  allele1="C57"
  allele2="DBA"
  target_gene="Myh6"
  target_gene_list <- c("Myh6", "Myh6_C57", "Myh6_DBA")
} else if (species == "human") {
  allele1="MUTANT"
  allele2="WT"
  target_gene="MYH7"
  target_gene_list <- c("MYH7", "MYH7_WT", "MYH7_MUTANT")
} else {
  stop("species should be mouse or human")
}

allele1_name=paste0(target_gene, "_", allele1)
allele2_name=paste0(target_gene, "_", allele2)

allele_expr <- allele_expr %>% 
  dplyr::select(all_of(allele1), all_of(allele2)) %>%
  rename(!!allele1_name := all_of(allele1),
         !!allele2_name := all_of(allele2)) %>%
  t() %>% 
  as.data.frame()

allele_expr <- rbind(expr, allele_expr)

dds <- DESeqDataSetFromMatrix(countData=allele_expr, 
                              colData=sampleTable, design= ~condition)
dds <- DESeq(dds)
saveRDS(dds, paste0(de.path, "dds.rds"))

# Construct a list containing all the pairs
comparisons <- combn(unique(as.character(sampleTable$condition)), 2, simplify = F)
comparisons <- lapply(comparisons, function(x) c("condition", x))

# Get the DEG results for all the pairs
diff_results_list <- list()
target_diff_list <- list()
source("/data02/hukaijie/EpiAllele/final_script/bulk_RNAseq/04_DE_function.R")
for (comp in comparisons) {
  diff_result <- perform_DE_analysis(dds, contrast = comp,
                                     target_geneset = target_gene_list, 
                                     results.path = de.path)
  target_diff <- diff_result %>% filter(gene %in% target_gene_list)
  group <- paste0(comp[2], ".vs.", comp[3])
  target_diff$comp <- group
  target_diff_list[[group]] <- target_diff
  diff_results_list[[group]] <- diff_result
}
target_diff <- do.call(rbind, target_diff_list)
rownames(target_diff) <- NULL
write.csv(target_diff, paste0(de.path, target_gene_list[1], ".diff.csv"))

pdf(paste0(de.path, "check.quality.pdf"), width = 6, height = 6)
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
dev.off()


# enrichment analysis ----
enrich.path <- paste0(results.path, "/enrichment/")
dir.create(enrich.path, showWarnings = F)

source("/data02/hukaijie/EpiAllele/final_script/bulk_RNAseq/04_enrichment_function.R")

for (comp in comparisons) {
  group <- paste0(comp[2], ".vs.", comp[3])
  DEG <- read.csv(paste0(de.path, group, ".DEresult.csv"))
  up_genes <- DEG %>% 
    filter(group == "UP") %>% 
    pull(gene)
  down_genes <- DEG %>% 
    filter(group == "DOWN") %>% 
    pull(gene)
  
  if (length(up_genes) > 20) {
    GO_enrichment(gene=up_genes, name = paste0(group, ".up"), 
                  title = paste0(group, " (up-regulated)"), 
                  results.path = enrich.path, species = species)
    KEGG_enrichment(up_genes, name = paste0(group, ".up"), 
                    title = paste0(group, " (up-regulated)"), 
                    results.path = enrich.path, species = species)
  } else {
    print(paste0(group, " up-regulated genes are less than 20"))
  }
  
  if (length(down_genes) > 20) {
    GO_enrichment(down_genes, name = paste0(group, ".down"),
                  title = paste0(group, " (down-regulated)"), 
                  results.path = enrich.path, species = species)
    KEGG_enrichment(down_genes, name = paste0(group, ".down"), 
                    title = paste0(group, " (down-regulated)"), 
                    results.path = enrich.path, species = species)
  } else {
    print(paste0(group, " down-regulated genes are less than 20"))
  }
}





# expression heatmap ----
expr <- read.table(paste0(expr.path, "counts.txt"), header = T, check.names = F)
tpm <- read.table(paste0(expr.path, "tpm.txt"), header = T, check.names = F)
DE_results <- list.files(de.path, pattern = "DEresult.csv", full.names = T)
geneset <- c()

`%notin%` <- Negate(`%in%`)
for (DE_result in DE_results) {
  group <- sub(".*/(.*).DEresult.csv", "\\1", DE_result)
  DEG <- read.csv(DE_result)
  DEG <- DEG %>% 
    filter(group %in% c("UP", "DOWN"),
           gene %notin% c("MYH7_WT", "MYH7_MUTANT", "Myh6_C57", "Myh6_DBA"))
  
  DEG_expr <- tpm[rownames(tpm) %in% DEG$gene,]
  
  pdf(paste0(expr.path, group, "DEG.expr.heatmap.pdf"), width = 10, height = 10)
  p1 <- pheatmap(DEG_expr, 
                 scale = "row",
                 cluster_rows = T, cluster_cols = F,
                 show_rownames = F,
                 color = colorRampPalette(c("#0000EF","white","red"))(100),
                 border_color = F,
                 annotation_col = sampleTable[,"condition", drop = F],
                 annotation_row = DEG[,"group",drop = F],
                 annotation_colors = list(condition = colorset$condition,
                                          group = c("UP" = "red", "DOWN" = "blue")),
                 heatmap_legend_param = list(
                   title = "Expression level"
                 ))
  print(p1)
  dev.off()
  
  geneset <- c(geneset, DEG$gene)
}


## DEGs ----
DEG_expr <- tpm[rownames(tpm) %in% geneset,]
pdf(paste0(expr.path, "allDEG.expr.heatmap.pdf"), width = 10, height = 10)
p1 <- pheatmap(DEG_expr, 
               scale = "row",
               cluster_rows = T, cluster_cols = F,
               show_rownames = F,
               color = colorRampPalette(c("#0000EF","white","red"))(100),
               border_color = F,
               annotation_col = sampleTable[,"condition", drop = F],
               annotation_colors = list(condition = colorset$condition),
               heatmap_legend_param = list(
                 title = "Expression level"
               ))
print(p1)
dev.off()

## all gene ----
pdf(paste0(expr.path, "allgene.expr.heatmap.pdf"), width = 10, height = 10)
p1 <- pheatmap(tpm, 
               scale = "row",
               cluster_rows = T, cluster_cols = T,
               show_rownames = F,
               color = colorRampPalette(c("#0000EF","white","red"))(100),
               border_color = F,
               annotation_col = sampleTable[,"condition", drop = F],
               annotation_colors = list(condition = colorset$condition),
               heatmap_legend_param = list(
                 title = "Expression level"
               ))
print(p1)
dev.off()
