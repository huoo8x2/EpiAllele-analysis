#!/usr/bin/env Rscript

## {{{ Install required libraries if required
for (p in c("tidyverse", "ggsci", "ggpubr", "stringr", "getopt",
            "viridis", "patchwork")) {
  if (!require(p, character.only = T)) {
    print(paste0("Please install ", p))
  }else{
    suppressMessages(library(p, quietly = T, character.only = T))
  }
}
## }}}

spec <- matrix(c(
  'result_path', 'r', 2, "character",
  'spikein', 's', 1, "logical"
), byrow=TRUE, ncol=4)
opt <- getopt(spec =spec)
result.path <- opt$result_path
spikein <- opt$spikein

file_list <- list.files(result.path, pattern = "group.info.csv")
sample_parent_group_list <- list()
for (file in file_list) {
  sample <- str_remove(file, ".filtered.reads.group.info.csv")
  info <- read.csv(paste0(result.path, "/", file), header = T, check.names = F)
  sample_parent_group_list[[sample]] <- table(info$group)
}
sample_parent_group <- bind_rows(sample_parent_group_list, .id = "sample") %>%
  column_to_rownames("sample") %>%
  mutate_all(~ifelse(is.na(.), 0, .))
head(sample_parent_group)
write.csv(sample_parent_group, paste0(result.path, "/allele_stat.csv"))

## plot a barplot to show the fraction of each allele in each sample
allele_fraction <- sample_parent_group %>%
  as.data.frame() %>%
  rownames_to_column("dataset") %>%
  gather(allele, value, -dataset) %>% 
  filter(allele != "unknown") %>% 
  group_by(dataset) %>% 
  mutate(fraction = value/sum(value))

width=length(allele_fraction$dataset)*1
pdf(paste0(result.path, "/allele_stat.pdf"), width = width, height = 5)
p <- ggplot(allele_fraction, aes(dataset, fraction, fill=allele))+
  geom_bar(position = "stack", stat = "identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))
print(p)
p <- ggplot(allele_fraction, aes(dataset, value, fill=allele))+
  geom_bar(position = "stack", stat = "identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))
print(p)
dev.off()

## boxplot for allele ratio, counts and normalized counts
sample_parent_group$sample <- rownames(sample_parent_group)
sample_parent_group <- sample_parent_group %>%
  mutate(Histone = sub("(.*)_rep.*", "\\1", sample), 
         Replicate = sub(".*_(rep.*)", "\\1", sample), 
         ratio = `C57`/`DBA`, sum=`C57`+`DBA`) 
scaling_factor <- read.csv(str_replace(result.path, "/allele/stat", "/stat_plot/scaling_factor.csv"), 
                           header = T, row.names = 1)


c57_fig = sample_parent_group %>% ggplot(aes(x = Histone, y = C57)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
  ylab("C57 allele counts") +
  xlab("Group")
dba_fig = sample_parent_group %>% ggplot(aes(x = Histone, y = DBA)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
  ylab("DBA allele counts") +
  xlab("Group")
sum_fig = sample_parent_group %>% ggplot(aes(x = Histone, y = sum)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
  ylab("Myh6 region counts") +
  xlab("Group")
alleleRatio_fig = sample_parent_group %>% ggplot(aes(x = Histone, y = ratio)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
  ylab("Allele Ratio") +
  xlab("Group")

if (spikein) {
  normDepth = inner_join(scaling_factor, sample_parent_group, by = c("Histone", "Replicate")) %>% 
    mutate(normDepth_spikein = MappedFragNum * scaleFactor_spikein,
           normDepth_seqdepth = MappedFragNum * scaleFactor_depth)
  
  normDepth_spikein_fig = normDepth %>% ggplot(aes(x = Histone, y = normDepth_spikein)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
    ylab("Target Gene Normalization Fragment Count (spikein scaled)") +
    xlab("Group") 
  
  normDepth_seqdepth_fig = normDepth %>% ggplot(aes(x = Histone, y = normDepth_seqdepth)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
    scale_y_continuous(labels = scales::comma) +
    ylab("Target Gene Normalization Fragment Count (seqDepth scaled)") +
    xlab("Group") 
  
  p <- wrap_plots(list(c57_fig, dba_fig, sum_fig, alleleRatio_fig, 
                       normDepth_spikein_fig, normDepth_seqdepth_fig), ncol = 2)
} else {
  normDepth = inner_join(scaling_factor, sample_parent_group, by = c("Histone", "Replicate")) %>% 
    mutate(normDepth = as.numeric(MappedFragNum * scaleFactor_depth))
  
  normDepth_fig = normDepth %>% ggplot(aes(x = Histone, y = normDepth)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
    scale_y_continuous(labels = scales::comma) +
    ylab("Target Gene Normalization Fragment Count (seqDepth scaled)") +
    xlab("Group") 
  
  p <- wrap_plots(list(c57_fig, dba_fig, sum_fig, alleleRatio_fig, normDepth_fig), ncol = 2)
}

pdf(paste0(result.path, "/allelestat_boxplot.pdf"), width = 10, height = 10)
print(p)
dev.off()





