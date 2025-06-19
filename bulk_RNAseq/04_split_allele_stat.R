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
  'species', 's', 2, "character"
), byrow=TRUE, ncol=4)
opt <- getopt(spec =spec)
result.path <- opt$result_path
species <- opt$species

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

width=length(allele_fraction$dataset)*0.5
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
if (species == "mouse") {
  allele1="C57"
  allele2="DBA"
  target_gene="Myh6"
} else if (species == "human") {
  allele1="MUTANT"
  allele2="WT"
  target_gene="MYH7"
} else if (species == "human_v2"){
  allele1="UNTARGET"
  allele2="TARGET"
  target_gene="MYH7"
} else {
  stop("species should be mouse or human")
}

sample_parent_group$sample <- rownames(sample_parent_group)
sample_parent_group <- sample_parent_group %>%
  mutate(split_sample = strsplit(as.character(sample), "-"),
         Group = sapply(split_sample, `[`, 1),
         Replicate = paste0("rep", sapply(split_sample, `[`, 2)),
         ratio = get(allele1) / get(allele2),
         sum = get(allele1) + get(allele2)) %>%
  dplyr::select(-split_sample)

allele1_fig = sample_parent_group %>% ggplot(aes(x = Group, y = get(allele1))) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
  ylab(paste0(allele1, " allele counts")) +
  xlab("Group")
allele2_fig = sample_parent_group %>% ggplot(aes(x = Group, y = get(allele2))) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
  ylab(paste0(allele2, " allele counts")) +
  xlab("Group")
sum_fig = sample_parent_group %>% ggplot(aes(x = Group, y = sum)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
  ylab(paste0(target_gene, " region counts")) +
  xlab("Group")
alleleRatio_fig = sample_parent_group %>% ggplot(aes(x = Group, y = ratio)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
  ylab(paste0(allele1, " / ", allele2, " ratio")) +
  xlab("Group")

p <- wrap_plots(list(allele1_fig, allele2_fig, sum_fig, alleleRatio_fig), ncol = 2)
pdf(paste0(result.path, "/allelestat_boxplot.pdf"), width = 10, height = 10)
print(p)
dev.off()





