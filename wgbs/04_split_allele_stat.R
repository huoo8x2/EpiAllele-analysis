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

width=length(allele_fraction$dataset)*0.5
pdf(paste0(result.path, "/allele_stat.pdf"), width = width, height = 5)
p <- ggplot(allele_fraction, aes(dataset, fraction, fill=allele))+
  geom_bar(position = "stack", stat = "identity")+
  theme_bw()+
  theme(axis.text.x = element_text(vjust=1, hjust=1, angle = 45))
print(p)
p <- ggplot(allele_fraction, aes(dataset, value, fill=allele))+
  geom_bar(position = "stack", stat = "identity")+
  theme_bw()+
  theme(axis.text.x = element_text(vjust=1, hjust=1, angle = 45))
print(p)
dev.off()