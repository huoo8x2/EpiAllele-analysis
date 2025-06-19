#!/usr/bin/env Rscript

## {{{ Install required libraries if required
for (p in c("tidyverse", "ggsci", "ggpubr", "stringr", "getopt", "viridis",
            "corrplot", "patchwork", "DESeq2", "GenomicRanges")) {
  if (!require(p, character.only = T)) {
    print(paste0("Please install ", p))
  }else{
    suppressMessages(library(p, quietly = T, character.only = T))
  }
}
## }}}


spec <- matrix(c(
  'results_path', 'r', 1, "character",
  'projPath', 'p', 1, "character",
  'spikein', 's', 1, "logical",
  'dup', 'd', 1, "logical"
), byrow=TRUE, ncol=4)

opt <- getopt(spec =spec)

results.path <- opt$results_path
if (! dir.exists(results.path)) {
  dir.create(results.path)
}
spikein=opt$spikein
dup=opt$dup

## Path to the project and histone list
projPath <- opt$projPath
sampleList <- list.dirs(projPath, full.names = FALSE, recursive = FALSE)
histList <- unique(sub("(.*)_rep.*", "\\1", sampleList))

# sequencing mapping summary ----
## Sequencing depth ----
## Collect the alignment results from the bowtie2 alignment summary files
alignResult = c()
for(sample in sampleList){
  alignRes = read.table(paste0(projPath, "/", sample, "/", sample, "_bowtie2_map.txt"), 
                        header = FALSE, fill = TRUE)
  alignRate = substr(alignRes$V1[6], 1, nchar(as.character(alignRes$V1[6]))-1)
  alignResult = data.frame(Histone = sub("(.*)_rep.*", "\\1", sample), 
                           Replicate = sub(".*_(rep.*)", "\\1", sample), 
                           SequencingDepth = alignRes$V1[1] %>% as.character %>% as.numeric, 
                           MappedFragNum = alignRes$V1[4] %>% as.character %>% as.numeric + alignRes$V1[5] %>% as.character %>% as.numeric, 
                           AlignmentRate = alignRate %>% as.numeric)  %>% rbind(alignResult, .)
}
alignResult$Histone = factor(alignResult$Histone, levels = histList)
alignResult %>% mutate(AlignmentRate = paste0(AlignmentRate, "%"))

## Spike-in alignment ----
if (spikein) {
  spikeAlign = c()
  for(sample in sampleList){
    spikeRes = read.table(paste0(projPath, "/", sample, "/", sample, "_bowtie2_spikeIn.txt"), 
                          header = FALSE, fill = TRUE)
    alignRate = substr(spikeRes$V1[6], 1, nchar(as.character(spikeRes$V1[6]))-1)
    spikeAlign = data.frame(Histone = sub("(.*)_rep.*", "\\1", sample), 
                            Replicate = sub(".*_(rep.*)", "\\1", sample), 
                            SequencingDepth = spikeRes$V1[1] %>% as.character %>% as.numeric, 
                            MappedFragNum_spikeIn = spikeRes$V1[4] %>% as.character %>% as.numeric + spikeRes$V1[5] %>% as.character %>% as.numeric, 
                            AlignmentRate_spikeIn = alignRate %>% as.numeric)  %>% rbind(spikeAlign, .)
  }
  spikeAlign$Histone = factor(spikeAlign$Histone, levels = histList)
  spikeAlign %>% mutate(AlignmentRate_spikeIn = paste0(AlignmentRate_spikeIn, "%"))
  write.csv(spikeAlign, paste0(results.path, "/spikein_alignment.csv"), row.names = FALSE)
  
  alignSummary = left_join(alignResult, spikeAlign, by = c("Histone", "Replicate", "SequencingDepth")) %>%
    mutate(AlignmentRate_hg38 = paste0(AlignmentRate, "%"), 
           AlignmentRate_spikeIn = paste0(AlignmentRate_spikeIn, "%"))
  write.csv(alignSummary, paste0(results.path, "/alignment.summary.csv"), row.names = FALSE)
  
} else{
  alignSummary = alignResult
  
  write.csv(alignSummary, paste0(results.path, "/alignment.summary.csv"), row.names = FALSE)
}


## Generate sequencing depth boxplot ----
seqdepth_fig = alignResult %>% ggplot(aes(x = Histone, y = SequencingDepth/1000000)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  ylab("Sequencing Depth per Million") +
  xlab("Group") + 
  ggtitle("Sequencing Depth")

mappednum_fig = alignResult %>% ggplot(aes(x = Histone, y = MappedFragNum/1000000)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  ylab("Mapped Fragments per Million") +
  xlab("Group") + 
  ggtitle("Alignable Fragment")

alignrate_fig = alignResult %>% ggplot(aes(x = Histone, y = AlignmentRate)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
  ylab("% of Mapped Fragments") +
  xlab("Group") + 
  ggtitle("Alignment Rate")

if (spikein) {
  spikein_fig = spikeAlign %>% ggplot(aes(x = Histone, y = MappedFragNum_spikeIn)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
    ylab("Spike-in Mapped Fragments Num") +
    xlab("Group") + 
    ggtitle("Alignable Fragment (E.coli)")
  p <- wrap_plots(seqdepth_fig, mappednum_fig, alignrate_fig, spikein_fig, ncol = 2, nrow=2)
} else {
  p <- wrap_plots(seqdepth_fig, mappednum_fig, alignrate_fig, ncol = 2, nrow=2)
}
pdf(paste0(results.path, "/sequencing_summary.pdf"), width = 10, height = 8)
print(p)
dev.off()

# Assess duplicated rates ----
if (dup) {
  ## Summarize the duplication information from the picard summary outputs ----
  dupResult = c()
  for(sample in sampleList){
    dupRes = read.table(paste0(projPath, "/", sample, "/", sample, ".picard.rmDup.txt"), header = TRUE, fill = TRUE)
    
    dupResult = data.frame(Histone = sub("(.*)_rep.*", "\\1", sample), 
                           Replicate = sub(".*_(rep.*)", "\\1", sample),  
                           MappedFragNum = dupRes$READ_PAIRS_EXAMINED[1] %>% 
                             as.character %>% as.numeric,
                           DuplicationRate = dupRes$PERCENT_DUPLICATION[1] %>% 
                             as.character %>% as.numeric * 100,
                           EstimatedLibrarySize = dupRes$ESTIMATED_LIBRARY_SIZE[1] %>% 
                             as.character %>% as.numeric) %>% 
      mutate(UniqueFragNum = MappedFragNum * (1-DuplicationRate/100))  %>% rbind(dupResult, .)
  }
  dupResult$Histone = factor(dupResult$Histone, levels = histList)
  alignDupSummary = left_join(alignSummary, dupResult, by = c("Histone", "Replicate", "MappedFragNum")) %>%
    mutate(DuplicationRate = paste0(DuplicationRate, "%"))
  alignDupSummary
  
  write.csv(alignDupSummary, paste0(results.path, "/alignDup.summary.csv"), row.names = FALSE)
  
  
  ## generate boxplot figure for the  duplication rate ----
  duprate_fig = dupResult %>% ggplot(aes(x = Histone, y = DuplicationRate)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
    ylab("Duplication Rate (*100%)") +
    xlab("Group") 
  
  
  libsize_fig = dupResult %>% ggplot(aes(x = Histone, y = EstimatedLibrarySize)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
    ylab("Estimated Library Size") +
    xlab("Group") 
  
  uniFrag_fig = dupResult %>% ggplot(aes(x = Histone, y = UniqueFragNum)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
    ylab("# of Unique Fragments") +
    xlab("Group") 
  
  p <- wrap_plots(duprate_fig, libsize_fig, uniFrag_fig, ncol = 2, nrow=2)
  pdf(paste0(results.path, "/duplication_summary.pdf"), width = 10, height = 8)
  print(p)
  dev.off()
}


# Assess mapped fragment size distribution ----
## Collect the fragment size information ----
fragLen = c()
for(sample in sampleList){
  fragLen = read.table(paste0(projPath, "/", sample, "/", sample, "_fragmentLen.txt"), 
                        header = FALSE, fill = TRUE) %>% 
    mutate(fragLen = V1 %>% as.numeric, 
           fragCount = V2 %>% as.numeric, 
           Weight = as.numeric(V2)/sum(as.numeric(V2)), 
           Histone = sub("(.*)_rep.*", "\\1", sample), 
           Replicate = sub(".*_(rep.*)", "\\1", sample), 
           sampleInfo = sample) %>%
    rbind(fragLen, .)
}
fragLen$sampleInfo = factor(fragLen$sampleInfo, levels = sampleList)
fragLen$Histone = factor(fragLen$Histone, levels = histList)

write.csv(fragLen, paste0(results.path, "/fragLen.csv"))

## Generate the fragment size density plot (violin plot) ----
vlnplot = fragLen %>% ggplot(aes(x = sampleInfo, y = fragLen, weight = Weight)) +
  geom_violin(fill="lightblue") +
  scale_y_continuous(breaks = seq(0, 800, 50)) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
  ylab("Fragment Length") +
  xlab("Group")

Den_plot = fragLen %>% ggplot(aes(x = fragLen, y = fragCount, color = Histone, 
                                  group = sampleInfo, linetype = Replicate)) +
  geom_line(size = 1) +
  scale_color_npg() +
  theme_bw() +
  xlab("Fragment Length") +
  ylab("Count")

p <- wrap_plots(vlnplot, Den_plot, ncol = 1, nrow=2)
pdf(paste0(results.path, "/fragLen_stat.pdf"), width = 8, height = 8)
print(p)
dev.off()

# Assess replicate reproducibility ----
reprod = c()
fragCount = NULL
for(sample in sampleList){
  if(is.null(fragCount)){
    fragCount = read.table(paste0(projPath, "/", sample, "/", sample, "_bowtie2.fragmentsCount.bin500.bed"), 
                           header = FALSE) 
    colnames(fragCount) = c("chrom", "bin", sample)
  }else{
    fragCountTmp = read.table(paste0(projPath, "/", sample, "/", sample, "_bowtie2.fragmentsCount.bin500.bed"), 
                              header = FALSE) 
    colnames(fragCountTmp) = c("chrom", "bin", sample)
    fragCount = full_join(fragCount, fragCountTmp, by = c("chrom", "bin"))
  }
}
M = cor(fragCount %>% select(-c("chrom", "bin")) %>% log2(), use = "complete.obs")

pdf(paste0(results.path, "/fragCounts_corplot.pdf"), width = 12, height = 9)
corrplot(M, method = "color", outline = T, addgrid.col = "darkgray", order="hclust", 
         addrect = 3, rect.col = "black", rect.lwd = 3,cl.pos = "b", tl.col = "indianred4", 
         tl.cex = 1, cl.cex = 1, 
         col = colorRampPalette(c("midnightblue","white","darkred"))(100))
dev.off()

# scaling factor stat ----
scaleFactor = c()
if (spikein) {
  for(sample in sampleList){
    spikeDepth = read.table(paste0(projPath, "/", sample, "/", sample, "_bowtie2_spikeIn.seqDepth"), 
                            header = FALSE, fill = TRUE)$V1[1]
    seqDepth = read.table(paste0(projPath, "/", sample, "/", sample, "_bowtie2_mapped.seqDepth"), 
                          header = FALSE, fill = TRUE)$V1[1]
    scaleFactor = data.frame(scaleFactor_spikein = 10000/spikeDepth, 
                             scaleFactor_depth = 10000000/seqDepth, 
                             Histone = sub("(.*)_rep.*", "\\1", sample), 
                             Replicate = sub(".*_(rep.*)", "\\1", sample)) %>% 
      rbind(scaleFactor, .)    
  }
} else{
  for(sample in sampleList){
    seqDepth = read.table(paste0(projPath, "/", sample, "/", sample, "_bowtie2_mapped.seqDepth"), 
                          header = FALSE, fill = TRUE)$V1[1]
    scaleFactor = data.frame(scaleFactor_depth = 10000000/seqDepth, 
                             Histone = sub("(.*)_rep.*", "\\1", sample), 
                             Replicate = sub(".*_(rep.*)", "\\1", sample)) %>% 
      rbind(scaleFactor, .)
  }
}

scaleFactor$Histone = factor(scaleFactor$Histone, levels = histList)

## Generate sequencing depth boxplot ----
if (spikein) {
  spikein_scale_fig = scaleFactor %>% ggplot(aes(x = Histone, y = scaleFactor_spikein)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
    ylab("Spike-in Scalling Factor") +
    xlab("Group")   

  seqdepth_scale_fig = scaleFactor %>% ggplot(aes(x = Histone, y = scaleFactor_depth)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
    ylab("SeqDepth Scalling Factor") +
    xlab("Group") 
  
  normDepth = inner_join(scaleFactor, alignResult, by = c("Histone", "Replicate")) %>% 
    mutate(normDepth_spikein = as.numeric(MappedFragNum * scaleFactor_spikein),
           normDepth_seqdepth = as.numeric(MappedFragNum * scaleFactor_depth))
  
  normDepth_spikein_fig = normDepth %>% ggplot(aes(x = Histone, y = normDepth_spikein)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
    ylab("Normalization Fragment Count (spikein scaled)") +
    xlab("Group") 
  
  normDepth_seqdepth_fig = normDepth %>% ggplot(aes(x = Histone, y = normDepth_seqdepth)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
    scale_y_continuous(labels = scales::comma) + 
    ylab("Normalization Fragment Count (seqDepth scaled)") +
    xlab("Group") 
  
  p <- wrap_plots(spikein_scale_fig, seqdepth_scale_fig,
                  normDepth_spikein_fig, normDepth_seqdepth_fig,
                  ncol = 2, nrow=2)
  
    write.csv(normDepth, paste0(results.path, "/scaling_factor.csv"))
    pdf(paste0(results.path, "/scaling_factor.pdf"), width = 10, height = 8)
    print(p)
    dev.off()

} else{
  scale_factor_fig = scaleFactor %>% ggplot(aes(x = Histone, y = scaleFactor_depth)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
    ylab("SeqDepth Scalling Factor") +
    xlab("Group") 
  
  normDepth = inner_join(scaleFactor, alignResult, by = c("Histone", "Replicate")) %>% 
    mutate(normDepth = as.numeric(MappedFragNum * scaleFactor_depth))
  
  normDepth_fig = normDepth %>% ggplot(aes(x = Histone, y = normDepth)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) +
    ylab("Normalization Fragment Count (seqDepth scaled)") +
    xlab("Group") 
  
  p <- wrap_plots(scale_factor_fig, normDepth_fig, ncol = 2, nrow=1)
  
  write.csv(normDepth, paste0(results.path, "/scaling_factor.csv"))
  pdf(paste0(results.path, "/scaling_factor.pdf"), width = 10, height = 4)
  print(p)
  dev.off()
}


